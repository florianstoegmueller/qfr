/*
 * This file is part of IIC-JKU QFR library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum/ for more information.
 */

#include "gtest/gtest.h"
#include "QuantumComputation.hpp"
#include <sstream>

class DynamicReorderingTest : public testing::TestWithParam<std::string> {

protected:
	std::unique_ptr<qc::QuantumComputation> qc;
	std::array<short, qc::MAX_QUBITS>       line{};
	dd::Edge                                e{}, in{};
	std::unique_ptr<dd::Package>            dd;
	qc::permutationMap                      varMap{};

	std::string circuit_dir = "./circuits/";

	void SetUp() override {
		dd = std::make_unique<dd::Package>();
		line.fill(qc::LINE_DEFAULT);
		qc = std::make_unique<qc::QuantumComputation>();
	}

	void TearDown() override {
		qc->reset();
	}

};

class DynamicReorderingTestVisualisation : public testing::TestWithParam<std::string> {

protected:
	std::unique_ptr<qc::QuantumComputation> qc;
	std::array<short, qc::MAX_QUBITS>       line{};
	dd::Edge                                e{}, in{};
	std::unique_ptr<dd::Package>            dd;
	qc::permutationMap varMapNone{}, varMapSifting{};
	dd::Edge none{}, sifting{};

	std::string circuit_dir = "./circuits/";
	std::string output_dir = "./output/";

	void SetUp() override {
		dd = std::make_unique<dd::Package>();
		line.fill(qc::LINE_DEFAULT);
		qc = std::make_unique<qc::QuantumComputation>();
		qc->import(circuit_dir + GetParam() + ".qasm");
	}

	void TearDown() override {
		qc->reset();
	}

};

TEST_F(DynamicReorderingTest, cx_exchange) {
	auto cx = qc::StandardOperation(2,qc::Control(1),0,qc::X);
	auto cx_rev = qc::StandardOperation(2,qc::Control(0),1,qc::X);

	auto cx_dd = cx.getDD(dd, line);
	dd->incRef(cx_dd);
	dd->printDD(cx_dd, 64);
	dd->printUniqueTable(2);

	auto cx_exg = dd->exchangeBaseCase(cx_dd, 0, 1);
	dd->printDD(cx_exg, 64);
	dd->printUniqueTable(2);

	auto cx_rev_dd = cx_rev.getDD(dd, line);
	dd->incRef(cx_rev_dd);
	dd->printDD(cx_rev_dd, 64);
	dd->printUniqueTable(2);

	EXPECT_TRUE(dd->equals(cx_exg, cx_rev_dd));
}

TEST_F(DynamicReorderingTest, cx_exchange_unique_table) {
	auto cx = qc::StandardOperation(2,qc::Control(1),0,qc::X);
	auto cx_rev = qc::StandardOperation(2,qc::Control(0),1,qc::X);

	auto cx_rev_dd = cx_rev.getDD(dd, line);
	dd->incRef(cx_rev_dd);
	dd->printDD(cx_rev_dd, 64);
	dd->printUniqueTable(2);

	auto cx_dd = cx.getDD(dd, line);
	dd->incRef(cx_dd);
	dd->printDD(cx_dd, 64);
	dd->printUniqueTable(2);

	auto cx_exg = dd->exchangeBaseCase(cx_dd, 0, 1);
	dd->printDD(cx_exg, 64);
	dd->printUniqueTable(2);

	EXPECT_TRUE(dd->equals(cx_exg, cx_rev_dd));
}

TEST_F(DynamicReorderingTest, toffoli_sifting) {
	std::stringstream ss{};
	ss
	<< ".numvars 3\n"
	<< ".variables a b c\n"
	<< ".begin\n"
	<< "t3 a b c\n"
	<< ".end\n";
	qc->import(ss, qc::Real);
	std::tie(e, varMap) = qc->buildFunctionality(dd, dd::Sifting);
	EXPECT_EQ(dd->size(e), 6);
}

TEST_F(DynamicReorderingTest, mct_sifting_small) {
	std::stringstream ss{};
	ss
	<< ".numvars 4\n"
	<< ".variables a b c d\n"
	<< ".begin\n"
	<< "t4 a b c d\n"
	<< ".end\n";
	qc->import(ss, qc::Real);
	std::tie(e, varMap) = qc->buildFunctionality(dd, dd::Sifting);
	EXPECT_EQ(dd->size(e), 8);
}

TEST_F(DynamicReorderingTest, mct_sifting_large) {
	// best case for mct gate is that the target is on the least significand qubit q0 (2*n nodes incl. terminal)
	// worst case for mct gate is that the target is on the most significand qubit qn-1
	// sifting should be able to turn the worst into the best case.
	std::stringstream ss{};
	ss
			<< ".numvars 16\n"
			<< ".variables a b c d e f g h i j k l m n o p\n"
			<< ".begin\n"
			<< "t16 a b c d e f g h i j k l m n o p\n"
			<< ".end\n";
	qc->import(ss, qc::Real);
	std::tie(e, varMap) = qc->buildFunctionality(dd, dd::Sifting);
	EXPECT_EQ(dd->size(e), 32);
}

INSTANTIATE_TEST_SUITE_P(SomeCircuits, DynamicReorderingTestVisualisation, testing::Values("bell", "grover", "test2", "test3", "test4"),
		[](const testing::TestParamInfo<DynamicReorderingTestVisualisation::ParamType>& info) {
			auto s = info.param;
			std::replace( s.begin(), s.end(), '.', '_');
			return s;
		});

TEST_P(DynamicReorderingTestVisualisation, simulationSize) {
	in = dd->makeZeroState(qc->getNqubits());
	std::tie(none, varMapNone) = qc->simulate(in, dd, dd::None);
	std::stringstream ss{};
	ss << output_dir << GetParam() << "_sim_none.dot";
	dd->export2Dot(none, ss.str(), true);
	auto sizeNone = dd->size(none);

	qc->reset();
	qc->import(circuit_dir + GetParam() + ".qasm");

	in = dd->makeZeroState(qc->getNqubits());
	std::tie(sifting, varMapSifting) = qc->simulate(in, dd, dd::Sifting);
	std::stringstream ss2{};
	ss2 << output_dir << GetParam() << "_sim_sifting.dot";
	dd->export2Dot(sifting, ss2.str(), true);
	auto sizeSifting = dd->size(sifting);
	for (const auto& var: varMapSifting) {
		if (var.first >= qc->getNqubits()) break;
		std::cout << var.first << ": " << var.second << std::endl;
	}
	std::cout << "sifting size: " << sizeSifting << " vs. orginal size: " << sizeNone << std::endl;
	EXPECT_LE(sizeSifting, sizeNone);
}

TEST_P(DynamicReorderingTestVisualisation, constructionSize) {
	std::tie(none, varMapNone) = qc->buildFunctionality(dd, dd::None);
	std::stringstream ss{};
	ss << output_dir << GetParam() << "_matrix_none.dot";
	dd->export2Dot(none, ss.str(), false);
	auto sizeNone = dd->size(none);

	qc->reset();
	qc->import(circuit_dir + GetParam() + ".qasm");

	std::tie(sifting, varMapSifting) = qc->buildFunctionality(dd, dd::Sifting);
	std::stringstream ss2{};
	ss2 << output_dir << GetParam() << "_matrix_sifting.dot";
	dd->export2Dot(sifting, ss2.str(), false);
	auto sizeSifting = dd->size(sifting);
	for (const auto& var: varMapSifting) {
		if (var.first >= qc->getNqubits()) break;
		std::cout << var.first << ": " << var.second << std::endl;
	}
	std::cout << "sifting size: " << sizeSifting << " vs. orginal size: " << sizeNone << std::endl;
	EXPECT_LE(sizeSifting, sizeNone);
}
