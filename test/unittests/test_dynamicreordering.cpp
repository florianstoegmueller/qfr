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
	e = qc->buildFunctionality(dd, dd::Sifting);
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
	e = qc->buildFunctionality(dd, dd::Sifting);
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
	e = qc->buildFunctionality(dd, dd::Sifting);
	EXPECT_EQ(dd->size(e), 32);
}
