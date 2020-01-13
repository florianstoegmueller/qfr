/*
 * This file is part of IIC-JKU QFR library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum/ for more information.
 */

#include "NonUnitaryOperation.hpp"

namespace qc {
    // Measurement constructor
	NonUnitaryOperation::NonUnitaryOperation(const unsigned short nq, const std::vector<unsigned short>& qubitRegister, const std::vector<unsigned short>& classicalRegister) 
		: NonUnitaryOperation(nq, classicalRegister, Measure) {
		//targets = qubitRegister;
		//targets = classicalRegister;
		assert(qubitRegister.size() == classicalRegister.size());
		// i-th qubit to be measured shall be measured into i-th classical register
		for (const auto& qubit: qubitRegister)
			controls.emplace_back(qubit);
	}

	// Snapshot constructor
	NonUnitaryOperation::NonUnitaryOperation(const unsigned short nq, const std::vector<unsigned short>& qubitRegister, int n) 
		: NonUnitaryOperation(nq, qubitRegister, Snapshot) {
		parameter[0] = n;
	}	

	// Snapshot constructor
	NonUnitaryOperation::NonUnitaryOperation(const unsigned short nq, const std::vector<unsigned short>& qubitRegister, Op op) {
		this->op = op;
		nqubits  = nq;
		targets  = qubitRegister;
		strcpy(name, opNames[op].c_str());
	}

    std::ostream& NonUnitaryOperation::print(std::ostream& os) const {
	    std::array<short, MAX_QUBITS> line{};
	    line.fill(LINE_DEFAULT);

		switch (op) {
			case Measure: 
				os << "Meas\t";
				for (int q = 0; q < controls.size(); ++q) {
					line[controls[q].qubit] = targets[q];
				}
				for (int i = 0; i < nqubits; ++i) {
					if (line[i] >= 0) {
						os << "\033[34m" << line[i] << "\t" << "\033[0m";
					} else {
						os << "|\t";
					}
				}
				break;
			case Reset: 
				os << "Rst \t";
				setLine(line);
				for (int i = 0; i < nqubits; ++i) {
					if (line[i] == LINE_TARGET) {
						os << "\033[31m" << "r\t" << "\033[0m";
					} else {
						os << "|\t";
					}
				}
				break;
			case Snapshot:
				os << "Snap\t";
				setLine(line);
				for (int i = 0; i < nqubits; ++i) {
					if (line[i] == LINE_TARGET) {
						os << "\033[33m" << "s\t" << "\033[0m";
					} else {
						os << "|\t";
					}
				}
				os << "\tp: " << targets.size() << " " << parameter[1];
				break;
			case ShowProbabilities: 
				os << "Show probabilities";
				break;
			case Barrier: 
				os << "Barrier \t";
				setLine(line);
				for (int i = 0; i < nqubits; ++i) {
					if (line[i] == LINE_TARGET) {
						os << "\033[31m" << "r\t" << "\033[0m";
					} else {
						os << "|\t";
					}
				}
				break;
		}
		return os;
	}

	void NonUnitaryOperation::dumpOpenQASM(std::ofstream& of, const regnames_t& qreg, const regnames_t& creg) const {
		switch (op) {
			case Measure: 
				if(isWholeQubitRegister(qreg, controls[0].qubit, controls.back().qubit) && 
				   isWholeQubitRegister(qreg, targets[0],        targets.back())) {
					of << "measure " << qreg[controls[0].qubit].first << " -> " << creg[targets[0]].first << ";" << std::endl;
				} else {
					for (int q = 0; q < controls.size(); ++q) {
						of << "measure " << qreg[controls[q].qubit].second << " -> " << creg[targets[q]].second << ";" << std::endl;
					}
				}
				break;
			case Reset: 
				if(isWholeQubitRegister(qreg, targets[0], targets.back())) {
					of << "reset " << qreg[targets[0]].first << ";" << std::endl;
				} else {
					for (auto target: targets) {
						of << "reset " << qreg[target].second << ";" << std::endl;
					}
				}
				break;
			case Snapshot: 
				if(targets.size() > 0) {
					of << "snapshot(" << parameter[0] << ") ";
					
					for (int q = 0; q < targets.size(); ++q) {
						if(q > 0) {
							of << ", ";
						}
						of << qreg[targets[q]].second;
					}
					of << ";" << std::endl;
				}
				break;
			case ShowProbabilities: 
				of << "show_probabilities;" << std::endl;
				break;
			case Barrier: 
				if(isWholeQubitRegister(qreg, targets[0],        targets.back())) {
					of << "barrier " << qreg[targets[0]].first << ";" << std::endl;
				} else {
					for (auto target: targets) {
						of << "barrier " << qreg[target].first << ";" << std::endl;
					}
				}
				break;
		}
	}
}