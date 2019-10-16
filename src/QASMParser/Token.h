//
// Created by Lukas Burgholzer on 22.10.19.
//

#ifndef INTERMEDIATEREPRESENTATION_TOKEN_H
#define INTERMEDIATEREPRESENTATION_TOKEN_H

#include <map>

namespace qasm {

	struct Token {

		enum class Kind {
			include,
			none,
			identifier,
			number,
			plus,
			semicolon,
			eof,
			lpar,
			rpar,
			lbrack,
			rbrack,
			lbrace,
			rbrace,
			comma,
			minus,
			times,
			nninteger,
			real,
			qreg,
			creg,
			ugate,
			cxgate,
			gate,
			pi,
			measure,
			openqasm,
			probabilities,
			sin,
			cos,
			tan,
			exp,
			ln,
			sqrt,
			div,
			power,
			string,
			gt,
			barrier,
			opaque,
			_if,
			eq,
			reset,
			snapshot
		};
		
		Kind kind = Kind::none;
		int line = 0;
		int col = 0;
		int val = 0;
		double valReal = 0.0;
		std::string str;

		Token() = default;
		Token(Kind kind, int line, int col) : kind(kind), line(line), col(col) { }
	};

	static std::map<Token::Kind, std::string> KindNames{
			{ Token::Kind::none,          "none" },
			{ Token::Kind::include,       "include" },
			{ Token::Kind::identifier,    "<identifier>" },
			{ Token::Kind::number,        "<number>" },
			{ Token::Kind::plus,          "+" },
			{ Token::Kind::semicolon,     ";" },
			{ Token::Kind::eof,           "EOF" },
			{ Token::Kind::lpar,          "(" },
			{ Token::Kind::rpar,          ")" },
			{ Token::Kind::lbrack,        "[" },
			{ Token::Kind::rbrack,        "]" },
			{ Token::Kind::lbrace,        "{" },
			{ Token::Kind::rbrace,        "}" },
			{ Token::Kind::comma,         "," },
			{ Token::Kind::minus,         "-" },
			{ Token::Kind::times,         "*" },
			{ Token::Kind::nninteger,     "<nninteger>" },
			{ Token::Kind::real,          "<real>" },
			{ Token::Kind::qreg,          "qreg" },
			{ Token::Kind::creg,          "creg" },
			{ Token::Kind::ugate,         "U" },
			{ Token::Kind::cxgate,        "CX" },
			{ Token::Kind::gate,          "gate" },
			{ Token::Kind::pi,            "pi" },
			{ Token::Kind::measure,       "measure" },
			{ Token::Kind::openqasm,      "openqasm" },
			{ Token::Kind::probabilities, "probabilities" },
			{ Token::Kind::opaque,        "opaque" },
			{ Token::Kind::sin,           "sin" },
			{ Token::Kind::cos,           "cos" },
			{ Token::Kind::tan,           "tan" },
			{ Token::Kind::exp,           "exp" },
			{ Token::Kind::ln,            "ln" },
			{ Token::Kind::sqrt,          "sqrt" },
			{ Token::Kind::div,           "/" },
			{ Token::Kind::power,         "^" },
			{ Token::Kind::string,        "string" },
			{ Token::Kind::gt,            ">" },
			{ Token::Kind::barrier,       "barrier" },
			{ Token::Kind::_if,           "if" },
			{ Token::Kind::eq,            "==" },
			{ Token::Kind::reset,         "reset" }
	};

}

#endif //INTERMEDIATEREPRESENTATION_TOKEN_H