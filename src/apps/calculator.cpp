#include "chalk/number.h"
#include "util/lexer.h"
#include <cassert>
#include <cctype>
#include <fmt/format.h>
#include <iterator>
#include <map>
#include <optional>
#include <readline/history.h>
#include <readline/readline.h>
#include <string>
using namespace chalk;

static char historyFile[] = ".chalk_history";

// trim whitespace from both ends of a string
std::string trimWhite(std::string_view s)
{
	size_t i = 0, j = s.size();
	while (i < s.size() && std::isspace(s[i]))
		++i;
	while (j > i && std::isspace(s[j - 1]))
		--j;
	return std::string{s.substr(i, j - i)};
}

/**
 * Get one line of user input from the terminal using readline library
 *   - trims whitespace from both ends
 *   - returns std::nullpt on EOF (i.e. CRTL-D), which is different from ""
 */
std::optional<std::string> getInput()
{
	char *s = readline(">> ");
	if (s == nullptr)
		return std::nullopt;
	auto line = trimWhite(s);
	if (!line.empty())
		add_history(s); // add non-trimmed to history, to preserve indentation
	free(s);
	return line;
}

template <typename T> class Parser
{
  public:
	util::Lexer lexer;

	// NOTE: the 'std::less<>' enables heterogeneous lookup
	using SymbolTable = std::map<std::string, T, std::less<>>;

	using Tok = util::Tok;
	using Token = util::Token;

  private:
	SymbolTable const *symbolTable = nullptr;
	T parsePrimary()
	{
		if (lexer.tryMatch(Tok::OpenParen))
		{
			auto r = parseSum();
			lexer.match(Tok::CloseParen);
			return r;
		}
		else if (auto tok = lexer.tryMatch(Tok::Int); tok)
		{
			return T(Integer(tok->value));
		}
		else if (auto tok = lexer.tryMatch(Tok::Float); tok)
		{
			return T(FloatingOctuple(tok->value));
		}
		else if (auto tok = lexer.tryMatch(Tok::Ident); tok)
		{
			if (lexer.tryMatch(Tok::OpenParen)) // function call
			{
				if (tok->value == "sqrt")
				{
					auto r = sqrt(parseSum());
					lexer.match(Tok::CloseParen);
					return r;
				}
				else
					throw util::ParseError(
					    fmt::format("unknown function '{}'", tok->value));
			}
			else // variable / constant
			{
				if (symbolTable != nullptr)
					if (auto it = symbolTable->find(tok->value);
					    it != symbolTable->end())
						return it->second;

				throw util::ParseError(
				    fmt::format("unknown variable '{}'", tok->value));
			}
		}
		else
			throw util::ParseError(
			    fmt::format("unexpected token '{}'", lexer->value));
	}

	T parsePower()
	{
		bool neg = lexer.tryMatch(Tok::Sub);
		T r = parsePrimary();
		if (lexer.tryMatch(Tok::Pow))
		{
			if (lexer.tryMatch(Tok::Sub))
				r = pow(r, -util::parse_int(lexer.match(Tok::Int).value));
			else
				r = pow(r, util::parse_int(lexer.match(Tok::Int).value));
		}
		return neg ? -r : r;
	}

	T parseProduct()
	{
		T r = parsePower();
		while (true)
			if (lexer.tryMatch(Tok::Mul))
				r = r * parsePower();
			else if (lexer.tryMatch(Tok::Div))
				r = r / parsePower();
			else
				return r;
	}

	T parseSum()
	{
		T r = parseProduct();
		while (true)
			if (lexer.tryMatch(Tok::Add))
				r = r + parseProduct();
			else if (lexer.tryMatch(Tok::Sub))
				r = r - parseProduct();
			else
				return r;
	}

  public:
	Parser() = default;
	Parser(std::string_view source, SymbolTable const *symbols)
	    : lexer(source), symbolTable(symbols){};

	T parseExpression() { return parseSum(); }
};

int main()
{
	// using Real = FloatingOctuple; // just (high-precision) floating point
	using Real = Number; // mixed exact/floating point calculations

	using_history();
	read_history(historyFile);
	Parser<Real>::SymbolTable symbolTable;
	symbolTable["pi"] = Real::pi();

	while (auto line = getInput())
	{
		try
		{
			auto parser = Parser<Real>(*line, &symbolTable);
			while (!parser.lexer.empty())
			{
				if (parser.lexer.peek(util::Tok::Ident, util::Tok::Assign))
				{
					auto varname = parser.lexer.match(util::Tok::Ident).value;
					parser.lexer.match(util::Tok::Assign);
					auto val = parser.parseExpression();
					symbolTable[std::string(varname)] = val;
				}
				else
				{
					auto val = parser.parseExpression();
					fmt::print("{}\n", val);
				}

				if (!parser.lexer.tryMatch(util::Tok::Semi))
					parser.lexer.match(util::Tok::None);
			}
		}
		catch (util::ParseError const &e)
		{
			fmt::print("error: {}\n", e.what());
		}
	}

	write_history(historyFile);
	fmt::print("\nquit\n");
}
