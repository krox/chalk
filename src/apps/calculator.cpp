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
	using Function = std::function<T(std::vector<T> const &)>;
	using VarTable = std::map<std::string, T, std::less<>>;
	using FunTable = std::map<std::string, Function, std::less<>>;

	using Tok = util::Tok;
	using Token = util::Token;

  private:
	VarTable const *varTable = nullptr;
	FunTable const *funTable = nullptr;

	T parsePrimary()
	{
		if (lexer.try_match("("))
		{
			auto r = parseSum();
			lexer.match(")");
			return r;
		}
		else if (auto tok = lexer.try_match(Tok::integer()); tok)
		{
			return T(Integer(tok->value));
		}
		else if (auto tok = lexer.try_match(Tok::floating()); tok)
		{
			return T(FloatingOctuple(tok->value));
		}
		else if (auto tok = lexer.try_match(Tok::ident()); tok)
		{
			if (lexer.try_match("(")) // function call
			{
				if (funTable != nullptr)
					if (auto it = funTable->find(tok->value);
					    it != funTable->end())
					{
						std::vector<T> args;
						while (!lexer.try_match(")"))
						{
							args.push_back(parseSum());
							if (!lexer.peek(")"))
								lexer.match(",");
						}
						return it->second(args);
					}

				throw util::ParseError(
				    fmt::format("unknown function '{}'", tok->value));
			}
			else // variable / constant
			{
				if (varTable != nullptr)
					if (auto it = varTable->find(tok->value);
					    it != varTable->end())
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
		bool neg = (bool)lexer.try_match("-");
		T r = parsePrimary();
		if (lexer.try_match("^"))
		{
			if (lexer.try_match("-"))
				r = pow(r, -util::parse_int(lexer.match(Tok::integer()).value));
			else
				r = pow(r, util::parse_int(lexer.match(Tok::integer()).value));
		}
		return neg ? -r : r;
	}

	T parseProduct()
	{
		T r = parsePower();
		while (true)
			if (lexer.try_match("*"))
				r = r * parsePower();
			else if (lexer.try_match("/"))
				r = r / parsePower();
			else
				return r;
	}

	T parseSum()
	{
		T r = parseProduct();
		while (true)
			if (lexer.try_match("+"))
				r = r + parseProduct();
			else if (lexer.try_match("-"))
				r = r - parseProduct();
			else
				return r;
	}

  public:
	Parser() = default;
	Parser(std::string_view source, VarTable const *vars, FunTable const *funs)
	    : lexer(source), varTable(vars), funTable(funs)
	{}

	T parseExpression() { return parseSum(); }
};

int main()
{
	// using Real = FloatingOctuple; // just (high-precision) floating point
	using Real = Number; // mixed exact/floating point calculations

	using_history();
	read_history(historyFile);

	Parser<Real>::VarTable varTable;
	varTable["pi"] = Real::pi();

	Parser<Real>::FunTable funTable;
	funTable["sqrt"] = [](std::vector<Real> const &args) {
		if (args.size() != 1)
			throw util::ParseError(
			    "wrong number of arguments for function sqrt");
		return sqrt(args[0]);
	};
	funTable["polyRoots"] = [](std::vector<Real> const &args) {
		auto poly = Polynomial<Real>(args);
		auto r = roots(poly);
		for (auto const &x : r)
			fmt::print("x = {}\n", x);
		return Real(0);
	};

	while (auto line = getInput())
	{
		try
		{
			auto parser = Parser<Real>(*line, &varTable, &funTable);
			while (!parser.lexer.empty())
			{
				if (parser.lexer.peek(util::Tok::ident(), "="))
				{
					auto varname = parser.lexer.match(util::Tok::ident()).value;
					parser.lexer.match("=");
					auto val = parser.parseExpression();
					varTable[std::string(varname)] = val;
				}
				else
				{
					auto val = parser.parseExpression();
					fmt::print("{}\n", val);
				}

				if (!parser.lexer.try_match(";"))
					parser.lexer.match(util::Tok::none());
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
