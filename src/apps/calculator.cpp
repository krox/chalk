#include "chalk/number.h"
#include "chalk/symbolic.h"
#include "util/lexer.h"
#include "util/string.h"
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

/**
 * Get one line of user input from the terminal using readline library
 *   - trims whitespace from both ends
 *   - returns std::nullpt on EOF (i.e. CRTL-D), which is different from ""
 *   - TODO: switch to 'linenoise' library.
 */
std::optional<std::string> getInput()
{
	char *s = readline(">> ");
	if (s == nullptr)
		return std::nullopt;
	auto line = std::string(util::trim_white(s));
	if (!line.empty())
		add_history(s); // add non-trimmed to history, to preserve indentation
	free(s);
	return line;
}

class Parser
{
  public:
	// yeah. the names 'Lexer' and 'Parser' are used very inconsistently..
	util::Parser lexer;

  private:
	Expression parse_primary()
	{
		if (lexer.match('('))
		{
			auto r = parse_sum();
			lexer.expect(')');
			return r;
		}
		else if (auto tok = lexer.integer(); !tok.empty())
		{
			return Expression::integer(tok);
		}
		/*else if (auto tok = lexer.try_match(Tok::floating()); tok)
		{
		    return T(FloatingOctuple(tok->value));
		}*/
		else if (auto ident = lexer.ident(); !ident.empty())
		{
			if (lexer.match('(')) // function call
			{
				/*if (funTable != nullptr)
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
				    }*/

				throw util::ParseError(
				    fmt::format("unknown function '{}'", ident));
			}
			else // variable / constant
			{
				/*if (varTable != nullptr)
				    if (auto it = varTable->find(tok->value);
				        it != varTable->end())
				        return it->second;*/

				return Expression::variable(ident);
			}
		}
		else
			throw util::ParseError(fmt::format("unexpected token in input"));
	}

	Expression parse_power()
	{
		bool neg = (bool)lexer.match('-');
		auto r = parse_primary();
		/*if (lexer.match('^') || lexer.match("**"))
		{
		    if (lexer.match('-'))
		        r = pow(r, -util::parse_int(lexer.match(Tok::integer()).value));
		    else
		        r = pow(r, util::parse_int(lexer.match(Tok::integer()).value));
		}*/
		return neg ? -r : r;
	}

	Expression parse_product()
	{
		auto r = parse_power();
		while (true)
			if (lexer.match('*'))
				r = r * parse_power();
			else if (lexer.match('/'))
				r = r / parse_power();
			else
				return r;
	}

	Expression parse_sum()
	{
		auto r = parse_product();
		while (true)
			if (lexer.match('+'))
				r = r + parse_product();
			else if (lexer.match('-'))
				r = r - parse_product();
			else
				return r;
	}

  public:
	Parser() : lexer({}) {}

	Expression operator()(std::string_view source)
	{
		lexer = util::Parser(source);
		auto r = parse_sum();
		if (!lexer.end())
			throw util::ParseError("unexpected trailing input");
		lexer = util::Parser(""); // clear lexer to avoid dangling references
		return r;
	}
};

int main()
{
	// using Real = FloatingOctuple; // just (high-precision) floating point

	auto expr = Expression(5) + Expression::variable("x") + Expression(7);
	fmt::print("expr = {}\n", expr);
	using_history();
	read_history(historyFile);

	// Parser<Real>::VarTable varTable;
	// varTable["pi"] = Real::pi();

	/*Parser<Real>::FunTable funTable;
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
	};*/
	Parser parser;

	while (auto line = getInput())
	{
		try
		{
			auto val = parser(*line);
			fmt::print("{}\n", val);
		}
		catch (util::ParseError const &e)
		{
			fmt::print("error: {}\n", e.what());
		}
	}

	write_history(historyFile);
	fmt::print("\nquit\n");
}
