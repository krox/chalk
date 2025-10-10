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
	util::hash_map<std::string, Expression> var_table;

  public:
	Parser() : lexer({}) {}

	Expression operator()(std::string_view source)
	{
		lexer = util::Parser(source);

		std::string_view id;

		if (lexer.ident("let"))
		{
			id = lexer.expect_ident();
			lexer.expect('=');
		}

		auto expr = Expression::parse(lexer);
		if (!lexer.end())
			throw util::ParseError("unexpected trailing input");
		expr = substitute(expr, var_table);

		if (!id.empty())
			var_table[id] = expr;

		// safety measure: remove (potentially dangling) reference
		lexer = util::Parser("");
		return expr;
	}
};

int main()
{
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
