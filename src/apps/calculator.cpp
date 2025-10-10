#include "CLI/CLI.hpp"
#include "chalk/number.h"
#include "chalk/symbolic.h"
#include <fmt/format.h>
#include <optional>
#include <readline/history.h>
#include <readline/readline.h>
#include <string>
using namespace chalk;

static char historyFile[] = ".chalk_history";

// Get one line of user input from the terminal using readline library
//   - trims whitespace from both ends
//   - returns std::nullpt on EOF (i.e. CRTL-D), which is different from ""
//   - TODO: switch to 'linenoise' library.
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

int main(int argc, char **argv)
{
	CLI::App app{"Chalk Calculator - Symbolic computation calculator"};

	std::string filename;
	app.add_option("-f,--file", filename,
	               "Execute statements from file before starting REPL");

	try
	{
		app.parse(argc, argv);
	}
	catch (CLI::ParseError const &e)
	{
		return app.exit(e);
	}

	using_history();
	read_history(historyFile);

	SymbolicContext context;

	// Execute file if provided
	if (!filename.empty())
	{
		try
		{
			context.execute_file(filename);
		}
		catch (util::ParseError const &e)
		{
			fmt::print("error: {}\n", e.what());
			return 1;
		}
		fmt::print("\n"); // Separator between file output and REPL
	}

	// Start REPL
	fmt::print("Welcome to Chalk\n");
	fmt::print("Type expressions, or 'quit' to exit.\n\n");

	while (auto line = getInput())
	{
		try
		{
			auto val = context(*line);
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
