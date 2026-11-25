#include "CLI/CLI.hpp"
#include "chalk/number.h"
#include "chalk/symbolic.h"
#include "util/linenoise.h"
#include "util/string.h"
#include <fmt/format.h>
#include <string>
using namespace chalk;

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

	util::Linenoise linenoise("chalk");

	// Start REPL
	fmt::print("Welcome to Chalk\n");
	fmt::print("Type expressions, or 'quit' to exit.\n\n");

	while (true)
	{
		auto line = linenoise(">> ");
		if (line.empty())
			break;

		line = std::string(util::trim_white(line));
		if (line == "quit" || line == "exit")
			break;

		try
		{
			auto val = context(line);
			fmt::print("{}\n", val);
		}
		catch (util::ParseError const &e)
		{
			fmt::print("error: {}\n", e.what());
		}
	}

	fmt::print("\nquit\n");
}
