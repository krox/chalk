#include "chalk/symbolic.h"
#include <cassert>
#include <cctype>
#include <fmt/format.h>
#include <optional>
#include <readline/history.h>
#include <readline/readline.h>
#include <string>
using namespace chalk;

/** trim whitespace from both ends of a string */
std::string trim_white(std::string_view s)
{
	size_t i = 0, j = s.size();
	while (i < s.size() && std::isspace(s[i]))
		++i;
	while (j > i && std::isspace(s[j - 1]))
		--j;
	return std::string{s.substr(i, j - i)};
}

/**
 * Get one line of user input from the terminal.
 *   - trims whitespace from both ends
 *   - returns std::nullpt on EOF (i.e. CRTL-D), which is different from ""
 *   - supports history using readline library, storing it in '~/.chalk_history'
 */
std::optional<std::string> get_input()
{
	char *s = readline(">> ");
	if (s == nullptr)
		return std::nullopt;
	auto line = trim_white(s);
	if (!line.empty())
		add_history(s); // add non-trimmed to history, to preserve indentation
	free(s);
	return line;
}

int main()
{
	while (auto line = get_input())
	{
		if (line->empty())
			continue;
		auto ex = Expression(*line);
		fmt::print("{}\n", ex);
	}
	fmt::print("\n");
}
