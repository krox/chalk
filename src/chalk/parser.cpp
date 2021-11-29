#include "chalk/parser.h"

#include <algorithm>
#include <cctype>

namespace chalk {

namespace {
bool isoperator(char c)
{
	switch (c)
	{
	case '(':
	case ')':
	case '+':
	case '-':
	case '*':
	case '/':
	case '%':
	case '^':
		return true;
	default:
		return false;
	}
}
} // namespace

std::vector<std::string_view> lex(std::string_view source)
{
	std::vector<std::string_view> tokens;

	while (true)
	{
		// skip whitespace and quit when nothing is left
		while (!source.empty() && std::isspace(source[0]))
			source.remove_prefix(1);
		if (source.empty())
			break;

		// munch one token
		size_t len = 1;
		if (std::isdigit(source[0])) // integer literals
		{
			while (len < source.size() && std::isdigit(source[len]))
				++len;
		}
		else if (std::isalpha(source[0]) || source[0] == '_') // symbols
		{
			while (len < source.size() &&
			       (std::isalnum(source[len]) || source[len] == '_'))
				++len;
		}
		else if (isoperator(source[0])) // operators / parentheses
		{}
		else
			throw std::runtime_error(
			    fmt::format("unexpected character '{}'", source[0]));

		tokens.push_back(source.substr(0, len));
		source.remove_prefix(len);
	}

	return tokens;
}

} // namespace chalk
