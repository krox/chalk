#ifndef CHALK_PARSER_H
#define CHALK_PARSER_H

#include "fmt/format.h"
#include <cassert>
#include <cctype>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace chalk {

namespace parser_detail {

template <typename T, typename F>
T parse_sum(std::vector<std::string_view> const &, size_t &, F &&);
template <typename T, typename F>
T parse_product(std::vector<std::string_view> const &, size_t &, F &&);
template <typename T, typename F>
T parse_primary(std::vector<std::string_view> const &, size_t &, F &&);

} // namespace parser_detail

std::vector<std::string_view> lex(std::string_view source);

/** apparently, std has no clean version of this */
int parse_int(std::string_view s);

template <typename T, typename F> T parse(std::string_view source, F &&f)
{
	using namespace parser_detail;
	auto tokens = lex(source);
	size_t pos = 0;
	auto ex = parse_sum<T>(tokens, pos, f);
	if (pos != tokens.size())
		throw std::runtime_error(
		    fmt::format("dont understand trailing token '{}'", tokens[pos]));
	return ex;
}

namespace parser_detail {

template <typename T, typename F>
T parse_sum(std::vector<std::string_view> const &tokens, size_t &pos, F &&f)
{
	auto ex = parse_product<T>(tokens, pos, f);
	while (pos != tokens.size())
	{
		if (tokens[pos] == "+")
		{
			++pos;
			ex += parse_product<T>(tokens, pos, f);
		}
		else if (tokens[pos] == "-")
		{
			++pos;
			ex -= parse_product<T>(tokens, pos, f);
		}
		else
			break;
	}
	return ex;
}

template <typename T, typename F>
T parse_product(std::vector<std::string_view> const &tokens, size_t &pos, F &&f)
{
	auto ex = parse_primary<T>(tokens, pos, f);
	while (pos != tokens.size())
	{
		if (tokens[pos] == "*")
		{
			++pos;
			ex *= parse_primary<T>(tokens, pos, f);
		}
		else if (tokens[pos] == "/")
		{
			++pos;
			ex /= parse_primary<T>(tokens, pos, f);
		}
		else
			break;
	}
	return ex;
}

template <typename T, typename F>
T parse_primary(std::vector<std::string_view> const &tokens, size_t &pos, F &&f)
{
	if (pos == tokens.size())
		throw std::runtime_error("unexpected end of input");
	assert(!tokens[pos].empty());

	T ex;
	if (std::isalpha(tokens[pos][0]) || std::isdigit(tokens[pos][0]) ||
	    tokens[pos][0] == '_')
		ex = f(tokens[pos++]);

	else if (tokens[pos] == "(")
	{
		++pos;
		ex = parse_sum<T>(tokens, pos, f);
		if (pos == tokens.size() || tokens[pos] != ")")
			throw std::runtime_error(
			    fmt::format("missing ')' after '{}'", tokens[pos - 1]));
		++pos;
	}

	if (pos != tokens.size())
		if (tokens[pos] == "^")
		{
			++pos;
			if (pos == tokens.size() || !std::isdigit(tokens[pos][0]))
				throw std::runtime_error("missing/invalid exponent after '^'");
			int b = parse_int(tokens[pos++]);
			ex = pow(ex, b);
		}

	return ex;
}

} // namespace parser_detail

} // namespace chalk

#endif
