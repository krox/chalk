#pragma once

#include "chalk/fraction.h"
#include "chalk/parser.h"
#include "fmt/format.h"
#include "util/functional.h"
#include "util/json.h" // for 'overloaded'. should be moved eventually
#include <memory>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace chalk {

class ExpressionVariant;

class Expression
{
	std::shared_ptr<const ExpressionVariant> ptr_;

  public:
	Expression();

	// low-level constructor for internal use
	explicit Expression(std::shared_ptr<const ExpressionVariant> p);

	// user-facing constructors
	explicit Expression(int value);
	explicit Expression(Rational const &value);
	static Expression variable(std::string_view name);
	static Expression integer(std::string_view s);
	static Expression sum(std::span<const Expression> terms);
	static Expression product(std::span<const Expression> terms);
	static Expression power(Expression const &base, int exp);

	template <class R = void, class Visitor> R visit(Visitor &&vis) const;

	ExpressionVariant const *data() const { return ptr_.get(); }
};

// might have dedicated classes later, but for now this is enough
using ConstantExpression = Rational;
using VariableExpression = std::string;

struct SumExpression
{
	std::vector<std::pair<Rational, Expression>> terms;
};

struct ProductExpression
{
	std::vector<Expression> terms;
};

struct PowerExpression
{
	Expression base;
	int exp;
};

struct FunctionExpression
{
	std::string name;
	std::vector<Expression> arguments;
};

class ExpressionVariant
    : public std::variant<Rational, std::string, SumExpression,
                          ProductExpression, PowerExpression,
                          FunctionExpression>
{
  public:
	using variant::variant;
};

template <class R, class Visitor> R Expression::visit(Visitor &&vis) const
{
	return std::visit(std::forward<Visitor>(vis), *ptr_);
}

Expression operator+(Expression const &a);
Expression operator-(Expression const &a);
Expression operator+(Expression const &a, Expression const &b);
Expression operator-(Expression const &a, Expression const &b);
Expression operator*(Expression const &a, Expression const &b);
Expression operator/(Expression const &a, Expression const &b);
Expression pow(Expression const &base, int exp);

bool is_zero(Expression const &);
bool is_one(Expression const &);
bool is_negative(Expression const &);
bool need_parens_product(Expression const &);
bool need_parens_power(Expression const &);
void print(Expression const &a,
           util::function_view<void(std::string_view)> sink,
           format_type format = {}, int depth = INT_MAX);
} // namespace chalk

template <> struct fmt::formatter<chalk::Expression>
{
  private:
	chalk::format_type type_ = chalk::format_type::text;
	int depth_ = INT_MAX;

  public:
	constexpr auto parse(auto &ctx) -> decltype(ctx.begin())
	{
		auto it = ctx.begin();
		auto end = ctx.end();

		if (it != end && '0' <= *it && *it <= '9')
		{
			depth_ = 0;
			while (it != end && '0' <= *it && *it <= '9')
			{
				depth_ = depth_ * 10 + (*it - '0');
				++it;
			}
		}

		if (it != end)
		{
			if (*it == 't')
			{
				type_ = chalk::format_type::text;
				++it;
			}
			else if (*it == 'd')
			{
				type_ = chalk::format_type::debug;
				++it;
			}
			/*else if (*it == 'l')
			{
			    type_ = format_type::latex;
			    ++it;
			}*/
		}
		// check that the end is reached
		if (it != end && *it != '}')
			throw fmt::format_error("invalid format");
		return it;
	}

	auto format(chalk::Expression const &a,
	            auto &ctx) const -> decltype(ctx.out())
	{
		auto it = ctx.out();
		auto sink = [&it](std::string_view s) {
			for (char c : s)
				*it++ = c;
		};
		print(a, sink, type_, depth_);
		return it;
	}
};
