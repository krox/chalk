#pragma once

#include "chalk/fraction.h"
#include "chalk/parser.h"
#include "fmt/format.h"
#include "util/functional.h"
#include "util/json.h" // for 'overloaded'. should be moved eventually
#include <compare>
#include <memory>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace chalk {

class Expression;

struct ConstantExpression
{
	Rational value;

	std::strong_ordering operator<=>(ConstantExpression const &other) const;
	bool operator==(ConstantExpression const &other) const;
};

struct VariableExpression
{
	std::string name;

	std::strong_ordering operator<=>(VariableExpression const &other) const;
	bool operator==(VariableExpression const &other) const;
};

struct SumExpression
{
	std::vector<std::pair<Rational, Expression>> terms;

	std::strong_ordering operator<=>(SumExpression const &other) const;
	bool operator==(SumExpression const &other) const;
};

struct ProductExpression
{
	std::vector<Expression> terms;

	std::strong_ordering operator<=>(ProductExpression const &other) const;
	bool operator==(ProductExpression const &other) const;
};

struct PowerExpression
{
	util::value_ptr<Expression> base; // not null
	int exp;

	std::strong_ordering operator<=>(PowerExpression const &other) const;

	bool operator==(PowerExpression const &other) const;

	PowerExpression(Expression const &base, int exp)
	    : base(util::make_value<Expression>(base)), exp(exp)
	{}
};

struct FunctionExpression
{
	std::string name;
	std::vector<Expression> arguments;

	std::strong_ordering operator<=>(FunctionExpression const &other) const;
	bool operator==(FunctionExpression const &other) const;
};
using ExpressionVariant =
    std::variant<ConstantExpression, VariableExpression, SumExpression,
                 ProductExpression, PowerExpression, FunctionExpression>;
class Expression : public ExpressionVariant
{
	Expression(ExpressionVariant const &v) : ExpressionVariant(v) {}
	Expression(ExpressionVariant &&v) : ExpressionVariant(std::move(v)) {}

  public:
	// default = ConstantExpression(0)
	Expression() = default;

	// convenience conversions
	explicit Expression(int value);
	explicit Expression(Rational const &value);

	// constructors of explicit types
	// note: due to normalization rules, the result might not actually be of
	// that type. E.g.:
	//   sum({0, X}) -> X
	static Expression variable(std::string_view name);
	static Expression integer(std::string_view s);
	static Expression sum(std::span<const Expression> terms);
	static Expression product(std::span<const Expression> terms);
	static Expression power(Expression const &base, int exp);
	static Expression function(std::string_view fun,
	                           std::span<const Expression> args);

	// this defines an arbitrary syntactic order on expressions, nothing
	// mathematically meaningful
	std::strong_ordering operator<=>(Expression const &other) const;

	// syntactic equality
	// (i.e. sufficient but not necessary for mathematical equivalence)
	bool operator==(Expression const &other) const;
};

// TODO: move this into 'Rational' itself. should be useful more generally, also
// related to finally getting a small-size-optimization for 'Integer'.
// inline const Rational one_rat = Rational(1);

// basic math operations
Expression operator+(Expression const &a);
Expression operator-(Expression const &a);
Expression operator+(Expression const &a, Expression const &b);
Expression operator-(Expression const &a, Expression const &b);
Expression operator*(Expression const &a, Expression const &b);
Expression operator/(Expression const &a, Expression const &b);
Expression pow(Expression const &base, int exp);

bool is_zero(Expression const &);
bool is_one(Expression const &);

// helpers for pretty printing
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
