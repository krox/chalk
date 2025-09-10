#pragma once

#include "chalk/fraction.h"
#include "chalk/parser.h"
#include "fmt/format.h"
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

} // namespace chalk

template <>
struct fmt::formatter<chalk::Expression> : fmt::formatter<std::string>
{
	template <typename FormatContext>
	auto format(chalk::Expression const &ex,
	            FormatContext &ctx) const -> decltype(ctx.out())
	{
		using namespace chalk;

		auto it = ctx.out();
		ex.visit(util::overloaded{
		    [&it](Rational const &x) { it = fmt::format_to(it, "{}", x); },
		    [&it](VariableExpression const &s) {
			    it = fmt::format_to(it, "{}", s);
		    },
		    [&it](SumExpression const &sum) {
			    it = fmt::format_to(it, "(");
			    for (size_t i = 0; i < sum.terms.size(); ++i)
			    {
				    if (i > 0)
					    it = fmt::format_to(it, " + ");
				    it = fmt::format_to(it, "{}*{}", sum.terms[i].first,
				                        sum.terms[i].second);
			    }
			    it = fmt::format_to(it, ")");
		    },
		    [&it](ProductExpression const &prod) {
			    it = fmt::format_to(it, "(");
			    for (size_t i = 0; i < prod.terms.size(); ++i)
			    {
				    if (i > 0)
					    it = fmt::format_to(it, " * ");
				    it = fmt::format_to(it, "{}", prod.terms[i]);
			    }
			    it = fmt::format_to(it, ")");
		    },
		    [&it](PowerExpression const &pow) {
			    it = fmt::format_to(it, "({})^{}", pow.base, pow.exp);
		    },
		    [&it](FunctionExpression const &fun) {
			    it = fmt::format_to(it, "{}(", fun.name);
			    for (size_t i = 0; i < fun.arguments.size(); ++i)
			    {
				    if (i > 0)
					    it = fmt::format_to(it, ", ");
				    it = fmt::format_to(it, "{}", fun.arguments[i]);
			    }
			    it = fmt::format_to(it, ")");
		    }

		});
		return it;
	}
};
