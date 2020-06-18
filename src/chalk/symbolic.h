#ifndef CHALK_SYMBOLIC_H
#define CHALK_SYMBOLIC_H

#include "chalk/fraction.h"
#include "chalk/parser.h"
#include "fmt/format.h"
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace chalk {

/**
 * - Internally everything is represented as 'shared_ptr<const BaseExpression>'
 *   with free functions ('make_constant', 'make_sum', 'make_function', ...)
 * - For the user, we provide a thin wrapper 'Expression' which
 *   provides operator overloading ('a+b', 'exp(a)', 'a+5', ...)
 */
class BaseExpression
{
  public:
	virtual ~BaseExpression(){};

	/**
	 * parens = 0 -> no parentheses
	 * parens = 1 -> parentheses around sums
	 * parens = 2 -> parentheses around sums/products
	 * parens = 3 -> parentheses around sums/products/powers
	 */
	virtual void to_string(fmt::memory_buffer &buf, int parens = 0) const = 0;
	// virtual void to_latex(fmt::memory_buffer &) const = 0;
};

using ExpressionPtr = std::shared_ptr<BaseExpression const>;

/**
 * Create expressions.
 *   - These already do basic simplifications (like 'x+0' -> 'x').
 *   - 'make_sum' can return arbitrary type of expression, not only
 *     'SumExpression'. Therefore these methods cant be written as constructors.
 */
ExpressionPtr make_variable(std::string_view);
ExpressionPtr make_constant(int);
ExpressionPtr make_constant(Rational const &);
ExpressionPtr make_sum(std::vector<ExpressionPtr> const &);
ExpressionPtr make_product(std::vector<ExpressionPtr> const &);
ExpressionPtr make_power(ExpressionPtr const &, ExpressionPtr const &);
ExpressionPtr make_function(std::string_view, ExpressionPtr const &);

/** This type should be used in user code */
struct Expression
{
	std::shared_ptr<const BaseExpression> ex;

	Expression() = default;
	Expression(std::shared_ptr<const BaseExpression> expr) : ex(std::move(expr))
	{}

	explicit Expression(int);
	explicit Expression(Rational const &);
	explicit Expression(std::string_view); // parses expressions

	std::string to_string() const;
};

inline Expression operator-(Expression const &a)
{
	return make_product({make_constant(-1), a.ex});
}
inline Expression operator+(Expression const &a, Expression const &b)
{
	return make_sum({a.ex, b.ex});
}
inline Expression operator-(Expression const &a, Expression const &b)
{
	return make_sum({a.ex, make_product({make_constant(-1), b.ex})});
}
inline Expression operator*(Expression const &a, Expression const &b)
{
	return make_product({a.ex, b.ex});
}
inline Expression operator/(Expression const &a, Expression const &b)
{
	return make_product({a.ex, make_power(b.ex, make_constant(-1))});
}
inline Expression pow(Expression const &a, Expression const &b)
{
	return make_power(a.ex, b.ex);
}

inline Expression pow(Expression const &a, int b)
{
	return make_power(a.ex, make_constant(b));
}

inline void operator+=(Expression &a, Expression const &b) { a = a + b; }
inline void operator-=(Expression &a, Expression const &b) { a = a - b; }
inline void operator*=(Expression &a, Expression const &b) { a = a * b; }
inline void operator/=(Expression &a, Expression const &b) { a = a / b; }

inline Expression::Expression(int value) : ex(make_constant(value)) {}
inline Expression::Expression(Rational const &value) : ex(make_constant(value))
{}
inline Expression::Expression(std::string_view s)
{
	auto f = [](std::string_view token) -> Expression {
		if (std::isdigit(token[0]))
			return make_constant(parse_int(token));
		else if (std::isalpha(token[0]))
			return make_variable(token);
		else
			throw std::runtime_error("cant process token");
	};
	ex = parse<Expression>(s, f).ex;
}

inline std::string Expression::to_string() const
{
	assert(ex);
	fmt::memory_buffer buf;
	ex->to_string(buf, 0);
	return fmt::to_string(buf);
}

} // namespace chalk

template <>
struct fmt::formatter<chalk::Expression> : fmt::formatter<std::string>
{
	template <typename FormatContext>
	auto format(chalk::Expression const &ex, FormatContext &ctx)
	{
		return format_to(ctx.out(), ex.to_string());
	}
};

#endif
