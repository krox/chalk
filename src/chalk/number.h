#pragma once

#include "chalk/floating.h"
#include "chalk/fraction.h"
#include "chalk/integer.h"
#include "chalk/quadratic.h"
#include <variant>

namespace chalk {

// Represents a real number, either in exact form or as a (high-precision)
// floating point number. Automatically converts to the latter if the result
// of some operation can not be represented exactly.
class Number
{
  public:
	using exact_type = Quadratic<Rational>;
	using approx_type = FloatingOctuple;

  private:
	std::variant<exact_type, approx_type> value_;

  public:
	Number() = default;
	Number(int value) : value_(exact_type(value)) {}
	Number(Integer value) : value_(exact_type(Rational(std::move(value)))) {}
	Number(exact_type value) : value_(std::move(value)) {}
	Number(double value) : value_(approx_type(value)) {}
	Number(approx_type value) : value_(std::move(value)) {}

	static Number pi() { return approx_type::pi(); }

	bool isExact() const { return std::holds_alternative<exact_type>(value_); }
	exact_type const &exact() const { return std::get<exact_type>(value_); }
	approx_type const &approx() const { return std::get<approx_type>(value_); }
	approx_type toApprox() const
	{
		if (isExact())
			return (approx_type)exact();
		else
			return approx();
	}

	double to_double() const
	{
		if (isExact())
			return (double)exact();
		else
			return approx().to_double();
	}
};

// unary Number

inline Number operator-(Number const &a)
{
	if (a.isExact())
		return -a.exact();
	else
		return -a.approx();
}

inline Number sqrt(Number const &a)
{
	if (a.isExact())
		if (auto r = trySqrt(a.exact()); r)
			return *r;
	return sqrt(a.toApprox());
}

inline Number pow(Number const &a, int e)
{
	if (a.isExact())
		return pow(a.exact(), e);
	else
		return pow(a.approx(), e);
}

// binary Number

#define CHALK_DEFINE_NUMBER_BINARY(fun, tryFun)                                \
	inline Number fun(Number const &a, Number const &b)                        \
	{                                                                          \
		if (a.isExact())                                                       \
		{                                                                      \
			if (b.isExact())                                                   \
			{                                                                  \
				if (auto r = tryFun(a.exact(), b.exact()); r)                  \
					return *r;                                                 \
				else                                                           \
					return fun(a.toApprox(), b.toApprox());                    \
			}                                                                  \
			else                                                               \
				return fun(a.toApprox(), b.approx());                          \
		}                                                                      \
		else                                                                   \
		{                                                                      \
			if (b.isExact())                                                   \
				return fun(a.approx(), b.toApprox());                          \
			else                                                               \
				return fun(a.approx(), b.approx());                            \
		}                                                                      \
	}

CHALK_DEFINE_NUMBER_BINARY(operator+, tryAdd)
CHALK_DEFINE_NUMBER_BINARY(operator-, trySub)
CHALK_DEFINE_NUMBER_BINARY(operator*, tryMul)
CHALK_DEFINE_NUMBER_BINARY(operator/, tryDiv)

} // namespace chalk

template <> struct fmt::formatter<chalk::Number>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Number &x, FormatContext &ctx)
	{
		if (x.isExact())
			return format_to(ctx.out(), "{}", x.exact());
		else
			return format_to(ctx.out(), "{}", x.approx().to_double());
	}
};
