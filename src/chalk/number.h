#pragma once

#include "chalk/floating.h"
#include "chalk/fraction.h"
#include "chalk/integer.h"
#include "chalk/numerics.h"
#include "chalk/polynomial.h"
#include "chalk/quadratic.h"
#include <variant>

namespace chalk {

// Represents a real number, either in exact form or as a (high-precision)
// floating point number. Automatically converts to the latter if the result
// of some operation can not be represented exactly.
// TODO: use more powerful 'exact' type, or even template over Exact/Approx
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

inline std::optional<Number> trySqrt(Number const &a)
{
	if (a.isExact())
	{
		if (auto r = trySqrt(a.exact()); r)
			return Number(*std::move(r));
		return trySqrt(a.toApprox());
	}

	return trySqrt(a.approx());
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

// comparisons

inline bool operator==(Number const &a, Number const &b)
{
	if (a.isExact() && b.isExact())
		return a.exact() == b.exact();

	// maybe 'return false' would be more appropriate mathematically
	return a.toApprox() == b.toApprox();
}

inline bool operator==(Number const &a, int b)
{
	if (a.isExact())
		return a.exact() == b;
	else
		return a.approx() == b;
}

// misc functions

// Polynomial root finding. Can give exact answers for some simple cases
inline std::vector<Number> roots(Polynomial<Number> const &poly)
{
	assert(poly.degree() >= 0); // cant compute roots of zero-polynomial
	if (poly.degree() == 0)     // constant (non-zero) -> no roots
		return {};

	// constant term (exact) zero -> exact zero root
	// NOTE: this is true even if the other coeffs are not exact
	if (poly[0].isExact() && poly[0] == 0)
	{
		auto tmp = std::vector<Number>(&poly.coefficients()[1],
		                               &*poly.coefficients().end());

		auto r = roots(Polynomial<Number>(std::move(tmp)));
		r.push_back(Number(0));
		return r;
	}

	// explicit formulas for linear and quadratic case
	// (these give exact results if the coeffs are all exact)
	if (poly.degree() == 1)
		return {-poly[0] / poly[1]};
	if (poly.degree() == 2)
	{
		auto p = poly[1] / poly[2];
		auto q = poly[0] / poly[2];
		auto d = trySqrt(p * p / 4 - q);
		if (!d)
			return {};

		// roots are very close together -> just return one
		// if (abs(*d) < 1.0e-12)
		//	return {-0.5 * p};

		// solutions are now -p/2 +- d. But it is numerically advantageous to
		// determine only one root this way, and the other by Vieta lemma
		if (p.toApprox() > 0)
		{
			auto x = -p / 2 - *d;
			return {x, q / x};
		}
		else
		{
			auto x = -p / 2 + *d;
			return {q / x, x};
		}
	}

	// otherwise, fall back to purely numerical solution
	using Real = Number::approx_type;
	using Complex = util::complex<Real>;
	std::vector<Complex> coeffs;
	coeffs.reserve(poly.degree() + 1);
	for (auto &c : poly.coefficients())
		coeffs.emplace_back(c.toApprox());
	auto eps = Real(1e-50);
	auto rc = rootsDurandKerner(Polynomial<Complex>(std::move(coeffs)), eps);
	std::vector<Number> r;
	for (auto &x : rc)
		if (abs(x.im) < eps)
			r.emplace_back(x.re);
	std::sort(r.begin(), r.end(), [](auto const &a, auto const &b) {
		return a.to_double() < b.to_double();
	});

	return r;
}

template <> struct RingTraits<Number> : RingTraitsSimple<Number>
{};

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
