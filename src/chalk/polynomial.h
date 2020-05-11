#ifndef CHALK_POLYNOMIAL_H
#define CHALK_POLYNOMIAL_H

#include "fmt/format.h"
#include <array>
#include <cassert>
#include <climits>
#include <utility>
#include <vector>

/** helpers to make Polynomial<double> work */
inline bool is_negative(double x) { return x < 0; }

namespace chalk {

template <typename R> class Polynomial;
template <typename R>
Polynomial<R> truncate(Polynomial<R> const &poly, int order);

/**
 * Univariate polynomial with dense storage and coefficients in R.
 * Can also be used truncated power series.
 */
template <typename R> class Polynomial
{
	std::vector<R> coefficients_;
	int max_order_ = INT_MAX;
	inline static const R zero_ = R(0); // stupid workaround...

	void cleanup()
	{
		assert(max_order_ >= 0);
		while (!coefficients_.empty() &&
		       ((int)coefficients_.size() - 1 > max_order_ ||
		        coefficients_.back() == R(0)))
			coefficients_.pop_back();
	}

  public:
	Polynomial() = default;
	Polynomial(int c) : coefficients_{R(c)} { cleanup(); }
	Polynomial(const R &c) : coefficients_{c} { cleanup(); }
	Polynomial(std::vector<R> coefficients, int max_order = INT_MAX)
	    : coefficients_(std::move(coefficients)), max_order_(max_order)
	{
		cleanup();
	}
	static Polynomial generator() { return Polynomial({R(0), R(1)}); }

	/** highest power with non-zero coefficient (convention: deg(0) = -1) */
	int degree() const { return (int)coefficients_.size() - 1; }

	int max_order() const { return max_order_; }

	/** highest is non-zero */
	std::vector<R> const &coefficients() const { return coefficients_; }

	/** array-like access to coefficients */
	R const &operator[](int i) const
	{
		if (i < 0 || i > degree())
			return zero_;
		return coefficients_[i];
	}

	/** evaluate the polynomial */
	R operator()(R const &x) const
	{
		R y = 0;
		for (int i = degree(); i >= 0; --i)
			y = y * x + coefficients_[i];
		return y;
	}

	void operator*=(R const &b)
	{
		for (auto &a : coefficients_)
			a *= b;
		cleanup();
	}
	void operator/=(R const &b)
	{
		for (auto &a : coefficients_)
			a /= b;
		cleanup();
	}
	void operator*=(int b)
	{
		for (auto &a : coefficients_)
			a *= b;
		cleanup();
	}
	void operator/=(int b)
	{
		for (auto &a : coefficients_)
			a /= b;
		cleanup();
	}

	friend Polynomial<R> truncate<R>(Polynomial<R> const &, int);
};

template <typename R>
Polynomial<R> truncate(Polynomial<R> const &poly, int order)
{
	assert(order >= 0);
	if (order >= poly.max_order())
		return poly;
	else
		return Polynomial<R>(poly.coefficients(), order);
}

template <typename R>
bool operator==(Polynomial<R> const &a, Polynomial<R> const &b)
{
	// assert(a.max_order() == INT_MAX);
	// assert(b.max_order() == INT_MAX);
	return a.coefficients == b.coefficients;
}

template <typename R> bool operator==(Polynomial<R> const &a, int b)
{
	// assert(a.max_order() == INT_MAX);
	if (b == 0)
		return a.coefficients().empty();
	return a.coefficients().size() == 1 && a.coefficients()[0] == b;
}

/** unary operations */
template <typename R> Polynomial<R> operator-(Polynomial<R> const &a)
{
	auto coeffs = a.coefficients();
	for (auto &c : coeffs)
		c = -c;
	return Polynomial<R>(std::move(coeffs), a.max_order());
}

/** binary operations (Polynomial <-> Polynomial) */
template <typename R>
Polynomial<R> operator+(Polynomial<R> const &a, Polynomial<R> const &b)
{
	auto c = a.coefficients();
	c.resize(std::max(a.degree(), b.degree()) + 1);
	for (int i = 0; i <= b.degree(); ++i)
		c[i] += b[i];
	return Polynomial<R>(std::move(c), std::min(a.max_order(), b.max_order()));
}
template <typename R>
Polynomial<R> operator-(Polynomial<R> const &a, Polynomial<R> const &b)
{
	auto c = a.coefficients();
	c.resize(std::max(a.degree(), b.degree()) + 1);
	for (int i = 0; i <= b.degree(); ++i)
		c[i] -= b[i];
	return Polynomial<R>(std::move(c), std::min(a.max_order(), b.max_order()));
}
template <typename R>
Polynomial<R> operator*(Polynomial<R> const &a, Polynomial<R> const &b)
{
	auto coeffs = std::vector<R>(a.degree() + b.degree() + 1);
	for (int i = 0; i <= a.degree(); ++i)
		for (int j = 0; j <= b.degree(); ++j)
			if (i + j <= std::min(a.max_order(), b.max_order()))
				coeffs[i + j] += a[i] * b[j];
	return Polynomial<R>(std::move(coeffs),
	                     std::min(a.max_order(), b.max_order()));
}

/** binary operations (Polynomial <-> R) */
template <typename R>
Polynomial<R> operator+(Polynomial<R> const &a, R const &b)
{
	auto coeffs = a.coefficients();
	if (coeffs.empty())
		coeffs.resize(1, b);
	else
		coeffs[0] += b;
	return Polynomial<R>(std::move(coeffs), a.max_order());
}
template <typename R>
Polynomial<R> operator-(Polynomial<R> const &a, R const &b)
{
	auto coeffs = a.coefficients();
	if (coeffs.empty())
		coeffs.resize(1, b);
	else
		coeffs[0] -= b;
	return Polynomial<R>(std::move(coeffs), a.max_order());
}
template <typename R>
Polynomial<R> operator*(Polynomial<R> const &a, R const &b)
{
	auto coeffs = std::vector<R>(a.degree() + 1);
	for (int i = 0; i <= a.degree(); ++i)
		coeffs[i] = a[i] * b;
	return Polynomial<R>(std::move(coeffs), a.max_order());
}
template <typename R>
Polynomial<R> operator/(Polynomial<R> const &a, R const &b)
{
	auto coeffs = std::vector<R>(a.degree() + 1);
	for (int i = 0; i <= a.degree(); ++i)
		coeffs[i] = a[i] / b;
	return Polynomial<R>(std::move(coeffs), a.max_order());
}

/** binary operations (Polynomial <-> int) */
template <typename R> Polynomial<R> operator+(Polynomial<R> const &a, int b)
{
	auto coeffs = a.coefficients();
	if (coeffs.empty())
		coeffs.resize(1, R(b));
	else
		coeffs[0] += b;
	return Polynomial<R>(std::move(coeffs), a.max_order());
}
template <typename R> Polynomial<R> operator-(Polynomial<R> const &a, int b)
{
	auto coeffs = a.coefficients();
	if (coeffs.empty())
		coeffs.resize(1, R(b));
	else
		coeffs[0] -= b;
	return Polynomial<R>(std::move(coeffs), a.max_order());
}
template <typename R> Polynomial<R> operator*(Polynomial<R> const &a, int b)
{
	auto coeffs = std::vector<R>(a.degree() + 1);
	for (int i = 0; i <= a.degree(); ++i)
		coeffs[i] = a[i] * b;
	return Polynomial<R>(std::move(coeffs), a.max_order());
}
template <typename R> Polynomial<R> operator/(Polynomial<R> const &a, int b)
{
	auto coeffs = std::vector<R>(a.degree() + 1);
	for (int i = 0; i <= a.degree(); ++i)
		coeffs[i] = a[i] / b;
	return Polynomial<R>(std::move(coeffs), a.max_order());
}

/** op-assigns for convenience (could be optimized...) */
template <typename R> void operator+=(Polynomial<R> &a, Polynomial<R> const &b)
{
	a = a + b;
}
template <typename R> void operator-=(Polynomial<R> &a, Polynomial<R> const &b)
{
	a = a - b;
}
template <typename R> void operator*=(Polynomial<R> &a, Polynomial<R> const &b)
{
	a = a * b;
}

} // namespace chalk

template <typename R> struct fmt::formatter<chalk::Polynomial<R>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Polynomial<R> &poly, FormatContext &ctx)
	    -> decltype(ctx.out())
	{
		auto it = ctx.out();

		// no terms -> output "0"
		if (poly.coefficients().empty())
			it = format_to(it, "0");

		// otherwise list the terms with " + " inbetween
		if (poly.max_order() == INT_MAX)
		{
			bool first = true;
			for (int i = poly.degree(); i >= 0; --i)
			{
				auto c = poly[i];
				if (c == 0)
					continue;

				if (is_negative(c))
				{
					if (first)
						it = format_to(it, "-");
					else
						it = format_to(it, " - ");
					c = -c;
				}
				else if (!first)
					it = format_to(it, " + ");
				first = false;
				if (i == 0)
					it = format_to(it, "{}", c);
				else if (!(c == 1))
					it = format_to(it, "{}*", c);

				if (i == 1)
					it = format_to(it, "x");
				else if (i >= 2)
					it = format_to(it, "x^{}", i);
			}
		}
		else
		{
			bool first = true;
			for (int i = 0; i <= poly.degree(); ++i)
			{
				auto c = poly[i];
				if (c == 0)
					continue;

				if (is_negative(c))
				{
					if (first)
						it = format_to(it, "-");
					else
						it = format_to(it, " - ");
					c = -c;
				}
				else if (!first)
					it = format_to(it, " + ");
				first = false;

				if (i == 0)
					it = format_to(it, "({})", c);
				else if (!(c == 1))
					it = format_to(it, "({})*", c);

				if (i == 1)
					it = format_to(it, "x");
				else if (i >= 2)
					it = format_to(it, "x^{}", i);
			}

			it = format_to(it, " + O(x^{})", poly.max_order() + 1);
		}

		return it;
	}
};

#endif
