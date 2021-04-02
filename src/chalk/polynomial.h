#pragma once

#include "chalk/rings.h"
#include "fmt/format.h"
#include <array>
#include <cassert>
#include <climits>
#include <utility>
#include <vector>

namespace chalk {

/**
 * Univariate polynomial with dense storage and coefficients in R.
 *   - can also be used truncated power series.
 *   - mixed-type operations are supported. for example
 *       Poly<int> * Poly<float> -> Poly<float>
 *       Poly<int> * float -> Poly<float>
 *   - in order to get things like
 *       Poly<T,x> * Poly<T,y> -> Poly<Poly<T,x>,y>
 *     some more work is needed. for example
 *       Poly<Poly<T,x>,y>(...) * Poly<T,y>
 *     or in some future maybe
 *       Scalar(Poly<T,x>) * Poly<T,y>
 */
template <typename R, char X = 'x'> class Polynomial
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
	explicit Polynomial(int c) : coefficients_{R(c)} { cleanup(); }
	explicit Polynomial(const R &c) : coefficients_{c} { cleanup(); }
	explicit Polynomial(std::vector<R> coefficients, int max_order = INT_MAX)
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
		R y = R(0);
		for (int i = degree(); i >= 0; --i)
			y = y * x + coefficients_[i];
		return y;
	}

	// these need private access for efficient implementations */
	template <typename R2, typename S, char X2>
	friend Polynomial<R2, X2> &operator+=(Polynomial<R2, X2> &a,
	                                      Polynomial<S, X2> const &b);
	template <typename R2, typename S, char X2>
	friend Polynomial<R2, X2> &operator-=(Polynomial<R2, X2> &a,
	                                      Polynomial<S, X2> const &b);
	template <typename R2, typename S, char X2>
	friend Polynomial<R2, X2> &operator*=(Polynomial<R2, X2> &a,
	                                      Polynomial<S, X2> const &b);
	template <typename R2, char X2, typename U>
	friend Polynomial<R2, X2> &operator+=(Polynomial<R2, X2> &a, U const b);
	template <typename R2, char X2, typename U>
	friend Polynomial<R2, X2> &operator-=(Polynomial<R2, X2> &a, U const b);
	template <typename R2, char X2, typename U>
	friend Polynomial<R2, X2> &operator*=(Polynomial<R2, X2> &a, U const b);
	template <typename R2, char X2, typename U>
	friend Polynomial<R2, X2> &operator/=(Polynomial<R2, X2> &a, U const b);
};

template <typename R, char X>
Polynomial<R, X> truncate(Polynomial<R, X> const &poly, int order)
{
	assert(order >= 0);
	if (order >= poly.max_order())
		return poly;
	else
		return Polynomial<R, X>(poly.coefficients(), order);
}

template <typename R, char X>
bool operator==(Polynomial<R, X> const &a, Polynomial<R, X> const &b)
{
	// assert(a.max_order() == INT_MAX);
	// assert(b.max_order() == INT_MAX);
	return a.coefficients() == b.coefficients();
}

template <typename R, char X, typename U>
bool operator==(Polynomial<R, X> const &a, U const &b)
{
	// assert(a.max_order() == INT_MAX);
	if (is_zero(b))
		return is_zero(a);
	return a.degree() == 0 && a[0] == b;
}

/** unary operations */
template <typename R, char X>
Polynomial<R, X> operator-(Polynomial<R, X> const &a)
{
	std::vector<R> r;
	r.reserve(a.degree() + 1);
	for (auto &coeff : a.coefficients())
		r.emplace_back(-coeff);
	return Polynomial<R, X>(std::move(r), a.max_order());
}

/** binary operations (Polynomial <-> Polynomial) */
template <typename R, typename S, char X>
auto operator+(Polynomial<R, X> const &a, Polynomial<S, X> const &b)
{
	using T = decltype(a[0] + b[0]);
	std::vector<T> r;
	int deg = std::max(a.degree(), b.degree());
	r.reserve(deg + 1);
	for (int i = 0; i <= deg; ++i)
		r.emplace_back(a[i] + b[i]);
	return Polynomial<T, X>(std::move(r),
	                        std::min(a.max_order(), b.max_order()));
}
template <typename R, typename S, char X>
auto operator-(Polynomial<R, X> const &a, Polynomial<S, X> const &b)
{
	using T = decltype(a[0] - b[0]);
	std::vector<T> r;
	int deg = std::max(a.degree(), b.degree());
	r.reserve(deg + 1);
	for (int i = 0; i <= deg; ++i)
		r.emplace_back(a[i] - b[i]);
	return Polynomial<T, X>(std::move(r),
	                        std::min(a.max_order(), b.max_order()));
}
template <typename R, typename S, char X>
auto operator*(Polynomial<R, X> const &a, Polynomial<S, X> const &b)
{
	using T = decltype(a[0] * b[0]);
	auto r = std::vector<T>(std::max(0, a.degree() + b.degree() + 1), T(0));
	for (int i = 0; i <= a.degree(); ++i)
		for (int j = 0; j <= b.degree(); ++j)
			if (i + j <= std::min(a.max_order(), b.max_order()))
				r[i + j] += a[i] * b[j];
	return Polynomial<T, X>(std::move(r),
	                        std::min(a.max_order(), b.max_order()));
}

/** binary operations (Polynomial <-> scalar) */
template <typename R, char X, typename U>
Polynomial<R, X> operator+(Polynomial<R, X> const &a, U const &b)
{
	auto coeffs = a.coefficients();
	if (coeffs.empty())
		coeffs.emplace_back(b);
	else
		coeffs[0] += b;
	return Polynomial<R, X>(std::move(coeffs), a.max_order());
}
template <typename R, char X, typename U>
Polynomial<R, X> operator-(Polynomial<R, X> const &a, U const &b)
{
	auto coeffs = a.coefficients();
	if (coeffs.empty())
		coeffs.emplace_back(-b);
	else
		coeffs[0] -= b;
	return Polynomial<R, X>(std::move(coeffs), a.max_order());
}
template <typename R, char X, typename U>
auto operator*(Polynomial<R, X> const &a, U const &b)
{
	using T = decltype(a[0] * b);
	std::vector<T> r;
	r.reserve(a.degree() + 1);
	for (int i = 0; i <= a.degree(); ++i)
		r.emplace_back(a[i] * b);
	return Polynomial<T, X>(std::move(r), a.max_order());
}
template <typename R, char X, typename U>
Polynomial<R, X> operator/(Polynomial<R, X> const &a, U const &b)
{
	using T = decltype(a[0] / b);
	std::vector<T> r;
	r.reserve(a.degree() + 1);
	for (int i = 0; i <= a.degree(); ++i)
		r.emplace_back(a[i] / b);
	return Polynomial<T, X>(std::move(r), a.max_order());
}
template <typename U, typename R, char X>
Polynomial<R, X> operator+(U const &a, Polynomial<R, X> const &b)
{
	return b + a;
}
template <typename U, typename R, char X>
Polynomial<R, X> operator-(U const &a, Polynomial<R, X> const &b)
{
	auto r = -b;
	r += a;
	return r;
}
template <typename U, typename R, char X>
auto operator*(U const &a, Polynomial<R, X> const &b)
{
	using T = decltype(a * b[0]);
	std::vector<T> r;
	r.reserve(b.degree() + 1);
	for (int i = 0; i <= b.degree(); ++i)
		r.emplace_back(a * b[i]);
	return Polynomial<T, X>(std::move(r), b.max_order());
}

/** op-assigns for convenience and speed */
template <typename R, typename S, char X>
Polynomial<R, X> &operator+=(Polynomial<R, X> &a, Polynomial<S, X> const &b)
{
	if (a.degree() < b.degree())
		a.coefficients_.resize(b.degree() + 1, R(0));
	for (int i = 0; i <= b.degree(); ++i)
		a.coefficients_[i] += b[i];
	a.max_order_ = std::min(a.max_order(), b.max_order());
	a.cleanup();
	return a;
}
template <typename R, typename S, char X>
Polynomial<R, X> &operator-=(Polynomial<R, X> &a, Polynomial<S, X> const &b)
{
	if (a.degree() < b.degree())
		a.coefficients_.resize(b.degree() + 1, R(0));
	for (int i = 0; i <= b.degree(); ++i)
		a.coefficients_[i] -= b[i];
	a.max_order_ = std::min(a.max_order(), b.max_order());
	a.cleanup();
	return a;
}
template <typename R, typename S, char X>
Polynomial<R, X> &operator*=(Polynomial<R, X> &a, Polynomial<S, X> const &b)
{
	a = a * b;
	return a;
}
template <typename R, char X, typename U>
Polynomial<R, X> &operator+=(Polynomial<R, X> &a, U const b)
{
	if (a.degree() == -1)
		a.coefficients_.emplace_back(b);
	else
		a.coefficients_[0] += b;
	a.cleanup();
	return a;
}
template <typename R, char X, typename U>
Polynomial<R, X> &operator-=(Polynomial<R, X> &a, U const b)
{
	if (a.degree() == -1)
		a.coefficients_.emplace_back(-b);
	else
		a.coefficients_[0] -= b;
	a.cleanup();
	return a;
}
template <typename R, char X, typename U>
Polynomial<R, X> &operator*=(Polynomial<R, X> &a, U const b)
{
	for (R &coeff : a.coefficients_)
		coeff *= b;
	a.cleanup();
	return a;
}
template <typename R, char X, typename U>
Polynomial<R, X> &operator/=(Polynomial<R, X> &a, U const b)
{
	for (R &coeff : a.coefficients_)
		coeff /= b;
	a.cleanup();
	return a;
}

template <typename R, char X, typename F>
auto map_coefficients(F f, Polynomial<R, X> const &a)
    -> Polynomial<decltype(f(a[0]))>
{
	std::vector<decltype(f(a[0]))> coeffs;
	coeffs.reserve(a.coefficients().size());
	for (auto &coeff : a.coefficients())
		coeffs.emplace_back(f(coeff));
	return Polynomial<decltype(f(a[0]))>(std::move(coeffs), a.max_order());
}

template <typename R, char X> struct RingTraits<Polynomial<R, X>>
{
	static bool is_zero(Polynomial<R, X> const &poly)
	{
		return poly.coefficients().empty();
	}
	static bool is_one(Polynomial<R, X> const &poly)
	{
		return poly.coefficients().size() == 1 &&
		       chalk::is_one(poly.coefficients()[0]);
	}
	static bool is_negative(Polynomial<R, X> const &poly)
	{
		return poly.degree() >= 0 &&
		       chalk::is_negative(poly.coefficients().back());
	}

	/** not exactly correct... */
	static bool need_parens_product(Polynomial<R, X> const &) { return true; }
	static bool need_parens_power(Polynomial<R, X> const &) { return true; }
};

} // namespace chalk

template <typename R, char X> struct fmt::formatter<chalk::Polynomial<R, X>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Polynomial<R, X> &poly, FormatContext &ctx)
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
					it = format_to(it, "{}", X);
				else if (i >= 2)
					it = format_to(it, "{}^{}", X, i);
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
					it = format_to(it, "{}", X);
				else if (i >= 2)
					it = format_to(it, "{}^{}", X, i);
			}

			it = format_to(it, " + O({}^{})", X, poly.max_order() + 1);
		}

		return it;
	}
};
