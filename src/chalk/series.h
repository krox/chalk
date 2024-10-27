#pragma once

#include "chalk/rings.h"
#include "fmt/format.h"
#include <array>
#include <cassert>
#include <climits>
#include <utility>
#include <vector>

namespace chalk {

// some forward declarations to make the compiler find everything...
template <typename R, int N, char X> class Series;
template <typename R, int N, char X, typename F>
auto mapCoefficients(F f, Series<R, N, X> const &a)
    -> Series<decltype(f(a[0])), N, X>;

/**
 * Truncated series.
 *   - order is compile-time
 *   - can be used for automatic differentiation (for a single variable)
 */
template <typename R, int N, char X = 'x'> class Series
{
	static_assert(N >= 1);

  public:
	// The coeffs dont fulfill any invariants, so they can be public.
	// But for consistency with other types, we keep the naming scheme.
	std::array<R, N + 1> coefficients_;

  public:
	Series() = default;
	explicit Series(int c)
	{
		coefficients_.fill(R(0));
		coefficients_[0] = R(c);
	}
	explicit Series(const R &c)
	{
		coefficients_.fill(R(0));
		coefficients_[0] = c;
	}

	static Series generator()
	{
		Series r;
		r.coefficients_.fill(R(0));
		r.coefficients_[1] = R(1);
		return r;
	}

	std::array<R, N + 1> &coefficients() { return coefficients_; }
	std::array<R, N + 1> const &coefficients() const { return coefficients_; }

	/** array-like access to coefficients */
	R &operator[](int i) { return coefficients_.at(i); }
	R const &operator[](int i) const { return coefficients_.at(i); }

	/** s.mapCoefficients(fun, args...) = mapCoefficients(fun(_, args...), s) */
	template <typename Fun, typename... Args>
	auto mapCoefficients(Fun fun, Args const &...args)
	{
		return chalk::mapCoefficients(
		    [&fun, &args...](R const &c) { return fun(c, args...); }, *this);
	}
};

template <typename R, int N, char X>
bool operator==(Series<R, N, X> const &a, Series<R, N, X> const &b)
{
	return a.coefficients() == b.coefficients();
}

template <typename R, int N, char X, typename U>
bool operator==(Series<R, N, X> const &a, Scalar<U> const &b)
{
	for (int i = 1; i <= N; ++i)
		if (!is_zero(a[i]))
			return false;
	return a[0] == b.value;
}

template <typename R, int N, char X>
bool operator==(Series<R, N, X> const &a, R const &b)
{
	for (int i = 1; i <= N; ++i)
		if (!is_zero(a[i]))
			return false;
	return a[0] == b;
}

template <typename R, int N, char X>
bool operator==(Series<R, N, X> const &a, int b)
{
	for (int i = 1; i <= N; ++i)
		if (!is_zero(a[i]))
			return false;
	return a[0] == b;
}

/** unary operations */
template <typename R, int N, char X>
Series<R, N, X> operator-(Series<R, N, X> const &a)
{
	Series<R, N, X> r;
	for (int i = 0; i < N; ++i)
		r[i] = -a[i];
	return r;
}

/** binary operations (Series <-> Series) */
template <typename R, int N, typename S, char X>
auto operator+(Series<R, N, X> const &a, Series<S, N, X> const &b)
{
	using T = decltype(a[0] + b[0]);
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = a[i] + b[i];
	return r;
}
template <typename R, int N, typename S, char X>
auto operator-(Series<R, N, X> const &a, Series<S, N, X> const &b)
{
	using T = decltype(a[0] - b[0]);
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = a[i] - b[i];
	return r;
}
template <typename R, int N, typename S, char X>
auto operator*(Series<R, N, X> const &a, Series<S, N, X> const &b)
{
	using T = decltype(a[0] * b[0]);
	auto r = Series<T, N, X>(0);
	for (int i = 0; i <= N; ++i)
		for (int j = 0; i + j <= N; ++j)
			r[i + j] += a[i] * b[j];
	return r;
}

/** binary operations (Series <-> scalar) */
template <typename R, int N, char X, typename U>
Series<R, N, X> operator+(Series<R, N, X> const &a, Scalar<U> const &b)
{
	auto r = a;
	r[0] += b.value;
	return r;
}
template <typename R, int N, char X, typename U>
Series<R, N, X> operator-(Series<R, N, X> const &a, Scalar<U> const &b)
{
	auto r = a;
	r[0] -= b.value;
	return r;
}
template <typename R, int N, char X, typename U>
auto operator*(Series<R, N, X> const &a, Scalar<U> const &b)
{
	using T = decltype(a[0] * b.value);
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = a[i] * b.value;
	return r;
}
template <typename R, int N, char X, typename U>
Series<R, N, X> operator/(Series<R, N, X> const &a, Scalar<U> const &b)
{
	using T = decltype(a[0] / b.value);
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = a[i] / b.value;
	return r;
}
template <typename R, int N, char X>
auto operator*(Series<R, N, X> const &a, R b)
{
	using T = decltype(a[0] * b);
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = a[i] * b;
	return r;
}
template <typename R, int N, char X>
auto operator*(R a, Series<R, N, X> const &b)
{
	using T = decltype(a * b[0]);
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = a * b[i];
	return r;
}
template <typename R, int N, char X>
Series<R, N, X> operator/(Series<R, N, X> const &a, R b)
{
	using T = decltype(a[0] / b);
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = a[i] / b;
	return r;
}
template <typename R, int N, char X>
auto operator*(Series<R, N, X> const &a, int b)
{
	using T = decltype(a[0] * b);
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = a[i] * b;
	return r;
}
template <typename R, int N, char X>
Series<R, N, X> operator/(Series<R, N, X> const &a, int b)
{
	using T = decltype(a[0] / b);
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = a[i] / b;
	return r;
}

/** op-assigns for convenience and speed */
template <typename R, int N, typename S, char X>
Series<R, N, X> &operator+=(Series<R, N, X> &a, Series<S, N, X> const &b)
{
	for (int i = 0; i <= N; ++i)
		a[i] += b[i];
	return a;
}
template <typename R, int N, typename S, char X>
Series<R, N, X> &operator-=(Series<R, N, X> &a, Series<S, N, X> const &b)
{
	for (int i = 0; i <= N; ++i)
		a[i] -= b[i];
	return a;
}
template <typename R, int N, typename S, char X>
Series<R, N, X> &operator*=(Series<R, N, X> &a, Series<S, N, X> const &b)
{
	a = a * b;
	return a;
}

/** op-assigns with scalar */
template <typename R, int N, char X, typename U>
Series<R, N, X> &operator+=(Series<R, N, X> &a, Scalar<U> const b)
{
	a[0] += b.value;
	return a;
}
template <typename R, int N, char X, typename U>
Series<R, N, X> &operator-=(Series<R, N, X> &a, Scalar<U> const b)
{
	a[0] -= b.value;
	return a;
}
template <typename R, int N, char X, typename U>
Series<R, N, X> &operator*=(Series<R, N, X> &a, Scalar<U> const &b)
{
	for (int i = 0; i <= N; ++i)
		a[i] *= b.value;
	return a;
}
template <typename R, int N, char X, typename U>
Series<R, N, X> &operator/=(Series<R, N, X> &a, Scalar<U> const &b)
{
	for (int i = 0; i <= N; ++i)
		a[i] /= b.value;
	return a;
}
template <typename R, int N, char X>
Series<R, N, X> &operator*=(Series<R, N, X> &a, R b)
{
	for (int i = 0; i <= N; ++i)
		a[i] *= b;
	return a;
}
template <typename R, int N, char X>
Series<R, N, X> &operator/=(Series<R, N, X> &a, R b)
{
	for (int i = 0; i <= N; ++i)
		a[i] /= b;
	return a;
}
template <typename R, int N, char X>
Series<R, N, X> &operator*=(Series<R, N, X> &a, int b)
{
	for (int i = 0; i <= N; ++i)
		a[i] *= b;
	return a;
}
template <typename R, int N, char X>
Series<R, N, X> &operator/=(Series<R, N, X> &a, int b)
{
	for (int i = 0; i <= N; ++i)
		a[i] /= b;
	return a;
}

template <typename R, int N, char X, typename F>
auto mapCoefficients(F f, Series<R, N, X> const &a)
    -> Series<decltype(f(a[0])), N, X>
{
	using T = decltype(f(a[0]));
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = f(a[i]);
	return r;
}

template <typename R, int N, char X, typename F>
auto mapCoefficientsNested(F f, Series<R, N, X> const &a)
    -> Series<decltype(mapCoefficients(f, a[0])), N, X>
{
	using T = decltype(mapCoefficients(f, a[0]));
	Series<T, N, X> r;
	for (int i = 0; i <= N; ++i)
		r[i] = mapCoefficients(f, a[i]);
	return r;
}

template <typename R, int N, char X> struct RingTraits<Series<R, N, X>>
{
	static bool is_zero(Series<R, N, X> const &s)
	{
		for (R const &c : s.coefficients())
			if (!is_zero(c))
				return false;
		return true;
	}
	static bool is_one(Series<R, N, X> const &s)
	{
		for (int i = 1; i <= N; ++i)
			if (!is_zero(s[i]))
				return false;
		return is_one(s[0]);
	}

	static bool is_negative(Series<R, N, X> const &) { return false; }
	static bool need_parens_product(Series<R, N, X> const &) { return true; }
	static bool need_parens_power(Series<R, N, X> const &) { return true; }
};

} // namespace chalk

template <typename R, int N, char X>
struct fmt::formatter<chalk::Series<R, N, X>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Series<R, N, X> &s,
	            FormatContext &ctx) const -> decltype(ctx.out())
	{
		auto it = ctx.out();

		bool first = true;
		bool printed_someting = false;
		for (int i = 0; i <= N; ++i)
		{
			auto c = s[i];
			if (c == 0)
				continue;
			printed_someting = true;

			if (is_negative(c))
			{
				if (first)
					it = fmt::format_to(it, "-");
				else
					it = fmt::format_to(it, " - ");
				c = -c;
			}
			else if (!first)
				it = fmt::format_to(it, " + ");
			first = false;

			if (i == 0)
				it = fmt::format_to(it, "({})", c);
			else if (!(c == 1))
				it = fmt::format_to(it, "({})*", c);

			if (i == 1)
				it = fmt::format_to(it, "{}", X);
			else if (i >= 2)
				it = fmt::format_to(it, "{}^{}", X, i);
		}

		if (!printed_someting)
			it = fmt::format_to(it, "0");

		it = fmt::format_to(it, " + O({}^{})", X, N + 1);

		return it;
	}
};
