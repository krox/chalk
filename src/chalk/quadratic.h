#pragma once

#include "chalk/rings.h"
#include "fmt/format.h"
#include <cassert>
#include <climits>
#include <optional>

namespace chalk {

// numbers of the form a + b√d with a,b,d ∈ R
//   - R assumed to be abelian (otherwise, inverse() is broken)
//   - d should be a non-square, we further normalize it to be be 'square-free'
//     (though depending on R, it might be tricky to define what that means)
//   - TODO: this implementation is quite inefficient. For example:
//         - many of the (simple) operators could actually use an 'unchecked'
//           version of the constructor that skips the 'removeSquareFactor' step
//         - for R = rational numbers, it would be more efficient to store this
//           as four integers instead of three rationals
//         - the op-assign operators should actually operate inplace
template <ring_element R> class Quadratic
{
	R a_, b_, d_; // d is undefined if b=0

  public:
	// constructrs
	Quadratic() = default;
	explicit Quadratic(int a) : a_(R(a)), b_(R(0)) {}
	explicit Quadratic(R a) : a_(std::move(a)), b_(0) {}
	Quadratic(R a, R b, R d)
	    : a_(std::move(a)), b_(std::move(b)), d_(std::move(d))
	{
		b_ *= removeSquareFactor(d_);
		if (is_one(d_))
		{
			a_ += b_;
			b_ = R(0);
		}
	}

	// attributes
	R const &real() const { return a_; }
	R const &imag() const { return b_; }
	R const &d() const { return d_; }

	// conversion (to floating point)
	template <class T> explicit operator T() const
	{
		return (T)a_ + (T)b_ * sqrt((T)d_);
	}
};

// unary Quadratic

template <class R> Quadratic<R> operator-(Quadratic<R> const &a)
{
	return {-a.real(), -a.imag(), a.d()};
}
template <class R> Quadratic<R> conj(Quadratic<R> const &a)
{
	// 'Galois conjugate'. This might or might not be the same as 'complex conj'
	return {a.real(), -a.imag(), a.d()};
}
template <class R> R norm2(Quadratic<R> const &a)
{
	// in number theory, this is just called "norm", not "squared norm".
	// but our naming is more useful when using Quadratic<Rational> as an
	// exact representation of (a subset of) real(/complex) numbers.
	return a.real() * a.real() - a.d() * a.imag() * a.imag();
}
template <class R> Quadratic<R> inverse(Quadratic<R> const &a)
{
	return conj(a) * inverse(norm2(a));
}

template <class R> std::optional<Quadratic<R>> trySqrt(Quadratic<R> const &a)
{
	if (is_zero(a.imag()))
		return Quadratic<R>{R(0), R(1), a.real()};

	auto n = trySqrt(norm2(a));
	if (!n)
		return {};

	auto x = trySqrt((a.real() + *n) / R(2));
	if (!x)
		x = trySqrt((a.real() - *n) / R(2));
	if (!x)
		return {};
	auto y = a.imag() / (*x * R(2));
	return Quadratic<R>{*std::move(x), std::move(y), a.d()};
}

template <class R> Quadratic<R> sqrt(Quadratic<R> const &a)
{
	return trySqrt(a).value();
}

template <class R> Quadratic<R> pow(Quadratic<R> a, int e)
{
	if (e < 0)
	{
		assert(e != INT_MIN);
		e = -e;
		a = inverse(a);
	}

	// TODO: this could be optimized (remove a lot of implicit 'd' checks)
	auto r = Quadratic<R>(1);
	for (; e != 0; e >>= 1, a *= a)
		if (e & 1)
			r *= a;
	return r;
}

// binary Quadratic <-> scalar

template <class R> Quadratic<R> operator+(Quadratic<R> const &a, R const &b)
{
	return {a.real() + b, a.imag(), a.d()};
}
template <class R> Quadratic<R> operator-(Quadratic<R> const &a, R const &b)
{
	return {a.real() - b, a.imag(), a.d()};
}
template <class R> Quadratic<R> operator*(Quadratic<R> const &a, R const &b)
{
	return {a.real() * b, a.imag() * b, a.d()};
}
template <class R> Quadratic<R> operator/(Quadratic<R> const &a, R const &b)
{
	return {a.real() / b, a.imag() / b, a.d()};
}
template <class R> Quadratic<R> operator+(Quadratic<R> const &a, int b)
{
	return a + R(b);
}
template <class R> Quadratic<R> operator-(Quadratic<R> const &a, int b)
{
	return a - R(b);
}
template <class R> Quadratic<R> operator*(Quadratic<R> const &a, int b)
{
	return a * R(b);
}
template <class R> Quadratic<R> operator/(Quadratic<R> const &a, int b)
{
	return a / R(b);
}

// binary Quadratic <-> Quadratic

template <class R>
std::optional<Quadratic<R>> tryAdd(Quadratic<R> const &a, Quadratic<R> const &b)
{
	if (is_zero(a.imag()))
		return Quadratic<R>{a.real() + b.real(), b.imag(), b.d()};
	if (is_zero(b.imag()))
		return Quadratic<R>{a.real() + b.real(), a.imag(), a.d()};
	if (!(a.d() == b.d()))
		return std::nullopt;
	return Quadratic<R>{a.real() + b.real(), a.imag() + b.imag(), a.d()};
}
template <class R>
std::optional<Quadratic<R>> trySub(Quadratic<R> const &a, Quadratic<R> const &b)
{
	if (is_zero(a.imag()))
		return Quadratic<R>{a.real() - b.real(), -b.imag(), b.d()};
	if (is_zero(b.imag()))
		return Quadratic<R>{a.real() - b.real(), a.imag(), a.d()};
	if (!(a.d() == b.d()))
		return std::nullopt;
	return Quadratic<R>{a.real() - b.real(), a.imag() - b.imag(), a.d()};
}
template <class R>
std::optional<Quadratic<R>> tryMul(Quadratic<R> const &a, Quadratic<R> const &b)
{
	if (is_zero(a.imag()))
		return Quadratic<R>{a.real() * b.real(), a.real() * b.imag(), b.d()};
	if (is_zero(b.imag()))
		return Quadratic<R>{a.real() * b.real(), a.imag() * b.real(), a.d()};
	if (!(a.d() == b.d()))
		return std::nullopt;
	return Quadratic<R>{a.real() * b.real() + a.d() * a.imag() * b.imag(),
	                    a.real() * b.imag() + a.imag() * b.real(), a.d()};
}

template <class R>
std::optional<Quadratic<R>> tryDiv(Quadratic<R> const &a, Quadratic<R> const &b)
{
	return tryMul(a, inverse(b));
}

// NOTE: .value() in noexcept function will std::terminate. this is intentional
template <class R>
Quadratic<R> operator+(Quadratic<R> const &a, Quadratic<R> const &b) noexcept
{
	return tryAdd(a, b).value();
}
template <class R>
Quadratic<R> operator-(Quadratic<R> const &a, Quadratic<R> const &b) noexcept
{
	return trySub(a, b).value();
}
template <class R>
Quadratic<R> operator*(Quadratic<R> const &a, Quadratic<R> const &b) noexcept
{
	return tryMul(a, b).value();
}
template <class R>
Quadratic<R> operator/(Quadratic<R> const &a, Quadratic<R> const &b) noexcept
{
	return tryDiv(a, b).value();
}

// op-assigns

template <class R, class U>
Quadratic<R> &operator+=(Quadratic<R> &a, U const &b)
{
	a = a + b;
	return a;
}
template <class R, class U>
Quadratic<R> &operator-=(Quadratic<R> &a, U const &b)
{
	a = a - b;
	return a;
}
template <class R, class U>
Quadratic<R> &operator*=(Quadratic<R> &a, U const &b)
{
	a = a * b;
	return a;
}
template <class R, class U>
Quadratic<R> &operator/=(Quadratic<R> &a, U const &b)
{
	a = a / b;
	return a;
}

// comparisons

template <class R> bool operator==(Quadratic<R> const &a, Quadratic<R> const &b)
{
	return a.real() == b.real() && a.imag() == b.imag() &&
	       (a.d() == b.d() || is_zero(a.imag()));
}
template <class R> bool operator==(Quadratic<R> const &a, R const &b)
{
	return is_zero(a.imag()) && a.real() == b;
}
template <class R> bool operator==(Quadratic<R> const &a, int b)
{
	return a == R(b);
}
template <class R, class U> bool operator!=(Quadratic<R> const &a, U const &b)
{
	return !(a == b);
}

template <class R> bool is_zero(Quadratic<R> const &a)
{
	return is_zero(a.real()) && is_zero(a.imag());
}
template <class R> bool is_one(Quadratic<R> const &a)
{
	return is_one(a.real()) && is_zero(a.imag());
}
template <class R> bool is_negative(Quadratic<R> const &a)
{
	return is_zero(a.imag()) && is_negative(a.real());
}
template <class R> bool need_parens_product(Quadratic<R> const &)
{
	return true;
}
template <class R> bool need_parens_power(Quadratic<R> const &) { return true; }

} // namespace chalk

template <class R> struct fmt::formatter<chalk::Quadratic<R>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <class FormatContext>
	auto format(const chalk::Quadratic<R> &x, FormatContext &ctx) const
	{
		using chalk::is_negative;
		using chalk::is_zero;
		if (is_zero(x.imag()))
			return format_to(ctx.out(), "{}", x.real());
		else if (is_zero(x.real()))
			return format_to(ctx.out(), "{}√{}", x.imag(), x.d());
		else if (is_negative(x.imag()))
			return format_to(ctx.out(), "{} - {}√{}", x.real(), -x.imag(),
			                 x.d());
		else
			return format_to(ctx.out(), "{} + {}√{}", x.real(), x.imag(),
			                 x.d());
	}
};
