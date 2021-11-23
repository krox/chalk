#pragma once

#include "chalk/integer.h"
#include "chalk/numtheory.h"
#include "chalk/rings.h"
#include "fmt/format.h"
#include <optional>

namespace chalk {

/**
 * Ring of fractions a/b with a,b ∈ R
 *   - assumes R is commutative
 *   - if R is an integral domain (no zero-divisors), this is a field
 */
template <typename R> class Fraction
{
	R num_, denom_;

	void reduce()
	{
		R g = gcd(num_, denom_);
		num_ /= g;
		denom_ /= g;
		if (isNegative(denom_))
		{
			num_ = -num_;
			denom_ = -denom_;
		}

		// makes sure that Fraction<int128_t> does not overflow
		// assert(INT64_MIN < num_ && num_ < INT64_MAX);
		// assert(INT64_MIN < denom_ && denom_ < INT64_MAX);
	}

  public:
	/** constructors */
	Fraction() : num_(0), denom_(1) {}
	explicit Fraction(int num) : num_(num), denom_(1) {}
	explicit Fraction(R num) : num_(std::move(num)), denom_(1) {}
	explicit Fraction(R num, R denom)
	    : num_(std::move(num)), denom_(std::move(denom))
	{
		reduce();
	}

	/** field access */
	R const &num() const { return num_; }
	R const &denom() const { return denom_; }

	template <typename T> explicit operator T() const
	{
		return T(num()) / T(denom());
	}
};

template <typename R> Fraction<R> operator-(Fraction<R> const &a)
{
	return Fraction<R>{-a.num(), a.denom()};
}
template <typename R> Fraction<R> inverse(Fraction<R> const &a)
{
	// assert(is_regular(a.num()));
	return Fraction<R>{a.denom(), a.num()};
}

/** binary operators (Fraction <-> Fraction) */
template <typename R>
Fraction<R> operator+(Fraction<R> const &a, Fraction<R> const &b)
{
	R num = a.num() * b.denom() + b.num() * a.denom();
	R denom = a.denom() * b.denom();
	return Fraction<R>{num, denom};
}
template <typename R>
Fraction<R> operator-(Fraction<R> const &a, Fraction<R> const &b)
{
	R num = a.num() * b.denom() - b.num() * a.denom();
	R denom = a.denom() * b.denom();
	return Fraction<R>{num, denom};
}
template <typename R>
Fraction<R> operator*(Fraction<R> const &a, Fraction<R> const &b)
{
	R num = a.num() * b.num();
	R denom = a.denom() * b.denom();
	return Fraction<R>{num, denom};
}
template <typename R>
Fraction<R> operator/(Fraction<R> const &a, Fraction<R> const &b)
{
	// assert(is_regular(b.num()));
	R num = a.num() * b.denom();
	R denom = a.denom() * b.num();
	return Fraction<R>{num, denom};
}

/** binary operators (Fraction <-> int) */
template <typename R> Fraction<R> operator+(Fraction<R> const &a, int b)
{
	return Fraction<R>{a.num() + b * a.denom(), a.denom()};
}
template <typename R> Fraction<R> operator-(Fraction<R> const &a, int b)
{
	return Fraction<R>{a.num() - b * a.denom(), a.denom()};
}
template <typename R> Fraction<R> operator*(Fraction<R> const &a, int b)
{
	return Fraction<R>{a.num() * b, a.denom()};
}
template <typename R> Fraction<R> operator/(Fraction<R> const &a, int b)
{
	// assert(is_regular(R(b)));
	return Fraction<R>{a.num(), a.denom() * b};
}
template <typename R> Fraction<R> operator+(int a, Fraction<R> const &b)
{
	return b + a;
}
template <typename R> Fraction<R> operator-(int a, Fraction<R> const &b)
{
	return Fraction<R>{a * b.denom() - b.num(), b.denom()};
}
template <typename R> Fraction<R> operator*(int a, Fraction<R> const &b)
{
	return b * a;
}
template <typename R> Fraction<R> operator/(int a, Fraction<R> const &b)
{
	// assert(is_regular(b));
	return Fraction<R>{a * b.denom(), b.num()};
}

/** binary operators assignments */
template <typename R> void operator+=(Fraction<R> &a, Fraction<R> const &b)
{
	a = a + b;
}
template <typename R> void operator-=(Fraction<R> &a, Fraction<R> const &b)
{
	a = a - b;
}
template <typename R> void operator*=(Fraction<R> &a, Fraction<R> const &b)
{
	a = a * b;
}
template <typename R> void operator/=(Fraction<R> &a, Fraction<R> const &b)
{
	a = a / b;
}
template <typename R> void operator+=(Fraction<R> &a, R const &b) { a = a + b; }
template <typename R> void operator-=(Fraction<R> &a, R const &b) { a = a - b; }
template <typename R> void operator*=(Fraction<R> &a, R const &b) { a = a * b; }
template <typename R> void operator/=(Fraction<R> &a, R const &b) { a = a / b; }
template <typename R> void operator+=(Fraction<R> &a, int b) { a = a + b; }
template <typename R> void operator-=(Fraction<R> &a, int b) { a = a - b; }
template <typename R> void operator*=(Fraction<R> &a, int b) { a = a * b; }
template <typename R> void operator/=(Fraction<R> &a, int b) { a = a / b; }

template <typename R> Fraction<R> removeSquareFactor(Fraction<R> &a)
{
	auto b = a.num() * a.denom();
	auto rDenom = a.denom() * a.denom();
	auto rNum = removeSquareFactor(b);
	a = Fraction<R>(b);
	return Fraction<R>(std::move(rNum), std::move(rDenom));
}

template <typename R> std::optional<Fraction<R>> trySqrt(Fraction<R> const &a)
{
	if (auto num = trySqrt(a.num()); num)
		if (auto denom = trySqrt(a.denom()); denom)
			return Fraction<R>{*std::move(num), *std::move(denom)};
	return {};
}

// comparisons

template <typename R>
bool operator==(Fraction<R> const &a, Fraction<R> const &b)
{
	// return a.num() * b.denom() == a.denom() * b.num();

	// TODO: check for which class for rings R this optimization is true
	return a.num() == b.num() && a.denom() == b.denom();
}
template <typename R>
bool operator!=(Fraction<R> const &a, Fraction<R> const &b)
{
	return a.num() != b.num() || a.denom() != b.denom();
}
template <typename R> bool operator<(Fraction<R> const &a, Fraction<R> const &b)
{
	return a.num() * b.denom() < b.num() * a.denom();
}
template <typename R>
bool operator<=(Fraction<R> const &a, Fraction<R> const &b)
{
	return a.num() * b.denom() <= b.num() * a.denom();
}
template <typename R> bool operator==(Fraction<R> const &a, int b)
{
	return a.denom() == 1 && a.num() == b;
}
template <typename R> bool operator!=(Fraction<R> const &a, int b)
{
	return !(a.denom() == 1 && a.num() == b);
}
template <typename R> bool operator<(Fraction<R> const &a, int b)
{
	return a.num() < b * a.denom();
}
template <typename R> bool operator<=(Fraction<R> const &a, int b)
{
	return a.num() <= b * a.denom();
}

template <typename R> struct RingTraits<Fraction<R>>
{
	static bool isZero(Fraction<R> const &value)
	{
		return chalk::isZero(value.num());
	}
	static bool isOne(Fraction<R> const &value)
	{
		return chalk::isOne(value.num()) && chalk::isOne(value.denom());
	}
	static bool isNegative(Fraction<R> const &value)
	{
		return chalk::isNegative(value.num());
	}

	static bool needParensProduct(Fraction<R> const &value)
	{
		if (chalk::isOne(value.denom()))
			return chalk::needParensProduct(value.num());
		else
			return false;
	}
	static bool needParensPower(Fraction<R> const &value)
	{
		if (chalk::isOne(value.denom()))
			return chalk::needParensPower(value.num());
		else
			return true;
	}
};

/** main usecase of fractions are of course rational numbers */
using Rational = Fraction<Integer>;
using Rational64 = Fraction<int64_t>;
using Rational128 = Fraction<int128_t>;

} // namespace chalk

template <typename R> struct fmt::formatter<chalk::Fraction<R>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Fraction<R> &x, FormatContext &ctx)
	{
		if (chalk::isOne(x.denom()))
			return format_to(ctx.out(), "{}", x.num());

		if (chalk::needParensProduct(x.num()))
		{
			if (chalk::needParensPower(x.denom()))
				return format_to(ctx.out(), "({})/({})", x.num(), x.denom());
			else
				return format_to(ctx.out(), "({})/{}", x.num(), x.denom());
		}
		else
		{
			if (chalk::needParensPower(x.denom()))
				return format_to(ctx.out(), "{}/({})", x.num(), x.denom());
			else
				return format_to(ctx.out(), "{}/{}", x.num(), x.denom());
		}
	}
};
