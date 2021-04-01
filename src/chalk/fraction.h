#pragma once

#include "chalk/integer.h"
#include "chalk/numtheory.h"
#include "chalk/rings.h"
#include "fmt/format.h"

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
		if (is_negative(denom_))
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

/** comparisons */
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
template <typename R> bool operator==(Fraction<R> const &a, int b)
{
	return a.denom() == 1 && a.num() == b;
}
template <typename R> bool operator!=(Fraction<R> const &a, int b)
{
	return !(a.denom() == 1 && a.num() == b);
}

template <typename R> struct RingTraits<Fraction<R>>
{
	static bool is_zero(Fraction<R> const &value)
	{
		return chalk::is_zero(value.num());
	}
	static bool is_one(Fraction<R> const &value)
	{
		return chalk::is_one(value.num()) && chalk::is_one(value.denom());
	}
	static bool is_negative(Fraction<R> const &value)
	{
		return chalk::is_negative(value.num());
	}

	static bool need_parens_product(Fraction<R> const &value)
	{
		if (chalk::is_one(value.denom()))
			return chalk::need_parens_product(value.num());
		else
			return false;
	}
	static bool need_parens_power(Fraction<R> const &value)
	{
		if (chalk::is_one(value.denom()))
			return chalk::need_parens_power(value.num());
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
		if (chalk::is_one(x.denom()))
			return format_to(ctx.out(), "{}", x.num());

		if (chalk::need_parens_product(x.num()))
		{
			if (chalk::need_parens_power(x.denom()))
				return format_to(ctx.out(), "({})/({})", x.num(), x.denom());
			else
				return format_to(ctx.out(), "({})/{}", x.num(), x.denom());
		}
		else
		{
			if (chalk::need_parens_power(x.denom()))
				return format_to(ctx.out(), "{}/({})", x.num(), x.denom());
			else
				return format_to(ctx.out(), "{}/{}", x.num(), x.denom());
		}
	}
};
