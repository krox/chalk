#ifndef CHALK_FRACTION_H
#define CHALK_FRACTION_H

#include "chalk/numtheory.h"
#include "fmt/format.h"

namespace chalk {

/** helper functions to make Fraction<int64_t> work */
bool is_negative(int64_t a) { return a < 0; }
bool is_invertible(int64_t a) { return a == 1 || a == -1; }
bool is_regular(int64_t a) { return a != 0; }
int64_t invertible_factor(int64_t a) { return a < 0 ? -1 : 1; }
// int64_t gcd(int64_t,int64_t); // defined in numtheory

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
		auto u = invertible_factor(denom_);
		num_ /= u;
		denom_ /= u;
	}

  public:
	/** constructors */
	Fraction() : num_(0), denom_(1) {}
	explicit Fraction(int num) : num_(num), denom_(1) {}
	explicit Fraction(R const &num) : num_(num), denom_(1) {}
	explicit Fraction(R const &num, R const &denom) : num_(num), denom_(denom)
	{
		reduce();
	}

	/** field access */
	R const &num() const { return num_; }
	R const &denom() const { return denom_; }
};

/** unary operators */
template <typename R> bool is_negative(Fraction<R> const &a)
{
	return is_negative(a.num());
}
template <typename R> bool is_invertible(Fraction<R> const &a)
{
	return is_regular(a.num());
}
template <typename R> bool is_regular(Fraction<R> const &a)
{
	return is_regular(a.num());
}
template <typename R> Fraction<R> invertible_factor(Fraction<R> const &a)
{
	return Fraction<R>{regular_factor(a.num()), a.denom()};
}
template <typename R> Fraction<R> operator-(Fraction<R> const &a)
{
	return Fraction<R>{-a.num(), a.denom()};
}
template <typename R> Fraction<R> inverse(Fraction<R> const &a)
{
	assert(is_regular(a.num()));
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
	assert(regular(b.num()));
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
	assert(is_regular(b));
	return Fraction<R>{a.num(), a.denom() * b};
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

} // namespace chalk

template <typename R> struct fmt::formatter<chalk::Fraction<R>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Fraction<R> &x, FormatContext &ctx)
	{
		if (x.denom() == 1)
			return format_to(ctx.out(), "{}", x.num());
		else
			return format_to(ctx.out(), "{}/{}", x.num(), x.denom());
	}
};

#endif
