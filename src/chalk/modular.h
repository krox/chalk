#ifndef CHALK_MODULAR_H
#define CHALK_MODULAR_H

/**
 * Basic modular arithmetic for 64-bit (and some 128-bit) integers.
 */

#include "fmt/format.h"
#include <algorithm>
#include <cassert>
#include <cstdint>

namespace chalk {

using int128_t = __int128_t;

/**
 * basic arithmetic operations
 *   - overflow safe
 *   - plenty of asserts
 *   - might be faster than naive expressions
 *     (turns out a branch is faster than a '%' at times)
 */
inline int64_t negmod(int64_t a, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if (a == 0)
		return 0;
	else
		return m - a;
}

inline int64_t invmod(int64_t a, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	int64_t a0 = m;
	int64_t a1 = a;
	int64_t b0 = 0;
	int64_t b1 = 1;

	while (a1 > 1)
	{
		int64_t q = a0 / a1;
		int64_t a2 = a0 - q * a1;
		int64_t b2 = b0 - q * b1;

		a0 = a1;
		a1 = a2;
		b0 = b1;
		b1 = b2;
	}

	if (b1 < 0)
		b1 += m;
	assert(0 <= b1 && b1 < m);
	return b1;
}

inline int64_t addmod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (b < m - a)
		return a + b;
	else
		return a + b - m;
}

inline int64_t submod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (a >= b)
		return a - b;
	else
		return a - b + m;
}

inline int64_t mulmod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	__int128_t a_ = a;
	__int128_t b_ = b;
	__int128_t m_ = m;
	return (int64_t)((a_ * b_) % m_);
}

inline int64_t divmod(int64_t a, int64_t b, int64_t m)
{
	b = invmod(b, m);
	return mulmod(a, b, m);
}

inline int64_t powmod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if (b < 0)
	{
		a = invmod(a, m);
		b = -b;
	}

	int64_t r = 1;
	for (; b != 0; b >>= 1, a = mulmod(a, a, m))
		if (b & 1)
			r = mulmod(r, a, m);
	return r;
}

inline __int128_t negmod(__int128_t a, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if (a == 0)
		return 0;
	else
		return m - a;
}

inline __int128_t invmod(__int128_t a, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	__int128_t a0 = m;
	__int128_t a1 = a;
	__int128_t b0 = 0;
	__int128_t b1 = 1;

	while (a1 > 1)
	{
		__int128_t q = a0 / a1;
		__int128_t a2 = a0 - q * a1;
		__int128_t b2 = b0 - q * b1;

		a0 = a1;
		a1 = a2;
		b0 = b1;
		b1 = b2;
	}

	if (b1 < 0)
		b1 += m;
	assert(0 <= b1 && b1 < m);
	return b1;
}

inline __int128_t addmod(__int128_t a, __int128_t b, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (b < m - a)
		return a + b;
	else
		return a + b - m;
}

inline __int128_t submod(__int128_t a, __int128_t b, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (a >= b)
		return a - b;
	else
		return a - b + m;
}

inline __int128_t mulmod(__int128_t a, __int128_t b, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (a > b)
		std::swap(a, b);

	__int128_t r = 0;
	for (; a != 0; a >>= 1, b = addmod(b, b, m))
		if (a & 1)
			r = addmod(r, b, m);
	return r;
}

inline __int128_t divmod(__int128_t a, __int128_t b, __int128_t m)
{
	b = invmod(b, m);
	return mulmod(a, b, m);
}

inline __int128_t powmod(__int128_t a, __int128_t b, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if (b < 0)
	{
		a = invmod(a, m);
		b = -b;
	}

	__int128_t r = 1;
	for (; b != 0; b >>= 1, a = mulmod(a, a, m))
		if (b & 1)
			r = mulmod(r, a, m);
	return r;
}

/**
 * Greatest common divisor and least common multiple.
 * Sign of a and b is ignored, result is always >= 0.
 * Convention: gcd(0,x) = gcd(x,0) = abs(x)
 *             lcm(0,x) = lcm(x,0) = 0
 */
int64_t gcd(int64_t a, int64_t b);
int64_t lcm(int64_t a, int64_t b);

/**
 * Residual classes [a] = a + mℤ with a,m ∈ ℤ, m > 0
 *   - is a field iff m is prime
 *   - might be refactored into general quotient-ring-class
 */
class IntMod
{
	int64_t a_ = 0, m_ = 0; // m_=0 corresponds to NaN state

  public:
	/** constructors */
	IntMod() = default;
	explicit IntMod(int64_t a, int64_t m) : a_(a), m_(m)
	{
		assert(m_ > 0);
		assert(0 <= a_ && a_ < m_);
	}

	int64_t a() const { return a_; }
	int64_t m() const { return m_; }

	static IntMod make(int64_t a, int64_t b)
	{
		if (b < 0)
			b = -b;
		assert(b != 0);
		a %= b;
		if (a < 0)
			a += b;
		return IntMod(a, b);
	}
};

/** properties of elements */
inline bool invertible(IntMod a) { return gcd(a.a(), a.m()) == 1; }
inline bool regular(IntMod a) { return gcd(a.a(), a.m()) == 1; }

/** unary operations */
inline IntMod operator-(IntMod a)
{
	return IntMod(negmod(a.a(), a.m()), a.m());
}
inline IntMod inverse(IntMod a) { return IntMod(invmod(a.a(), a.m()), a.m()); }

/** binary operations */
inline IntMod operator+(IntMod a, IntMod b)
{
	assert(a.m() == b.m());
	return IntMod(addmod(a.a(), b.a(), a.m()), a.m());
}
inline IntMod operator-(IntMod a, IntMod b)
{
	assert(a.m() == b.m());
	return IntMod(submod(a.a(), b.a(), a.m()), a.m());
}
inline IntMod operator*(IntMod a, IntMod b)
{
	assert(a.m() == b.m());
	return IntMod(mulmod(a.a(), b.a(), a.m()), a.m());
}
inline IntMod operator/(IntMod a, IntMod b)
{
	assert(a.m() == b.m());
	return IntMod(divmod(a.a(), b.a(), a.m()), a.m());
}
inline IntMod pow(IntMod a, int b)
{
	return IntMod(powmod(a.a(), b, a.m()), a.m());
}

/** binary assigns */
inline void operator+=(IntMod &a, IntMod b) { a = a + b; }
inline void operator-=(IntMod &a, IntMod b) { a = a - b; }
inline void operator*=(IntMod &a, IntMod b) { a = a * b; }
inline void operator/=(IntMod &a, IntMod b) { a = a / b; }

/** comparisons */
inline bool operator==(IntMod a, IntMod b)
{
	assert(a.m() == b.m());
	return a.a() == b.a();
}
inline bool operator!=(IntMod a, IntMod b)
{
	assert(a.m() == b.m());
	return a.a() == b.a();
}
inline bool operator==(IntMod a, int b) { return a == IntMod::make(b, a.m()); }
inline bool operator!=(IntMod a, int b) { return a != IntMod::make(b, a.m()); }

} // namespace chalk

namespace fmt {

template <> struct formatter<chalk::IntMod>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::IntMod &a, FormatContext &ctx)
	    -> decltype(ctx.out())
	{
		return format_to(ctx.out(), "[{}]", a.a());
	}
};

} // namespace fmt

#endif
