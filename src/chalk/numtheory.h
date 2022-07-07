#pragma once

// Modular arithmetic, some number theory and some combinatorics for "small"
// (i.e. non GMP) numbers. Mostly (signed) 64 bit.

#include "fmt/format.h"
#include "util/linalg.h"
#include <cassert>
#include <cstdint>
#include <vector>

namespace chalk {

using int128_t = __int128;
using uint128_t = unsigned __int128;

// basic arithmetic operations
//   - overflow safe
//   - plenty of asserts
//   - might be faster than naive expressions
//     (turns out a branch can be faster than a '%' at times)
inline constexpr int64_t negmod(int64_t a, int64_t m) noexcept
{
	assert(m > 0);
	assert(0 <= a && a < m);
	return a ? m - a : 0;
}

inline constexpr int64_t invmod(int64_t a, int64_t m) noexcept
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

inline constexpr int64_t addmod(int64_t a, int64_t b, int64_t m) noexcept
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);
	return b < m - a ? a + b : a + (b - m);
}

inline constexpr int64_t submod(int64_t a, int64_t b, int64_t m) noexcept
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);
	return a >= b ? a - b : a - b + m;
}

inline constexpr int64_t mulmod(int64_t a, int64_t b, int64_t m) noexcept
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);
	return (int128_t)a * (int128_t)b % m;
}

inline constexpr int64_t fmamod(int64_t a, int64_t b, int64_t c,
                                int64_t m) noexcept
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);
	assert(0 <= c && c < m);
	return ((int128_t)a * b + c) % m;
}

inline constexpr int64_t fmmamod(int64_t a, int64_t b, int64_t c, int64_t d,
                                 int64_t m) noexcept
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);
	assert(0 <= c && c < m);
	assert(0 <= d && d < m);
	return ((int128_t)a * b + (int128_t)c * d) % m;
}

inline constexpr int64_t divmod(int64_t a, int64_t b, int64_t m) noexcept
{
	return mulmod(a, invmod(b, m), m);
}

inline constexpr int64_t powmod(int64_t a, int64_t b, int64_t m) noexcept
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

// Greatest common divisor and least common multiple.
// Sign of a and b is ignored, result is always >= 0.
// Convention: gcd(0,x) = gcd(x,0) = abs(x)
//             lcm(0,x) = lcm(x,0) = 0

inline constexpr int64_t gcd(int64_t a, int64_t b) noexcept
{
	while (true)
	{
		if (a == 0)
			return b >= 0 ? b : -b;
		b %= a;

		if (b == 0)
			return a >= 0 ? a : -a;
		a %= b;
	}
}

inline constexpr int64_t lcm(int64_t a, int64_t b) noexcept
{
	if (a == 0 || b == 0)
		return 0;
	else
		return a / gcd(a, b) * b;
}

inline constexpr void remove_common_factor(int64_t &a, int64_t &b) noexcept
{
	auto d = gcd(a, b);
	a /= d;
	b /= d;
	if (a < 0)
	{
		a = -a;
		b = -b;
	}
}

// returns (s, a-s*s) with s = floor(sqrt(a)), or s=0 for negative a
inline constexpr std::pair<int64_t, int64_t> sqrt_rem(int64_t a) noexcept
{
	if (a <= 0)
		return {0, a};

	// fp approximation + 1 newton step. careful testing required to avoid
	// off-by-1 errors (newton can oscillate +-1 instead of converging exactly)
	auto s = (int64_t)std::sqrt((double)a + 1);
	s = (s + a / s) / 2;
	return {s, a - s * s};
}

inline constexpr int64_t isqrt(int64_t a) noexcept { return sqrt_rem(a).first; }

// tests if a is the square of an integer
// (both 0 and 1 are considered square, same as GMP's definition)
inline constexpr bool is_square(int64_t a) noexcept
{
	// there are only 12 squares mod 64, so this cheap check
	// catches ~81% of all numbers
	if ((18302063659313855980u >> (a & 63)) & 1)
		return false;
	return sqrt_rem(a).second == 0;
}

inline constexpr int64_t ipow(int64_t a, int n) noexcept
{
	assert(n >= 0);
	int64_t r = 1;
	while (n--)
		r *= a;
	return r;
}

inline constexpr std::pair<int64_t, int64_t> root_rem(int64_t a, int n) noexcept
{
	assert(n >= 1);
	assert(a >= 0);
	if (n == 1)
		return {a, 0};
	if (n == 2)
		return sqrt_rem(a);
	if (a < n)
		return {0, a};

	int64_t x = a / n;
	while (true)
	{
		int64_t axn1 = a; // = a / ipow(x, n - 1), but that expression overflows
		for (int i = 0; i < n - 1; ++i)
			axn1 /= x;
		auto y = ((n - 1) * x + axn1) / n;
		if (y >= x)
			break;
		x = y;
	}
	return {x, a - ipow(x, n)};
}

inline constexpr int64_t iroot(int64_t a, int n) noexcept
{
	return root_rem(a, n).first;
}

// true if a can be written as b^e with e>1
// (0 and 1 are considered powers, just as in GMP, negatives not supported)
inline constexpr bool is_power(int64_t a) noexcept
{
	assert(a >= 0);

	// 0, 1, 4, 8, 16, ...
	if ((a & (a - 1)) == 0)
		return a != 2;
	if (is_square(a))
		return true;
	for (auto n : {3, 5, 7, 11, 13, 17})
		if (root_rem(a, n).second == 0)
			return true;

	switch (a)
	{
	case ipow(3LL, 19):
	case ipow(5LL, 19):
	case ipow(6LL, 19):
	case ipow(7LL, 19):
	case ipow(8LL, 19):
	case ipow(3LL, 23):
	case ipow(5LL, 23):
	case ipow(6LL, 23):
	case ipow(3LL, 29):
	case ipow(3LL, 31):
	case ipow(3LL, 37):
		return true;
	default:
		return false;
	}
}

// compute all primes up to n (inclusive)
std::vector<int64_t> primes(int64_t n);

// tests if n (>= 3 ) is a strong probable prime to base a (>= 0)
bool is_sprp(int64_t a, int64_t n);

// tests if n is prime number
bool is_prime(int64_t n);

// find a non-trivial factor of n
int64_t find_factor(int64_t n);

// Factor n (>= 1) into primes. Result is in ascending order.
std::vector<int64_t> factor(int64_t n);

// Jacobi symbol (a/n) defined for any odd integer n
int jacobi(int64_t a, int64_t n);

// Euler's totient function
int64_t phi(int64_t a);

// order of a mod m. returns 0 if gcd(a,b) != 1
int64_t ordermod(int64_t a, int64_t m);

// factorial(n) = n! = n * (n-1) * ... * 2 * 1
int64_t factorial(int64_t n);

// binomial(n,k) = (n choose k) = n! / (k! * (n-k)!)
int64_t binomial(int64_t n, int64_t k);

// Residual classes [a] = a + mℤ with a,m ∈ ℤ, m > 0
//   - is a field iff m is prime
//   - might be refactored into general quotient-ring-class
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

// properties of elements
inline bool invertible(IntMod a) { return gcd(a.a(), a.m()) == 1; }
inline bool regular(IntMod a) { return gcd(a.a(), a.m()) == 1; }

// unary operations
inline IntMod operator-(IntMod a)
{
	return IntMod(negmod(a.a(), a.m()), a.m());
}
inline IntMod inverse(IntMod a) { return IntMod(invmod(a.a(), a.m()), a.m()); }

// binary operations
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

// binary assigns
inline void operator+=(IntMod &a, IntMod b) { a = a + b; }
inline void operator-=(IntMod &a, IntMod b) { a = a - b; }
inline void operator*=(IntMod &a, IntMod b) { a = a * b; }
inline void operator/=(IntMod &a, IntMod b) { a = a / b; }

// comparisons
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
