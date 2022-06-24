#pragma once

/**
 * C++ wrapper for arbitrary precision integers.
 * In contrast to the c++ wrapper provided by GMP itself, this one does not use
 * expression-templates.
 */

#include "chalk/rings.h"
#include "fmt/format.h"
#include <cassert>
#include <gmp.h>
#include <optional>
#include <string>

namespace chalk {

/**
 * As far as possible a drop-in replacement for 'int' (or int128_t).
 * Differences include:
 *   - no overflows, no undefined behaviour (thats the whole point of this)
 *   - always initialized to 0 (even local variables)
 *   - expensive to copy (even after default-initialization)
 *   - 'a % b' ignores the sign of b, result is always non-negative
 */
class Integer
{
  public:
	// TODO: i think a small-size optimization could be very worthwhile here
	//       (i.e. use int64_t or int128_t as long as possible). After all,
	//       that was one of the reasons I wrote this wrapper instead of just
	//       using the 'mpz_class' type of GMP.

	mpz_t z_; // semi-private: use with care!

	// constructors / destructor / moves

	Integer() { mpz_init(z_); }
	explicit Integer(int value) { mpz_init_set_si(z_, value); }
	explicit Integer(std::string const &value)
	{
		mpz_init_set_str(z_, value.c_str(), 0);
	}
	explicit Integer(std::string_view value) : Integer(std::string(value)) {}
	~Integer() { mpz_clear(z_); }
	Integer(Integer const &other) { mpz_init_set(z_, other.z_); }
	Integer(Integer &&other)
	{
		mpz_init(z_);
		mpz_swap(z_, other.z_);
	}
	void operator=(Integer const &other) { mpz_set(z_, other.z_); }
	void operator=(Integer &&other) { mpz_swap(z_, other.z_); }
	void operator=(int value) { mpz_set_si(z_, value); }

	// convert to int/double/string

	bool fits_int() const { return mpz_fits_sint_p(z_); }
	int to_int() const
	{
		assert(mpz_fits_sint_p(z_));
		return (int)mpz_get_si(z_);
	}
	double to_double() const { return (double)mpz_get_d(z_); }
	std::string to_string() const
	{
		char *tmp = mpz_get_str(nullptr, 10, z_);
		auto s = std::string(tmp);
		void (*freefunc)(void *, size_t);
		mp_get_memory_functions(nullptr, nullptr, &freefunc);
		freefunc(tmp, s.size() + 1);
		return s;
	}
	explicit operator int() const { return to_int(); }
	explicit operator double() const { return to_double(); }
};

// functions that are missing from GMP to make the interface more uniform

inline void mpz_add_si(mpz_ptr r, mpz_srcptr a, long b)
{
	if (b >= 0)
		mpz_add_ui(r, a, (unsigned long)b);
	else
		mpz_sub_ui(r, a, -(unsigned long)b);
}
inline void mpz_sub_si(mpz_ptr r, mpz_srcptr a, long b)
{
	if (b >= 0)
		mpz_sub_ui(r, a, (unsigned long)b);
	else
		mpz_add_ui(r, a, -(unsigned long)b);
}
inline void mpz_mod_si(mpz_ptr r, mpz_srcptr a, long b)
{
	mpz_mod_ui(r, a, b >= 0 ? (unsigned long)b : -(unsigned long)b);
}
inline void mpz_divexact_si(mpz_ptr r, mpz_srcptr a, long b)
{
	assert(b >= 0);
	mpz_mod_ui(r, a, (unsigned long)b);
}
inline void mpz_cdiv_q_si(mpz_ptr r, mpz_srcptr a, long b)
{
	assert(b >= 0);
	mpz_cdiv_q_ui(r, a, (unsigned long)b);
}

// basic arithmetic

#define CHALK_DEFINE_INTEGER_UNARY(fun, gmp_fun)                               \
	inline Integer fun(Integer const &a)                                       \
	{                                                                          \
		Integer r;                                                             \
		gmp_fun(r.z_, a.z_);                                                   \
		return r;                                                              \
	}

#define CHALK_DEFINE_INTEGER_BINARY(fun, funInplace, gmp_fun)                  \
	inline Integer fun(Integer const &a, Integer const &b)                     \
	{                                                                          \
		Integer r;                                                             \
		gmp_fun(r.z_, a.z_, b.z_);                                             \
		return r;                                                              \
	}                                                                          \
	inline Integer fun(Integer const &a, int b)                                \
	{                                                                          \
		Integer r;                                                             \
		gmp_fun##_si(r.z_, a.z_, b);                                           \
		return r;                                                              \
	}                                                                          \
	inline Integer &funInplace(Integer &a, Integer const &b)                   \
	{                                                                          \
		gmp_fun(a.z_, a.z_, b.z_);                                             \
		return a;                                                              \
	}                                                                          \
	inline Integer &funInplace(Integer &a, int b)                              \
	{                                                                          \
		gmp_fun##_si(a.z_, a.z_, b);                                           \
		return a;                                                              \
	}

CHALK_DEFINE_INTEGER_UNARY(operator-, mpz_neg)
CHALK_DEFINE_INTEGER_UNARY(nextPrime, mpz_nextprime)

CHALK_DEFINE_INTEGER_BINARY(operator+, operator+=, mpz_add)
CHALK_DEFINE_INTEGER_BINARY(operator-, operator-=, mpz_sub)
CHALK_DEFINE_INTEGER_BINARY(operator*, operator*=, mpz_mul)
CHALK_DEFINE_INTEGER_BINARY(operator/, operator/=, mpz_cdiv_q)
CHALK_DEFINE_INTEGER_BINARY(operator%, operator%=, mpz_mod)
CHALK_DEFINE_INTEGER_BINARY(divExact, divExactInplace, mpz_divexact)

#undef CHALK_DEFINE_INTEGER_UNARY
#undef CHALK_DEFINE_INTEGER_BINARY

// comparisons

inline bool operator==(Integer const &a, Integer const &b)
{
	return mpz_cmp(a.z_, b.z_) == 0;
}
inline bool operator<(Integer const &a, Integer const &b)
{
	return mpz_cmp(a.z_, b.z_) < 0;
}
inline bool operator==(Integer const &a, int b)
{
	return mpz_cmp_si(a.z_, b) == 0;
}
inline bool operator<(Integer const &a, int b)
{
	return mpz_cmp_si(a.z_, b) < 0;
}
inline bool operator<(int a, Integer const &b)
{
	return mpz_cmp_si(b.z_, a) > 0;
}

inline int sign(Integer const &a) { return mpz_sgn(a.z_); } // 1 / -1 / 0

inline std::optional<Integer> trySqrt(Integer const &a)
{
	if (!mpz_perfect_square_p(a.z_))
		return {};
	Integer r;
	mpz_sqrt(r.z_, a.z_);
	return r;
}

inline Integer operator+(int a, Integer const &b) { return b + a; }
inline Integer operator*(int a, Integer const &b) { return b * a; }

// greatest common divisor
//    - signs of input are ignored, result is always non-negative
//    - convention: gcd(0,x) = abs(x)
inline Integer gcd(Integer const &a, Integer const &b)
{
	Integer r;
	mpz_gcd(r.z_, a.z_, b.z_);
	return r;
}

inline Integer gcd(Integer const &a, Integer const &b, Integer const &c)
{
	Integer r;
	mpz_gcd(r.z_, a.z_, b.z_);
	mpz_gcd(r.z_, r.z_, c.z_);
	return r;
}

inline bool divisible(Integer const &n, Integer const &d)
{
	return mpz_divisible_p(n.z_, d.z_) != 0;
}

// divide a,b by their greatest common divisor. Also makes a non-negative.
inline Integer removeCommonFactor(Integer &a, Integer &b)
{
	auto d = gcd(a, b);
	if (sign(a) < 0)
		mpz_neg(d.z_, d.z_);
	mpz_divexact(a.z_, a.z_, d.z_);
	mpz_divexact(b.z_, b.z_, d.z_);
	return d;
}

inline Integer removeCommonFactor(Integer &a, Integer &b, Integer &c)
{
	auto d = gcd(a, b, c);
	if (sign(a) < 0)
		mpz_neg(d.z_, d.z_);
	mpz_divexact(a.z_, a.z_, d.z_);
	mpz_divexact(b.z_, b.z_, d.z_);
	mpz_divexact(c.z_, c.z_, d.z_);
	return d;
}

inline bool is_prime(Integer const &a)
{
	// 2 = definitely prime
	// 1 = probably prime ( error rate ~ (1/4)^reps )
	// 0 = definitely composite
	int reps = 32;
	return mpz_probab_prime_p(a.z_, reps) != 0;
}

inline Integer removeSquareFactor(Integer &a)
{
	auto r = Integer(1);
	for (auto p = Integer(2); !(a < p * p); p = nextPrime(p))
		if (divisible(a, p * p))
		{
			a = divExact(a, p * p);
			r *= p;
			p -= 1;
		}
	return r;
}

template <> struct RingTraits<Integer> : RingTraitsSimple<Integer>
{
	static bool isZero(Integer const &a) { return mpz_sgn(a.z_) == 0; }
	static bool isNegative(Integer const &a) { return mpz_sgn(a.z_) < 0; }
};

} // namespace chalk

template <> struct fmt::formatter<chalk::Integer>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Integer &x, FormatContext &ctx)
	{
		return format_to(ctx.out(), "{}", x.to_string());
	}
};
