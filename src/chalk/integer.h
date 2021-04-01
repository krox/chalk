#ifndef CHALK_INTEGER_H
#define CHALK_INTEGER_H

/**
 * C++ wrapper for arbitrary precision integers.
 * In contrast to the c++ wrapper provided by GMP itself, this one does not use
 * expression-templates.
 */

#include "chalk/rings.h"
#include "fmt/format.h"
#include <cassert>
#include <gmp.h>
#include <string>

namespace chalk {

/**
 * As far as possible a drop-in replacement for 'int' (or int128_t).
 * Differences include:
 *   - no overflows, no undefined behaviour (thats the whole point of this)
 *   - always initialized to 0 (even local variables)
 *   - expensive to copy (even after default-initialization)
 *   - 'a % b' ignores the sign of b, result is always non-negative
 *   - no divison operator (dont need it, rounding might lead to bugs)
 */
class Integer
{
  public:
	// TODO: i think a small-size optimization could be very worthwhile here
	//       (i.e. use int64_t or int128_t as long as possible). After all,
	//       that was one of the reasons I wrote this wrapper instead of just
	//       using the 'mpz_class' type of GMP.

	mpz_t z_; // semi-private: use with care!

	/** constructors / destructor / moves */
	Integer() { mpz_init(z_); }
	explicit Integer(int value) { mpz_init_set_si(z_, value); }
	explicit Integer(std::string const &value)
	{
		mpz_init_set_str(z_, value.c_str(), 0);
	}
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

	/** convert to double/string */
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

	/** arithmetic (unary)*/
	Integer operator+() const { return *this; }
	Integer operator-() const
	{
		Integer r;
		mpz_neg(r.z_, z_);
		return r;
	}

	/** arithmetic (binary Integer <-> Integer) */
	Integer operator+(Integer const &b) const
	{
		Integer r;
		mpz_add(r.z_, z_, b.z_);
		return r;
	}
	Integer operator-(Integer const &b) const
	{
		Integer r;
		mpz_sub(r.z_, z_, b.z_);
		return r;
	}
	Integer operator*(Integer const &b) const
	{
		Integer r;
		mpz_mul(r.z_, z_, b.z_);
		return r;
	}
	Integer operator/(Integer const &b) const
	{
		Integer r;
		mpz_cdiv_q(r.z_, z_, b.z_);
		return r;
	}
	Integer operator%(Integer const &b) const
	{
		Integer r;
		mpz_mod(r.z_, z_, b.z_);
		return r;
	}

	/** arithmetic (binary assigment Integer <-> Integer) */
	void operator+=(Integer const &b) { mpz_add(z_, z_, b.z_); }
	void operator-=(Integer const &b) { mpz_sub(z_, z_, b.z_); }
	void operator*=(Integer const &b) { mpz_mul(z_, z_, b.z_); }
	void operator/=(Integer const &b) { mpz_cdiv_q(z_, z_, b.z_); }
	void operator%=(Integer const &b) { mpz_mod(z_, z_, b.z_); }

	/** arithmetic (binary Integer <-> int) */
	Integer operator+(int b) const
	{
		Integer r;
		if (b >= 0)
			mpz_add_ui(r.z_, z_, b);
		else
		{
			// there is a subtle bug if b==INT_MIN on platforms with int==long
			assert(b > LONG_MIN);
			mpz_sub_ui(r.z_, z_, -(long)b);
		}
		return r;
	}

	Integer operator-(int b) const
	{
		Integer r;
		if (b >= 0)
			mpz_sub_ui(r.z_, z_, b);
		else
		{
			assert(b > LONG_MIN);
			mpz_add_ui(r.z_, z_, -(long)b);
		}
		return r;
	}
	Integer operator*(int b) const
	{
		Integer r;
		mpz_mul_si(r.z_, z_, b);
		return r;
	}
	int operator%(int b) const
	{
		assert(b > 0);
		Integer r;
		return mpz_mod_ui(r.z_, z_, b);
	}

	/** arithmetic (binary assign Integer <-> int) */
	void operator+=(int b)
	{
		if (b >= 0)
			mpz_add_ui(z_, z_, b);
		else
		{
			assert(b > LONG_MIN);
			mpz_sub_ui(z_, z_, -(long)b);
		}
	}
	void operator-=(int b)
	{
		if (b >= 0)
			mpz_sub_ui(z_, z_, b);
		else
		{
			assert(b > LONG_MIN);
			mpz_add_ui(z_, z_, -(long)b);
		}
	}
	void operator*=(int b) { mpz_mul_si(z_, z_, b); }
	void operator/=(int b)
	{
		assert(b > 0);
		mpz_cdiv_q_ui(z_, z_, b);
	}
	void operator%=(int b)
	{
		assert(b > 0);
		mpz_mod_ui(z_, z_, b);
	}

	friend Integer gcd(Integer const &, Integer const &);

	/** comaprisons */
	bool operator==(Integer const &b) const { return mpz_cmp(z_, b.z_) == 0; }
	bool operator<(Integer const &b) const { return mpz_cmp(z_, b.z_) < 0; }
	bool operator==(int b) const { return mpz_cmp_si(z_, b) == 0; }
	bool operator<(int b) const { return mpz_cmp_si(z_, b) < 0; }
	int sign() const { return mpz_sgn(z_); } // 1 / -1 / 0
};

inline Integer operator+(int a, Integer const &b) { return b + a; }
inline Integer operator*(int a, Integer const &b) { return b * a; }

inline Integer gcd(Integer const &a, Integer const &b)
{
	Integer r;
	mpz_gcd(r.z_, a.z_, b.z_);
	return r;
}

inline bool isPrime(Integer const &a)
{
	// 2 = definitely prime
	// 1 = probably prime ( error rate ~ (1/4)^reps )
	// 0 = definitely composite
	int reps = 32;
	return mpz_probab_prime_p(a.z_, reps) != 0;
}

inline Integer nextPrime(Integer const &a)
{
	Integer r;
	mpz_nextprime(r.z_, a.z_);
	return r;
}

template <> struct RingTraits<Integer> : RingTraitsSimple<Integer>
{};

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

#endif
