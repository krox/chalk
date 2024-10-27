#pragma once

/**
 * C++ wrapper for arbitrary precision integers.
 * In contrast to the c++ wrapper provided by GMP itself, this one does not use
 * expression-templates.
 */

#include "chalk/parser.h"
#include "chalk/rings.h"
#include "fmt/format.h"
#include <cassert>
#include <gmp.h>
#include <optional>
#include <span>
#include <string>

namespace chalk {

// As far as possible a drop-in replacement for 'int' with arbitrary precision
//   - implicitly convertible from int
//   - For compatibility, the division operator '/' has the same rounding
//     semantics as builin 'int'. If you know that the result is exact, you
//     should use the faster 'div_exact' instead.
//   - most functions are marked 'noexcept' because GMP can not handle
//     out-of-memory gracefully anyway.
//   - NOTE on GMP internals: the data pointer of mpz_t is always non-null and
//     some routines read one limb uncontitionally, even if its content is
//     undefined for 'size'=0. mpz_init thus sets the pointer to some global
//     dummy value to avoid allocation.
class Integer
{
  public:
	// TODO: i think a small-size optimization could be very worthwhile here
	//       (i.e. use int64_t or int128_t as long as possible). After all,
	//       that was one of the reasons I wrote this wrapper instead of just
	//       using the 'mpz_class' type of GMP.
	//       On Second thought, using the underlying 'mpn' functions directly,
	//       i.e., reimplementing the 'mpz' layer completely could be a cool
	//       project, and might lead to the most effective
	//       small-size-optimization.

	mpz_t z_; // semi-private: use with care!

	// constructors / destructor / moves

	// default constructor initializes to zero
	// (note: does not allocate, though it does some initialization)
	Integer() noexcept { mpz_init(z_); }

	// construct from integer
	Integer(int value) noexcept { mpz_init_set_si(z_, value); }

	// takes value as simple string.
	//     * defaults to decimal, but also supports hex (0x...)
	//     * throws on ill-formatted strings
	explicit Integer(char const *s)
	{
		assert(s);
		int base = 10;
		if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X'))
		{
			base = 16;
			s += 2;
		}
		// NOTE: GMP ignores whitespace, and can also do its own base-detection,
		// but we dont want to expose that really.
		int r = mpz_init_set_str(z_, s, base);
		if (r != 0)
		{
			// NOTE: mpz_init_set_str returns -1 on ill-formatted strings, but
			// still initializes the mpz_t, so we have to clear it manually.
			mpz_clear(z_);
			throw std::invalid_argument("invalid integer string");
		}
	}
	explicit Integer(std::string const &value) : Integer(value.c_str()) {}
	explicit Integer(std::string_view value) : Integer(std::string(value)) {}

	// pseudo-constructor that understands basic arithmetic expressions such as
	// "123+456" or "2^1000"
	static Integer parse(std::string_view expr)
	{
		return chalk::parse<Integer>(
		    expr, [](std::string_view s) { return Integer(s); });
	}

	~Integer() { mpz_clear(z_); }
	Integer(Integer const &other) noexcept { mpz_init_set(z_, other.z_); }
	Integer(Integer &&other) noexcept
	{
		mpz_init(z_);
		mpz_swap(z_, other.z_);
	}
	void operator=(Integer const &other) noexcept { mpz_set(z_, other.z_); }
	void operator=(Integer &&other) noexcept { mpz_swap(z_, other.z_); }
	void operator=(int value) noexcept { mpz_set_si(z_, value); }

	void swap(Integer &other) noexcept { mpz_swap(z_, other.z_); }
	friend void swap(Integer &a, Integer &b) noexcept { a.swap(b); }

	// convert to int/double/string

	bool fits_int() const noexcept { return mpz_fits_sint_p(z_); }
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
	explicit operator bool() const { return mpz_sgn(z_); }

	// low-level access to the data
	//     * sign is stored separately, raw data is positive
	//     * returns empty slice for zero, otherwise highest limb is non-zero
	std::span<const mp_limb_t> limbs() const noexcept
	{
		return std::span(mpz_limbs_read(z_), mpz_size(z_));
	}

	// uniform random integr in [0, m]
	// (GMP has its rng function, but I like to use my own generators)
	template <class Rng> static Integer uniform(Integer const &m, Rng &rng)
	{
		// only need the simple case when the rng already produces limbs
		static_assert(sizeof(typename Rng::result_type) == 8 &&
		              Rng::min() == 0 && Rng::max() == UINT64_MAX);
		static_assert(sizeof(mp_limb_t) == 8 && GMP_LIMB_BITS == 64 &&
		              GMP_NUMB_BITS == 64);

		if (mpz_sgn(m.z_) <= 0)
			return Integer(0);

		auto s = mpz_size(m.z_);
		Integer r;
		auto ptr = mpz_limbs_write(r.z_, s);
		auto m_ptr = mpz_limbs_read(m.z_);

		// Acceptance probability > 50%. Worst case (m is a power 2^64 or
		// slightly above) could be optimized as most rejections happen in the
		// second limb before generating eveything.
		while (true)
		{
			ptr[s - 1] = rng.uniform(m_ptr[s - 1]);
			for (size_t i = 0; i < s - 1; ++i)
				ptr[i] = rng();
			if (mpn_cmp(ptr, m_ptr, s) <= 0)
				break;
		}

		mpz_limbs_finish(r.z_, s);
		return r;
	}
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
inline void mpz_tdiv_q_si(mpz_ptr r, mpz_srcptr a, long b)
{
	assert(b >= 0);
	mpz_tdiv_q_ui(r, a, (unsigned long)b);
}
inline void mpz_tdiv_r_si(mpz_ptr r, mpz_srcptr a, long b)
{
	assert(b >= 0);
	mpz_tdiv_r_ui(r, a, (unsigned long)b);
}

// basic arithmetic

#define CHALK_DEFINE_INTEGER_UNARY(fun, gmp_fun)                               \
	inline Integer fun(Integer const &a) noexcept                              \
	{                                                                          \
		Integer r;                                                             \
		gmp_fun(r.z_, a.z_);                                                   \
		return r;                                                              \
	}                                                                          \
	inline Integer fun(Integer &&a) noexcept                                   \
	{                                                                          \
		gmp_fun(a.z_, a.z_);                                                   \
		return std::move(a);                                                   \
	}

#define CHALK_DEFINE_INTEGER_BINARY(fun, funInplace, gmp_fun)                  \
	inline Integer fun(Integer const &a, Integer const &b) noexcept            \
	{                                                                          \
		Integer r;                                                             \
		gmp_fun(r.z_, a.z_, b.z_);                                             \
		return r;                                                              \
	}                                                                          \
	inline Integer fun(Integer &&a, Integer const &b) noexcept                 \
	{                                                                          \
		gmp_fun(a.z_, a.z_, b.z_);                                             \
		return std::move(a);                                                   \
	}                                                                          \
	inline Integer fun(Integer const &a, Integer &&b) noexcept                 \
	{                                                                          \
		gmp_fun(b.z_, a.z_, b.z_);                                             \
		return std::move(b);                                                   \
	}                                                                          \
	inline Integer fun(Integer &&a, Integer &&b) noexcept                      \
	{                                                                          \
		gmp_fun(a.z_, a.z_, b.z_);                                             \
		return std::move(a);                                                   \
	}                                                                          \
	inline Integer fun(Integer const &a, int b) noexcept                       \
	{                                                                          \
		Integer r;                                                             \
		gmp_fun##_si(r.z_, a.z_, b);                                           \
		return r;                                                              \
	}                                                                          \
	inline Integer fun(Integer &&a, int b) noexcept                            \
	{                                                                          \
		gmp_fun##_si(a.z_, a.z_, b);                                           \
		return std::move(a);                                                   \
	}                                                                          \
	inline Integer &funInplace(Integer &a, Integer const &b) noexcept          \
	{                                                                          \
		gmp_fun(a.z_, a.z_, b.z_);                                             \
		return a;                                                              \
	}                                                                          \
	inline Integer &funInplace(Integer &a, int b) noexcept                     \
	{                                                                          \
		gmp_fun##_si(a.z_, a.z_, b);                                           \
		return a;                                                              \
	}

CHALK_DEFINE_INTEGER_UNARY(operator-, mpz_neg)
CHALK_DEFINE_INTEGER_UNARY(next_prime, mpz_nextprime)

CHALK_DEFINE_INTEGER_BINARY(operator+, operator+=, mpz_add)
CHALK_DEFINE_INTEGER_BINARY(operator-, operator-=, mpz_sub)
CHALK_DEFINE_INTEGER_BINARY(operator*, operator*=, mpz_mul)
CHALK_DEFINE_INTEGER_BINARY(operator/, operator/=, mpz_tdiv_q)
CHALK_DEFINE_INTEGER_BINARY(div_exact, div_exact_assign, mpz_divexact)
CHALK_DEFINE_INTEGER_BINARY(operator%, operator%=, mpz_tdiv_r)
CHALK_DEFINE_INTEGER_BINARY(mod, mod_assign, mpz_mod)

#undef CHALK_DEFINE_INTEGER_UNARY
#undef CHALK_DEFINE_INTEGER_BINARY

inline Integer pow(Integer const &a, unsigned long b) noexcept
{
	Integer r;
	mpz_pow_ui(r.z_, a.z_, b);
	return r;
}

inline Integer pow(Integer const &a, int b) noexcept
{
	assert(b >= 0);
	return pow(a, (unsigned long)b);
}

// comparisons

inline bool operator==(Integer const &a, Integer const &b) noexcept
{
	return mpz_cmp(a.z_, b.z_) == 0;
}

inline std::strong_ordering operator<=>(Integer const &a,
                                        Integer const &b) noexcept
{
	int r = mpz_cmp(a.z_, b.z_);
	if (r == 0)
		return std::strong_ordering::equal;
	if (r < 0)
		return std::strong_ordering::less;
	return std::strong_ordering::greater;
}

inline bool operator==(Integer const &a, int b) noexcept
{
	return mpz_cmp_si(a.z_, b) == 0;
}

inline std::strong_ordering operator<=>(Integer const &a, int b) noexcept
{
	int r = mpz_cmp_si(a.z_, b);
	if (r == 0)
		return std::strong_ordering::equal;
	if (r < 0)
		return std::strong_ordering::less;
	return std::strong_ordering::greater;
}

inline int sign(Integer const &a) noexcept
{
	return mpz_sgn(a.z_);
} // 1 / -1 / 0

inline std::optional<Integer> trySqrt(Integer const &a) noexcept
{
	if (!mpz_perfect_square_p(a.z_))
		return {};
	Integer r;
	mpz_sqrt(r.z_, a.z_);
	return r;
}

inline Integer operator+(int a, Integer const &b) noexcept { return b + a; }
inline Integer operator*(int a, Integer const &b) noexcept { return b * a; }

// greatest common divisor
//    - signs of input are ignored, result is always non-negative
//    - convention: gcd(0,x) = abs(x)
inline Integer gcd(Integer const &a, Integer const &b) noexcept
{
	Integer r;
	mpz_gcd(r.z_, a.z_, b.z_);
	return r;
}

inline Integer gcd(Integer const &a, Integer const &b,
                   Integer const &c) noexcept
{
	Integer r;
	mpz_gcd(r.z_, a.z_, b.z_);
	mpz_gcd(r.z_, r.z_, c.z_);
	return r;
}

inline bool divisible(Integer const &n, Integer const &d) noexcept
{
	return mpz_divisible_p(n.z_, d.z_) != 0;
}

// divide a,b by their greatest common divisor. Also makes a non-negative.
inline Integer remove_common_factor(Integer &a, Integer &b) noexcept
{
	auto d = gcd(a, b);
	if (sign(a) < 0)
		mpz_neg(d.z_, d.z_);
	mpz_divexact(a.z_, a.z_, d.z_);
	mpz_divexact(b.z_, b.z_, d.z_);
	return d;
}

inline Integer remove_common_factor(Integer &a, Integer &b, Integer &c) noexcept
{
	auto d = gcd(a, b, c);
	if (sign(a) < 0)
		mpz_neg(d.z_, d.z_);
	mpz_divexact(a.z_, a.z_, d.z_);
	mpz_divexact(b.z_, b.z_, d.z_);
	mpz_divexact(c.z_, c.z_, d.z_);
	return d;
}

inline bool is_prime(Integer const &a) noexcept
{
	// 2 = definitely prime
	// 1 = probably prime ( error rate ~ (1/4)^reps )
	// 0 = definitely composite
	int reps = 32;
	return mpz_probab_prime_p(a.z_, reps) != 0;
}

inline Integer removeSquareFactor(Integer &a) noexcept
{
	auto r = Integer(1);
	for (auto p = Integer(2); !(a < p * p); p = next_prime(p))
		if (divisible(a, p * p))
		{
			div_exact_assign(a, p * p);
			r *= p;
			p -= 1;
		}
	return r;
}

inline Integer powmod(Integer const &b, Integer const &e,
                      Integer const &m) noexcept
{
	Integer r;
	mpz_powm(r.z_, b.z_, e.z_, m.z_);
	return r;
}

inline Integer powmod(Integer const &b, unsigned long e,
                      Integer const &m) noexcept
{
	Integer r;
	mpz_powm_ui(r.z_, b.z_, e, m.z_);
	return r;
}

inline Integer invmod(Integer const &a, Integer const &m) noexcept
{
	// m=1 is supported as inverse(0)=0
	// m=0 is undefined
	// returns 0 if inverse does not exist
	Integer r;
	[[maybe_unused]] auto success = mpz_invert(r.z_, a.z_, m.z_);
	assert(success);
	return r;
}

template <> struct RingTraits<Integer> : RingTraitsSimple<Integer>
{
	static bool is_zero(Integer const &a) noexcept
	{
		return mpz_sgn(a.z_) == 0;
	}
	static bool is_negative(Integer const &a) noexcept
	{
		return mpz_sgn(a.z_) < 0;
	}
};

// makes Integer usable on the command line using CLI11 library (found by ADL)
inline bool lexical_cast(std::string const &input, Integer &v)
{
	v = Integer(input);
	return true;
}

} // namespace chalk

template <> struct fmt::formatter<chalk::Integer>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Integer &x, FormatContext &ctx) const
	{
		return fmt::format_to(ctx.out(), "{}", x.to_string());
	}
};
