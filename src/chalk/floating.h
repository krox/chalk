#ifndef CHALK_FLOATING_H
#define CHALK_FLOATING_H

/**
 * C++ wrapper for arbitrary precision floating point numbers using GNU MPFR.
 */

#include <stdio.h> // needs to be before mpfr.h

#include "fmt/format.h"
#include <cassert>
#include <mpfr.h>
#include <stdexcept>
#include <string>

namespace chalk {

/**
 * As far as possible a drop-in replacement for 'double'. New values are
 * computed with precision determined by the global Floating::precision.
 * All computations use "round to nearest" (MPFR_RNDN).
 */
class Floating
{
  public:
	// internal handle is public, but use with care
	// (in particular, mpfr_init/clear is called by constructor/destructor)
	mpfr_t f_;

	// precision (in mantissa bit) to compute everything in
	// (should be set at the beginning of any calculation)
	inline static mpfr_prec_t precision = 237; // ieee-754 "octuple precision"

	/** constructors / destructor / moves */
	Floating() { mpfr_init2(f_, precision); }
	explicit Floating(double value)
	{
		mpfr_init2(f_, precision);
		mpfr_set_d(f_, value, MPFR_RNDN);
	}
	explicit Floating(std::string const &value)
	{
		mpfr_init2(f_, precision);
		mpfr_set_str(f_, value.c_str(), 0, MPFR_RNDN);
	}
	~Floating() { mpfr_clear(f_); }
	Floating(Floating const &other)
	{
		mpfr_init2(f_, precision);
		mpfr_set(f_, other.f_, MPFR_RNDN);
	}
	Floating(Floating &&other)
	{
		mpfr_init2(f_, precision);
		mpfr_swap(f_, other.f_);
	}
	void operator=(Floating const &other) { mpfr_set(f_, other.f_, MPFR_RNDN); }
	void operator=(Floating &&other) { mpfr_swap(f_, other.f_); }
	void operator=(double value) { mpfr_set_d(f_, value, MPFR_RNDN); }

	/** convert to double/string */
	double to_double() const { return (double)mpfr_get_d(f_, MPFR_RNDN); }
	std::string to_string() const
	{
		char *tmp = nullptr;
		if (mpfr_asprintf(&tmp, "%Re", f_) < 0)
			throw std::runtime_error("error when calling mpfr_asprintf()");
		auto s = std::string(tmp);
		mpfr_free_str(tmp);
		return s;
	}
	explicit operator double() const { return to_double(); }

	/** mathematical constants */
	static Floating pi()
	{
		Floating r;
		mpfr_const_pi(r.f_, MPFR_RNDN);
		return r;
	}
	static Floating log2()
	{
		Floating r;
		mpfr_const_log2(r.f_, MPFR_RNDN);
		return r;
	}
	static Floating euler()
	{
		Floating r;
		mpfr_const_euler(r.f_, MPFR_RNDN);
		return r;
	}
	static Floating catalan()
	{
		Floating r;
		mpfr_const_catalan(r.f_, MPFR_RNDN);
		return r;
	}
};

/** unary arithmetic */
inline Floating operator+(Floating const &a) { return a; }
inline Floating operator-(Floating const &a)
{
	Floating r;
	mpfr_neg(r.f_, a.f_, MPFR_RNDN);
	return r;
}

/** binary arithmetic (Floating <-> Floating) */
inline Floating operator+(Floating const &a, Floating const &b)
{
	Floating r;
	mpfr_add(r.f_, a.f_, b.f_, MPFR_RNDN);
	return r;
}
inline Floating operator-(Floating const &a, Floating const &b)
{
	Floating r;
	mpfr_sub(r.f_, a.f_, b.f_, MPFR_RNDN);
	return r;
}
inline Floating operator*(Floating const &a, Floating const &b)
{
	Floating r;
	mpfr_mul(r.f_, a.f_, b.f_, MPFR_RNDN);
	return r;
}
inline Floating operator/(Floating const &a, Floating const &b)
{
	Floating r;
	mpfr_div(r.f_, a.f_, b.f_, MPFR_RNDN);
	return r;
}

/** assign arithmetic (Floating <-> Floating) */
inline void operator+=(Floating &a, Floating const &b)
{
	mpfr_add(a.f_, a.f_, b.f_, MPFR_RNDN);
}
inline void operator-=(Floating &a, Floating const &b)
{
	mpfr_sub(a.f_, a.f_, b.f_, MPFR_RNDN);
}
inline void operator*=(Floating &a, Floating const &b)
{
	mpfr_mul(a.f_, a.f_, b.f_, MPFR_RNDN);
}
inline void operator/=(Floating &a, Floating const &b)
{
	mpfr_div(a.f_, a.f_, b.f_, MPFR_RNDN);
}

/** binary arithmetic (Floating <-> double) */
inline Floating operator+(Floating const &a, double b)
{
	Floating r;
	mpfr_add_d(r.f_, a.f_, b, MPFR_RNDN);
	return r;
}
inline Floating operator-(Floating const &a, double b)
{
	Floating r;
	mpfr_sub_d(r.f_, a.f_, b, MPFR_RNDN);
	return r;
}
inline Floating operator*(Floating const &a, double b)
{
	Floating r;
	mpfr_mul_d(r.f_, a.f_, b, MPFR_RNDN);
	return r;
}
inline Floating operator/(Floating const &a, double b)
{
	Floating r;
	mpfr_div_d(r.f_, a.f_, b, MPFR_RNDN);
	return r;
}
inline Floating operator+(double a, Floating const &b)
{
	Floating r;
	mpfr_add_d(r.f_, b.f_, a, MPFR_RNDN);
	return r;
}
inline Floating operator-(double a, Floating const &b)
{
	Floating r;
	mpfr_d_sub(r.f_, a, b.f_, MPFR_RNDN);
	return r;
}
inline Floating operator*(double a, Floating const &b)
{
	Floating r;
	mpfr_mul_d(r.f_, b.f_, a, MPFR_RNDN);
	return r;
}
inline Floating operator/(double a, Floating const &b)
{
	Floating r;
	mpfr_d_div(r.f_, a, b.f_, MPFR_RNDN);
	return r;
}

/** assign arithmetic (Floating <-> double) */
inline void operator+=(Floating &a, double b)
{
	mpfr_add_d(a.f_, a.f_, b, MPFR_RNDN);
}
inline void operator-=(Floating &a, double b)
{
	mpfr_sub_d(a.f_, a.f_, b, MPFR_RNDN);
}
inline void operator*=(Floating &a, double b)
{
	mpfr_mul_d(a.f_, a.f_, b, MPFR_RNDN);
}
inline void operator/=(Floating &a, double b)
{
	mpfr_div_d(a.f_, a.f_, b, MPFR_RNDN);
}

/** comparison Floating <-> Floating (false if unordered, just as IEEE)*/
inline bool operator>(Floating const &a, Floating const &b)
{
	return mpfr_greater_p(a.f_, b.f_);
}
inline bool operator>=(Floating const &a, Floating const &b)
{
	return mpfr_greaterequal_p(a.f_, b.f_);
}
inline bool operator<(Floating const &a, Floating const &b)
{
	return mpfr_less_p(a.f_, b.f_);
}
inline bool operator<=(Floating const &a, Floating const &b)
{
	return mpfr_lessequal_p(a.f_, b.f_);
}
inline bool operator==(Floating const &a, Floating const &b)
{
	return mpfr_equal_p(a.f_, b.f_);
}

/** comparison Floating <-> double  */
inline bool operator<(Floating const &a, double b)
{
	return mpfr_cmp_d(a.f_, b) < 0; // 0 if either is NaN
}

/** mathematical functions of one parameter */
#define CHALK_DEFINE_FLOATING_FUNCTION(name)                                   \
	inline Floating name(Floating const &a)                                    \
	{                                                                          \
		Floating r;                                                            \
		mpfr_##name(r.f_, a.f_, MPFR_RNDN);                                    \
		return r;                                                              \
	}

// logarithms / exponentials
CHALK_DEFINE_FLOATING_FUNCTION(log)
CHALK_DEFINE_FLOATING_FUNCTION(log2)
CHALK_DEFINE_FLOATING_FUNCTION(log10)
CHALK_DEFINE_FLOATING_FUNCTION(log1p)
CHALK_DEFINE_FLOATING_FUNCTION(exp)
CHALK_DEFINE_FLOATING_FUNCTION(exp2)
CHALK_DEFINE_FLOATING_FUNCTION(exp10)
CHALK_DEFINE_FLOATING_FUNCTION(expm1)

// trigonometric / hyperbolic functions
CHALK_DEFINE_FLOATING_FUNCTION(sin)
CHALK_DEFINE_FLOATING_FUNCTION(cos)
CHALK_DEFINE_FLOATING_FUNCTION(tan)
CHALK_DEFINE_FLOATING_FUNCTION(sec)
CHALK_DEFINE_FLOATING_FUNCTION(csc)
CHALK_DEFINE_FLOATING_FUNCTION(cot)
CHALK_DEFINE_FLOATING_FUNCTION(asin)
CHALK_DEFINE_FLOATING_FUNCTION(acos)
CHALK_DEFINE_FLOATING_FUNCTION(atan)
CHALK_DEFINE_FLOATING_FUNCTION(sinh)
CHALK_DEFINE_FLOATING_FUNCTION(cosh)
CHALK_DEFINE_FLOATING_FUNCTION(tanh)
CHALK_DEFINE_FLOATING_FUNCTION(sech)
CHALK_DEFINE_FLOATING_FUNCTION(csch)
CHALK_DEFINE_FLOATING_FUNCTION(coth)
CHALK_DEFINE_FLOATING_FUNCTION(asinh)
CHALK_DEFINE_FLOATING_FUNCTION(acosh)
CHALK_DEFINE_FLOATING_FUNCTION(atanh)

// special functions
CHALK_DEFINE_FLOATING_FUNCTION(gamma)
CHALK_DEFINE_FLOATING_FUNCTION(eint)
CHALK_DEFINE_FLOATING_FUNCTION(lngamma)
CHALK_DEFINE_FLOATING_FUNCTION(digamma)
CHALK_DEFINE_FLOATING_FUNCTION(zeta)
CHALK_DEFINE_FLOATING_FUNCTION(erf)
CHALK_DEFINE_FLOATING_FUNCTION(erfc)

} // namespace chalk

template <> struct fmt::formatter<chalk::Floating>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Floating &x, FormatContext &ctx)
	{
		return format_to(ctx.out(), "{}", x.to_string());
	}
};

#endif
