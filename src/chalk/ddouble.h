#pragma once

#include "chalk/floating.h"
#include "util/span.h"
#include <cassert>
#include <cmath>
#include <random>

namespace chalk {

/**
 * "double double" type which is a (non-IEEE) floating point type implemented
 * as the sum of two doubles and thus can represet about 107 bits of effective
 * mantissa. In practice, this is this fastest possible datatype with more
 * precision than 'double'. (Much faster than libraries like MPFR, which start
 * to shine at much higher precision).
 *
 * Notes on correctness:
 *    - This module will break if compiled with '-ffast-math'. Just dont. Some
 *      weaker flags are okay though, e.g.:
 *      -fno-signed-zeros -fno-math-errno -ffinite-math-only
 *    - It will also probably break if FP rounding mode is changed away
 *      from 'nearest'. But who would ever consider that?
 *    - Correct handling of inf/nan/signed-zero/subnormals is not really tested.
 *      (And probably wont be fixed if broken. Dont want branches.)
 *
 * Notes on performance:
 *    - Make sure you compiler/standard-library emits actual fused multiply-add
 *      instructions for 'std::fma()' (supported since AVX2 on x86 CPUs).
 *      This will speed up all operations quite a lot.
 *    - Only use ddouble where necessary. Mixed operations between double
 *      and ddouble are all supported and should be prefered. In particular,
 *      use static member function like 'ddouble::sum(a,b)' to get ddouble
 *      results from double paramters.
 *    - For the last bit of performance, there, there are some specialized
 *      instructions that are (very) slightly faster than the naive expression
 *      (The compiler is probably less helpful here than it is with 'double')
 *            2 * x  ->  times2(x)     (same for times4, divide2, divide4)
 *          512 * x  ->  ldexp(x, 9)
 *            x * x  ->  sqr(x)        (is this true?)
 *            1 / x  ->  inverse(x)    (or this?)
 *
 * Notes on implementation:
 *    - When in doubt, choose speed over precision. Also dont really care about
 *      correct handling of NaN/infinity/signed-zero/subnormals
 *    - As few branches as possible, to make vectorization straightforward.
 */
class ddouble
{
	double high_;
	double low_;

	constexpr ddouble(double high, double low) : high_(high), low_(low)
	{
		assert(high_ + low_ == high);
	}

  public:
	explicit ddouble(double a) : high_(a), low_(0) {}

	static ddouble unchecked(double high, double low)
	{
		return ddouble(high, low);
	}

	/** should not really be needed in user code */
	double high() const { return high_; }
	double low() const { return low_; }

	/** constant naming as in std::numbers, though not exactly the same list */
	// generated with high precision MPFR, so should be precise in all digits
	// clang-format off
	static constexpr ddouble e()          { return ddouble(0x1.5bf0a8b145769p+1,  0x1.4d57ee2b1013ap-53); } // 2.718...
	static constexpr ddouble inv_e()      { return ddouble(0x1.78b56362cef38p-2, -0x1.ca8a4270fadf5p-57); }
	static constexpr ddouble egamma()     { return ddouble(0x1.2788cfc6fb619p-1, -0x1.6cb90701fbfabp-58); } // 0.577...
	static constexpr ddouble inv_egamma() { return ddouble(0x1.bb8226f502bf8p+0, -0x1.7abec73926687p-56); }
	static constexpr ddouble pi()         { return ddouble(0x1.921fb54442d18p+1,  0x1.1a62633145c07p-53); } // 3.141...
	static constexpr ddouble inv_pi()     { return ddouble(0x1.45f306dc9c883p-2, -0x1.6b01ec5417056p-56); }
	static constexpr ddouble sqrt2()      { return ddouble(0x1.6a09e667f3bcdp+0, -0x1.bdd3413b26456p-54); } // 1.414...
	static constexpr ddouble inv_sqrt2()  { return ddouble(0x1.6a09e667f3bcdp-1, -0x1.bdd3413b26456p-55); }
	static constexpr ddouble sqrt3()      { return ddouble(0x1.bb67ae8584caap+0,  0x1.cec95d0b5c1e3p-54); } // 1.732...
	static constexpr ddouble inv_sqrt3()  { return ddouble(0x1.279a74590331cp-1,  0x1.34863e0792bedp-55); }
	static constexpr ddouble ln2()        { return ddouble(0x1.62e42fefa39efp-1,  0x1.abc9e3b39803fp-56); } // 0.693...
	static constexpr ddouble log2e()      { return ddouble(0x1.71547652b82fep+0,  0x1.777d0ffda0d24p-56); }
	static constexpr ddouble ln10()       { return ddouble(0x1.26bb1bbb55516p+1, -0x1.f48ad494ea3e9p-53); } // 2.302...
	static constexpr ddouble log10e()     { return ddouble(0x1.bcb7b1526e50ep-2,  0x1.95355baaafad3p-57); }
	static constexpr ddouble phi()        { return ddouble(0x1.9e3779b97f4a8p+0, -0x1.f506319fcfd19p-55); } // 1.618...
	static constexpr ddouble inv_phi()    { return ddouble(0x1.3c6ef372fe950p-1, -0x1.f506319fcfd19p-55); }
	// clang-format on

	template <typename RNG> static ddouble random(RNG &rng)
	{
		auto dist = std::uniform_real_distribution<double>(0.0, 1.0);
		return sum_quick(dist(rng), ldexp(dist(rng), -53));
	}
	/** Compute a+b. Result is exact (except of overflows). No rounding. */
	static ddouble sum(double a, double b)
	{
		double high = a + b;
		double v = high - a;
		double low = (a - (high - v)) + (b - v);
		return ddouble(high, low);
	}

	/** Compute a+b, assuming abs(a) >= abs(b). */
	static ddouble sum_quick(double a, double b)
	{
		double high = a + b;
		double low = b - (high - a);
		return ddouble(high, low);
	}

	static ddouble mul(double a, double b)
	{
		double high = a * b;
		double low = fma(a, b, -high);
		return ddouble(high, low);
	}

	static ddouble div(double a, double b)
	{
		double high = a / b;
		double low = fma(-high, b, a) / b;
		return ddouble(high, low);
	}
};

/** unary operators */
inline ddouble operator-(ddouble a)
{
	return ddouble::unchecked(-a.high(), -a.low());
}

inline ddouble times2(ddouble a)
{
	return ddouble::unchecked(a.high() * 2, a.low() * 2);
}

inline ddouble times4(ddouble a)
{
	return ddouble::unchecked(a.high() * 4, a.low() * 4);
}

inline ddouble divide2(ddouble a)
{
	return ddouble::unchecked(a.high() / 2, a.low() / 2);
}

inline ddouble divide4(ddouble a)
{
	return ddouble::unchecked(a.high() / 4, a.low() / 4);
}

/** set r = a * 2^e */
inline ddouble ldexp(ddouble a, int e)
{
	return ddouble::unchecked(std::ldexp(a.high(), e), std::ldexp(a.low(), e));
}

/** decompose a = r * 2^e with 0.5 <= |r| < 1 */
inline ddouble frexp(ddouble a, int *e)
{
	double high = std::frexp(a.high(), e);
	double low = std::ldexp(a.low(), -*e);
	return ddouble::unchecked(high, low);
}

inline int ilogb(ddouble a) { return std::ilogb(a.high()); }

/** binary double <-> ddouble operators */
inline ddouble operator+(ddouble a, double b)
{
	ddouble tmp = ddouble::sum(a.high(), b);
	return ddouble::sum_quick(tmp.high(), tmp.low() + a.low());
}
inline ddouble operator+(double a, ddouble b) { return b + a; }
inline ddouble operator-(ddouble a, double b) { return a + (-b); }
inline ddouble operator-(double a, ddouble b) { return (-b) + a; }
inline ddouble operator*(ddouble a, double b)
{
	auto tmp = ddouble::mul(a.high(), b);
	return ddouble::sum_quick(tmp.high(), tmp.low() + a.low() * b);
}
inline ddouble operator*(double a, ddouble b) { return b * a; }
inline ddouble operator/(ddouble a, double b)
{
	double high = a.high() / b;
	double low = (fma(-high, b, a.high()) + a.low()) / b;
	return ddouble::sum_quick(high, low);
}

/** binary ddouble <-> ddouble operators */
inline ddouble operator+(ddouble a, ddouble b)
{
	ddouble tmp = ddouble::sum(a.high(), b.high());
	return ddouble::sum_quick(tmp.high(), tmp.low() + a.low() + b.low());
}
inline ddouble operator-(ddouble a, ddouble b) { return a + (-b); }
inline ddouble operator*(ddouble a, ddouble b)
{
	auto tmp = ddouble::mul(a.high(), b.high());
	return ddouble::sum_quick(tmp.high(), tmp.low() + a.high() * b.low() +
	                                          a.low() * b.high());
}
inline ddouble operator/(ddouble a, ddouble b)
{
	double high = a.high() / b.high();
	double low =
	    (fma(-high, b.high(), a.high()) + a.low() - high * b.low()) / b.high();
	return ddouble::sum_quick(high, low);
}

/** convenience operators */
inline void operator+=(ddouble &a, double b) { a = a + b; }
inline void operator-=(ddouble &a, double b) { a = a - b; }
inline void operator*=(ddouble &a, double b) { a = a * b; }
inline void operator/=(ddouble &a, double b) { a = a / b; }
inline void operator+=(ddouble &a, ddouble b) { a = a + b; }
inline void operator-=(ddouble &a, ddouble b) { a = a - b; }
inline void operator*=(ddouble &a, ddouble b) { a = a * b; }
inline void operator/=(ddouble &a, ddouble b) { a = a / b; }

namespace ddouble_detail {
ddouble taylor(util::span<const ddouble> c, ddouble x)
{
	ddouble r = c[0] + x * c[1];
	ddouble xi = x;

	for (size_t i = 2; i < c.size(); ++i)
	{
		xi *= x;
		r += c[i] * xi;
	}
	return r;
}

// taylor coefficients of exp
static ddouble coeffs_exp[] = {
    ddouble::unchecked(0x1.0000000000000p+0, 0x0.0000000000000p+0),
    ddouble::unchecked(0x1.0000000000000p+0, 0x0.0000000000000p+0),
    ddouble::unchecked(0x1.0000000000000p-1, 0x0.0000000000000p+0),
    ddouble::unchecked(0x1.5555555555555p-3, 0x1.5555555555555p-57),
    ddouble::unchecked(0x1.5555555555555p-5, 0x1.5555555555555p-59),
    ddouble::unchecked(0x1.1111111111111p-7, 0x1.1111111111111p-63),
    ddouble::unchecked(0x1.6c16c16c16c17p-10, -0x1.f49f49f49f49fp-65),
    ddouble::unchecked(0x1.a01a01a01a01ap-13, 0x1.a01a01a01a01ap-73),
    ddouble::unchecked(0x1.a01a01a01a01ap-16, 0x1.a01a01a01a01ap-76),
    ddouble::unchecked(0x1.71de3a556c734p-19, -0x1.c154f8ddc6c00p-73),
    ddouble::unchecked(0x1.27e4fb7789f5cp-22, 0x1.cbbc05b4fa99ap-76),
    ddouble::unchecked(0x1.ae64567f544e4p-26, -0x1.c062e06d1f209p-80),
    ddouble::unchecked(0x1.1eed8eff8d898p-29, -0x1.2aec959e14c06p-83),
    // ddouble::unchecked(0x1.6124613a86d09p-33, 0x1.f28e0cc748ebep-87),
    // ddouble::unchecked(0x1.93974a8c07c9dp-37, 0x1.05d6f8a2efd1fp-92),
    // ddouble::unchecked(0x1.ae7f3e733b81fp-41, 0x1.1d8656b0ee8cbp-097),
    // ddouble::unchecked(0x1.ae7f3e733b81fp-45, 0x1.1d8656b0ee8cbp-101),
};

// taylor coefficients of sin (odd only)
static ddouble coeffs_sin[] = {
    ddouble::unchecked(0x1.0000000000000p+0, 0x0.0000000000000p+0),
    ddouble::unchecked(-0x1.5555555555555p-3, -0x1.5555555555555p-57),
    ddouble::unchecked(0x1.1111111111111p-7, 0x1.1111111111111p-63),
    ddouble::unchecked(-0x1.a01a01a01a01ap-13, -0x1.a01a01a01a01ap-73),
    ddouble::unchecked(0x1.71de3a556c734p-19, -0x1.c154f8ddc6c00p-73),
    ddouble::unchecked(-0x1.ae64567f544e4p-26, 0x1.c062e06d1f209p-80),
    ddouble::unchecked(0x1.6124613a86d09p-33, 0x1.f28e0cc748ebep-87),
    ddouble::unchecked(-0x1.ae7f3e733b81fp-41, -0x1.1d8656b0ee8cbp-97),
    ddouble::unchecked(0x1.952c77030ad4ap-49, 0x1.ac981465ddc6cp-103),
    ddouble::unchecked(-0x1.2f49b46814157p-57, -0x1.2650f61dbdcb4p-112),
    ddouble::unchecked(0x1.71b8ef6dcf572p-66, -0x1.d043ae40c4647p-120),
    ddouble::unchecked(-0x1.761b41316381ap-75, 0x1.3423c7d91404fp-130),
    ddouble::unchecked(0x1.3f3ccdd165fa9p-84, -0x1.58ddadf344487p-139),
    ddouble::unchecked(-0x1.d1ab1c2dccea3p-94, -0x1.054d0c78aea14p-149),
    ddouble::unchecked(0x1.259f98b4358adp-103, 0x1.eaf8c39dd9bc5p-157),
    ddouble::unchecked(-0x1.434d2e783f5bcp-113, -0x1.0b87b91be9affp-167),
    ddouble::unchecked(0x1.3981254dd0d52p-123, -0x1.2b1f4c8015a2fp-177),
    ddouble::unchecked(-0x1.0dc59c716d91fp-133, -0x1.419e3fad3f031p-188),
    ddouble::unchecked(0x1.9ec8d1c94e85bp-144, -0x1.670e9d4784ec6p-201),
    ddouble::unchecked(-0x1.1e99449a4bacep-154, 0x1.fefbb89514b3cp-210),
    ddouble::unchecked(0x1.65e61c39d0241p-165, -0x1.c0ed181727269p-220),
    ddouble::unchecked(-0x1.95db45257e512p-176, -0x1.6e5d72b6f79b9p-231),
    ddouble::unchecked(0x1.a3cb872220648p-187, -0x1.c7f4e85b8e6cdp-241),
    ddouble::unchecked(-0x1.8da8e0a127ebap-198, 0x1.21d2eac9d275cp-252),
    ddouble::unchecked(0x1.5a42f0dfeb086p-209, -0x1.35ae015f78f6ep-264),
    ddouble::unchecked(-0x1.161872bf7b823p-220, -0x1.bb96c8e2e8897p-275),
    ddouble::unchecked(0x1.9d4f1058674dfp-232, 0x1.03c81b6914d59p-286),
    ddouble::unchecked(-0x1.1d008faac5c50p-243, -0x1.50348ded2636fp-298),
    ddouble::unchecked(0x1.6db793c887b97p-255, -0x1.966963ad60539p-314),
    ddouble::unchecked(-0x1.b5bfc17fa97d3p-267, 0x1.ff5794693c028p-321),
};

// taylor coefficients of cos (even only)
static ddouble coeffs_cos[] = {
    ddouble::unchecked(0x1.0000000000000p+0, 0x0.0000000000000p+0),
    ddouble::unchecked(-0x1.0000000000000p-1, 0x0.0000000000000p+0),
    ddouble::unchecked(0x1.5555555555555p-5, 0x1.5555555555555p-59),
    ddouble::unchecked(-0x1.6c16c16c16c17p-10, 0x1.f49f49f49f49fp-65),
    ddouble::unchecked(0x1.a01a01a01a01ap-16, 0x1.a01a01a01a01ap-76),
    ddouble::unchecked(-0x1.27e4fb7789f5cp-22, -0x1.cbbc05b4fa99ap-76),
    ddouble::unchecked(0x1.1eed8eff8d898p-29, -0x1.2aec959e14c06p-83),
    ddouble::unchecked(-0x1.93974a8c07c9dp-37, -0x1.05d6f8a2efd1fp-92),
    ddouble::unchecked(0x1.ae7f3e733b81fp-45, 0x1.1d8656b0ee8cbp-101),
    ddouble::unchecked(-0x1.6827863b97d97p-53, -0x1.eec01221a8b0bp-107),
    ddouble::unchecked(0x1.e542ba4020225p-62, 0x1.ea72b4afe3c2fp-120),
    ddouble::unchecked(-0x1.0ce396db7f853p-70, 0x1.aebcdbd20331cp-124),
    ddouble::unchecked(0x1.f2cf01972f578p-80, -0x1.9ada5fcc1ab14p-135),
    ddouble::unchecked(-0x1.88e85fc6a4e5ap-89, 0x1.71c37ebd16540p-143),
    ddouble::unchecked(0x1.0a18a2635085dp-98, 0x1.b9e2e28e1aa54p-153),
    ddouble::unchecked(-0x1.3932c5047d60ep-108, -0x1.832b7b530a627p-162),
    ddouble::unchecked(0x1.434d2e783f5bcp-118, 0x1.0b87b91be9affp-172),
    ddouble::unchecked(-0x1.2710231c0fd7ap-128, -0x1.3f8a2b4af9d6bp-184),
    ddouble::unchecked(0x1.df983290c2ca9p-139, 0x1.5835c6895393bp-194),
    ddouble::unchecked(-0x1.5d4acb9c0c3abp-149, 0x1.6ec2c8f5b13b2p-205),
    ddouble::unchecked(0x1.ca8ed42a12ae3p-160, 0x1.a07244abad2abp-224),
    ddouble::unchecked(-0x1.10af527530de8p-170, -0x1.b626c912ee5c8p-225),
    ddouble::unchecked(0x1.272b1b03fec6ap-181, 0x1.3f67cc9f9fdb8p-235),
    ddouble::unchecked(-0x1.240804f659510p-192, -0x1.8b291b93c9718p-246),
    ddouble::unchecked(0x1.091b406b6ff26p-203, 0x1.e973637973b18p-257),
    ddouble::unchecked(-0x1.bb36f6e12cd78p-215, -0x1.02f85029a29b0p-270),
    ddouble::unchecked(0x1.56457989358c9p-226, -0x1.e3792533eafc8p-282),
    ddouble::unchecked(-0x1.e9d8f6ed83eaap-238, 0x1.be25ac1066519p-293),
    ddouble::unchecked(0x1.45b77f9e98e12p-249, 0x1.e4b05119ccb1bp-303),
    ddouble::unchecked(-0x1.938cc661b03f6p-261, -0x1.c4da1977e56d6p-318),
};

} // namespace ddouble_detail

/** returns square of a. slightly faster than 'a * a' */
inline ddouble sqr(ddouble a)
{
	auto tmp = ddouble::mul(a.high(), a.high());
	return ddouble::sum_quick(tmp.high(),
	                          tmp.low() + std::ldexp(a.high() * a.low(), 1));
}

/*inline ddouble pow(ddouble a, int b)
{
    (high + low) ^ b = high ^ b + b * low * high ^ (b - 1)
}*/

inline ddouble inverse(ddouble a)
{
	double high = 1 / a.high();
	double low = (fma(-high, a.high(), 1.0) - high * a.low()) / a.high();
	return ddouble::sum_quick(high, low);
}

inline ddouble sqrt(ddouble a)
{
	auto high = std::sqrt(a.high());
	ddouble r = ldexp(high + a / high, -1);
	return r;
}

/** inverse square root '1/sqrt(a)' (probably faster than sqrt(a) itself) */
inline ddouble rec_sqrt(ddouble a)
{
	// TODO: would be nice to get the SSE instruction 'rsqrt' as a first
	// approximation (only gives ~11 correct bits, but very fast)
	double r = 1.0 / std::sqrt(a.high());
	return (0.5 * r) * (3 - a * ddouble::mul(r, r)); // one newton step
}

inline ddouble cbrt(ddouble a)
{
	auto high = std::cbrt(a.high());
	ddouble r = (2 * high + a / (ddouble::mul(high, high))) / 3;
	return r;
}

inline ddouble exp(ddouble a)
{
	using namespace ddouble_detail;

	// Idea: exp(k*log(2) + r) = exp(r/16)^16 * 2^k
	// (the repeated half/square-trick looses a bit of accuracy, which is fine)
	int k = (int)std::round(a.high() * (1.0 / M_LN2));
	a -= k * ddouble::ln2();
	a = ldexp(a, -4);

	// now it is |a| <= ln(2)/32 = 0.02166..., so the series converges quickly
	assert(std::abs(a.high()) < 0.022);
	ddouble r = taylor(coeffs_exp, a);

	// restore full answer
	r = sqr(sqr(sqr(sqr(r))));
	return ldexp(r, k);
}

inline ddouble log(ddouble a)
{
	double high = std::log(a.high());
	return high + a * exp(-ddouble(high)) - 1.0;
}

ddouble sin(ddouble a)
{
	using namespace ddouble_detail;

	// reduce mod pi
	int k = (int)std::round(a.high() * (1 / M_PI));
	a -= k * ddouble::pi();

	// now |a| <= pi/2 = 1.57...
	assert(std::abs(a.high()) < 1.58);

	// remaining a is small -> Taylor
	if (k & 1)
		return -a * taylor(coeffs_sin, a * a);
	else
		return a * taylor(coeffs_sin, a * a);
}

ddouble cos(ddouble a)
{
	using namespace ddouble_detail;

	// reduce mod pi
	int k = (int)std::round(a.high() * (1 / M_PI));
	a -= k * ddouble::pi();

	// now |a| <= pi/2 = 1.57...
	assert(std::abs(a.high()) < 1.58);

	// remaining a is small -> Taylor
	if (k & 1)
		return -taylor(coeffs_cos, a * a);
	else
		return taylor(coeffs_cos, a * a);
}

ddouble tan(ddouble a) { return sin(a) / cos(a); }
ddouble cot(ddouble a) { return cos(a) / sin(a); }
ddouble sec(ddouble a) { return inverse(cos(a)); }
ddouble csc(ddouble a) { return inverse(sin(a)); }

inline bool operator<(ddouble a, ddouble b)
{
	return a.high() < b.high() || (a.high() == b.high() && a.low() < b.low());
}
inline bool operator<=(ddouble a, ddouble b)
{
	return a.high() < b.high() || (a.high() == b.high() && a.low() <= b.low());
}
inline bool operator==(ddouble a, ddouble b)
{
	return a.high() == b.high() && a.low() == b.low();
}

inline bool operator<(ddouble a, double b)
{
	return a.high() < b || (a.high() == b && a.low() < 0);
}
inline bool operator<=(ddouble a, double b)
{
	return a.high() < b || (a.high() == b && a.low() <= 0);
}
inline bool operator==(ddouble a, double b)
{
	return a.high() == b && a.low() == 0;
}
inline bool operator<(double a, ddouble b)
{
	return a < b.high() || (a == b.high() && 0 < b.low());
}
inline bool operator<=(double a, ddouble b)
{
	return a < b.high() || (a == b.high() && 0 <= b.low());
}
inline bool operator==(double a, ddouble b) { return b == a; }

inline bool operator>(ddouble a, ddouble b) { return b < a; }
inline bool operator>=(ddouble a, ddouble b) { return b <= a; }
inline bool operator!=(ddouble a, ddouble b) { return !(a == b); }
inline bool operator>(ddouble a, double b) { return b < a; }
inline bool operator>=(ddouble a, double b) { return b <= a; }
inline bool operator!=(ddouble a, double b) { return !(a == b); }
inline bool operator>(double a, ddouble b) { return b < a; }
inline bool operator>=(double a, ddouble b) { return b <= a; }
inline bool operator!=(double a, ddouble b) { return !(a == b); }

template <> struct RingTraits<ddouble>
{
	/** potentially faster than full comparison */
	static bool is_zero(ddouble const &value) { return value.high() == 0.0; }
	static bool is_one(ddouble const &value) { return value == 1; }

	/** just set to false if signs dont make sense for T */
	static bool is_negative(ddouble const &value) { return value.high() < 0; }

	/** do we need parentheses if T occurs in a product/power? */
	static bool need_parens_product(ddouble const &) { return false; }
	static bool need_parens_power(ddouble const &value)
	{
		return is_negative(value);
	}
};

} // namespace chalk

namespace Eigen {

template <typename T> class NumTraits;

template <> struct NumTraits<chalk::ddouble>
{
	using Real = chalk::ddouble;
	using NonInteger = chalk::ddouble;
	using Nested = chalk::ddouble;
	using Literal = double; // is this correct?

	static inline Real epsilon() { return Real(ldexp(1.0, -107)); }
	static inline Real dummy_precision() { return Real(ldexp(1.0, -99)); }
	static inline int digits10() { return 29; }
	// missing: highest(), lowest()

	enum
	{
		IsComplex = 0,
		IsInteger = 0,
		IsSigned = 1,
		RequireInitialization = 1,

		// not tuned at all
		ReadCost = 1,
		AddCost = 3,
		MulCost = 10
	};
};

} // namespace Eigen
