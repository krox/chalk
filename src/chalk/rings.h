#pragma once

#include "fmt/format.h"
#include <complex>

namespace chalk {

using int128_t = __int128;

template <typename T> struct RingTraits;

/** implementation of RingTraits for int/float/Integer/Floating/... */
template <typename T> struct RingTraitsSimple
{
	/** potentially faster than full comparison */
	static bool is_zero(T const &value) { return value == 0; }
	static bool is_one(T const &value) { return value == 1; }

	/** just set to false if signs dont make sense for T */
	static bool is_negative(T const &value) { return value < 0; }

	/** do we need parentheses if T occurs in a product/power? */
	static bool need_parens_product(T const &) { return false; }
	static bool need_parens_power(T const &value) { return is_negative(value); }
};

template <> struct RingTraits<int32_t> : RingTraitsSimple<int32_t>
{};
template <> struct RingTraits<int64_t> : RingTraitsSimple<int32_t>
{};
template <> struct RingTraits<int128_t> : RingTraitsSimple<int128_t>
{};
template <> struct RingTraits<float> : RingTraitsSimple<float>
{};
template <> struct RingTraits<double> : RingTraitsSimple<double>
{};

template <typename T> struct RingTraits<std::complex<T>>
{
	static bool is_zero(std::complex<T> const &value) { return value == T(0); }
	static bool is_one(std::complex<T> const &value) { return value == T(1); }
	static bool is_negative(std::complex<T> const &value) { return false; }
	static bool need_parens_product(std::complex<T> const &) { return true; }
	static bool need_parens_power(std::complex<T> const &) { return true; }
};

template <typename T> bool is_zero(T const &value)
{
	return RingTraits<T>::is_zero(value);
}
template <typename T> bool is_one(T const &value)
{
	return RingTraits<T>::is_one(value);
}
template <typename T> bool is_negative(T const &value)
{
	return RingTraits<T>::is_negative(value);
}
template <typename T> bool need_parens_product(T const &value)
{
	return RingTraits<T>::need_parens_product(value);
}
template <typename T> bool need_parens_power(T const &value)
{
	return RingTraits<T>::need_parens_power(value);
}

template <typename T> struct Scalar
{
	T value;

	Scalar() = default;
	explicit Scalar(T v) : value(std::move(v)) {}
};

template <typename T> struct RingTraits<Scalar<T>>
{
	static bool is_zero(Scalar<T> const &value) { return is_zero(value.value); }
	static bool is_one(Scalar<T> const &value) { return is_one(value.value); }
	static bool is_negative(Scalar<T> const &value)
	{
		return is_negative(value.value);
	}
	static bool need_parens_product(Scalar<T> const &value)
	{
		return needs_parens_product(value.value);
	}
	static bool need_parens_power(Scalar<T> const &value)
	{
		return need_paren_power(value.value);
	}
};

using std::pow;
using std::sqrt;

// helper to pass overloaded functions as templated parameters.
// for details see: https://florianjw.de/en/passing_overloaded_functions.html
#define CHALK_LIFT(...)                                                        \
	[](auto &&... args) -> decltype(auto) {                                    \
		return __VA_ARGS__(std::forward<decltype(args)>(args)...);             \
	}

} // namespace chalk
