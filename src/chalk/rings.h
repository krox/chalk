#pragma once

#include "fmt/format.h"

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

} // namespace chalk
