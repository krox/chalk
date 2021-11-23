#pragma once

#include "fmt/format.h"
#include <complex>

namespace chalk {

using int128_t = __int128;

template <typename T> struct RingTraits;

// implementation of RingTraits for int/float/Integer/Floating/...
template <typename T> struct RingTraitsSimple
{
	// potentially faster than full comparison
	static bool isZero(T const &value) { return value == 0; }
	static bool isOne(T const &value) { return value == 1; }

	// hints for nice formatting of compound types
	static bool isNegative(T const &value) { return value < 0; }
	static bool needParensProduct(T const &) { return false; }
	static bool needParensPower(T const &value) { return isNegative(value); }
};

template <> struct RingTraits<int32_t> : RingTraitsSimple<int32_t>
{};
template <> struct RingTraits<int64_t> : RingTraitsSimple<int64_t>
{};
template <> struct RingTraits<int128_t> : RingTraitsSimple<int128_t>
{};
template <> struct RingTraits<float> : RingTraitsSimple<float>
{};
template <> struct RingTraits<double> : RingTraitsSimple<double>
{};

template <typename T> struct RingTraits<std::complex<T>>
{
	static bool isZero(std::complex<T> const &value) { return value == T(0); }
	static bool isOne(std::complex<T> const &value) { return value == T(1); }
	static bool isNegative(std::complex<T> const &value) { return false; }
	static bool needParensProduct(std::complex<T> const &) { return true; }
	static bool needParensPower(std::complex<T> const &) { return true; }
};

template <typename T> bool isZero(T const &value)
{
	return RingTraits<T>::isZero(value);
}
template <typename T> bool isOne(T const &value)
{
	return RingTraits<T>::isOne(value);
}
template <typename T> bool isNegative(T const &value)
{
	return RingTraits<T>::isNegative(value);
}
template <typename T> bool needParensProduct(T const &value)
{
	return RingTraits<T>::needParensProduct(value);
}
template <typename T> bool needParensPower(T const &value)
{
	return RingTraits<T>::needParensPower(value);
}

template <typename T> struct Scalar
{
	T value;

	Scalar() = default;
	explicit Scalar(T v) : value(std::move(v)) {}
};

template <typename T> struct RingTraits<Scalar<T>>
{
	static bool isZero(Scalar<T> const &value) { return isZero(value.value); }
	static bool isOne(Scalar<T> const &value) { return isOne(value.value); }
	static bool isNegative(Scalar<T> const &value)
	{
		return isNegative(value.value);
	}
	static bool needParensProduct(Scalar<T> const &value)
	{
		return needParensProduct(value.value);
	}
	static bool needParensPower(Scalar<T> const &value)
	{
		return needParensPower(value.value);
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
