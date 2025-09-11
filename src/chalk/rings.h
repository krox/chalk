#pragma once

#include "fmt/format.h"
#include "util/complex.h"
#include <complex>

namespace chalk {

using int128_t = __int128;

// abstract representation of "0" and "1", however that maps to any particular
// ring
struct zero_t
{};
struct one_t
{};
inline static constexpr zero_t zero;
inline static constexpr one_t one;

enum class format_type
{
	text,
	debug
	// TODO: latex
};

// concept to capture int{32,64,128}_t, float, double
template <class T>
concept builtin_number =
    std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
    std::is_same_v<T, int128_t> || std::is_same_v<T, float> ||
    std::is_same_v<T, double>;

bool is_zero(builtin_number auto const &value) { return value == 0; }
bool is_one(builtin_number auto const &value) { return value == 1; }
bool is_negative(builtin_number auto const &value) { return value < 0; }
bool need_parens_product(builtin_number auto const &) { return false; }
bool need_parens_power(builtin_number auto const &value)
{
	return is_negative(value);
}

template <class T> bool is_zero(std::complex<T> const &value)
{
	return value == T(0);
}
template <class T> bool is_one(std::complex<T> const &value)
{
	return value == T(1);
}
template <class T> bool is_negative(std::complex<T> const &) { return false; }
template <class T> bool need_parens_product(std::complex<T> const &)
{
	return true;
}
template <class T> bool need_parens_power(std::complex<T> const &)
{
	return true;
}

template <class T> bool is_zero(util::complex<T> const &value)
{
	return value == T(0);
}
template <class T> bool is_one(util::complex<T> const &value)
{
	return value == T(1);
}
template <class T> bool is_negative(util::complex<T> const &) { return false; }
template <class T> bool need_parens_product(util::complex<T> const &)
{
	return true;
}
template <class T> bool need_parens_power(util::complex<T> const &)
{
	return true;
}

template <typename T> struct Scalar
{
	T value;

	Scalar() = default;
	explicit Scalar(T v) : value(std::move(v)) {}
};
/*
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
        return need_parens_product(value.value);
    }
    static bool need_parens_power(Scalar<T> const &value)
    {
        return need_parens_power(value.value);
    }
};*/

using std::pow;
using std::sqrt;

// helper to pass overloaded functions as templated parameters.
// for details see: https://florianjw.de/en/passing_overloaded_functions.html
#define CHALK_LIFT(...)                                                        \
	[](auto &&...args) -> decltype(auto) {                                     \
		return __VA_ARGS__(std::forward<decltype(args)>(args)...);             \
	}

} // namespace chalk
