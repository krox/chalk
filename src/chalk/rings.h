#pragma once

#include "fmt/format.h"
#include "util/complex.h"
#include <complex>
#include <concepts>
#include <cstdint>

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

template <class T> struct Scalar
{
	T value;

	Scalar() = default;
	explicit Scalar(T v) : value(std::move(v)) {}
};

// A type of ring *elements* (Integer, Fraction<R>, Polynomial<R>, …), not
// the ring as a structure. Customization is by operators / ADL free
// functions (is_zero, is_one). Mixed comparison with int is not required:
// T(int) together with is_zero/is_one is enough. Compound assignment may
// return anything.
template <class T>
concept ring_element = requires(T a, T b, int n) {
	T{};
	T(n);
	{ a + b } -> std::convertible_to<T>;
	{ a - b } -> std::convertible_to<T>;
	{ (a * b) } -> std::convertible_to<T>;
	{ -a } -> std::convertible_to<T>;
	a += b;
	a -= b;
	a *= b;
	{ a == b } -> std::convertible_to<bool>;
	{ is_zero(a) } -> std::convertible_to<bool>;
	{ is_one(a) } -> std::convertible_to<bool>;
};

// Visible here so gcd_domain<int64_t> does not depend on include order.
// Definitions live in numtheory.h.
constexpr int64_t gcd(int64_t a, int64_t b) noexcept;
constexpr void remove_common_factor(int64_t &a, int64_t &b) noexcept;

// What Fraction<R> actually needs to reduce a/b. gcd's return type must be T;
// remove_common_factor may return anything (Integer returns the gcd, int64_t
// returns void).
template <class T>
concept gcd_domain = ring_element<T> && requires(T a, T b) {
	{ gcd(a, b) } -> std::convertible_to<T>;
	remove_common_factor(a, b);
};

using std::pow;
using std::sqrt;

// helper to pass overloaded functions as templated parameters.
// for details see: https://florianjw.de/en/passing_overloaded_functions.html
#define CHALK_LIFT(...)                                                        \
	[](auto &&...args) -> decltype(auto) {                                     \
		return __VA_ARGS__(std::forward<decltype(args)>(args)...);             \
	}

} // namespace chalk
