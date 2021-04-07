#pragma once

#include "fmt/format.h"
#include "util/span.h"
#include <array>
#include <bitset>
#include <cassert>
#include <climits>
#include <optional>
#include <utility>
#include <vector>

namespace chalk {

/** a single term of a polynomial: a*x1^e1*...*xn^en */
template <typename R, size_t rank> struct Monomial
{
	R coefficient = {};
	std::array<int, rank> exponent = {};
};

template <typename R, size_t rank>
Monomial<R, rank> operator*(Monomial<R, rank> const &a,
                            Monomial<R, rank> const &b)
{
	std::array<int, rank> ex;
	for (size_t i = 0; i < rank; ++i)
		ex[i] = a.exponent[i] + b.exponent[i];
	return Monomial<R, rank>{a.coefficient * b.coefficient, ex};
}

template <typename R, size_t rank>
std::optional<Monomial<R, rank>> operator/(Monomial<R, rank> const &a,
                                           Monomial<R, rank> const &b)
{
	std::array<int, rank> ex;
	for (size_t i = 0; i < rank; ++i)
		if (a.exponent[i] >= b.exponent[i])
			ex[i] = a.exponent[i] - b.exponent[i];
		else
			return std::nullopt;
	return Monomial<R, rank>{a.coefficient / b.coefficient, ex};
}

template <typename R, size_t rank>
bool operator==(Monomial<R, rank> const &a, Monomial<R, rank> const &b)
{
	return a.coefficient == b.coefficient && a.exponent == b.exponent;
}

template <typename R, size_t rank>
bool order_degrevlex(Monomial<R, rank> const &a, Monomial<R, rank> const &b,
                     std::array<int, rank> const &ws)
{
	int totA = 0, totB = 0;
	for (size_t i = 0; i < rank; ++i)
	{
		totA += a.exponent[i] * ws[i];
		totB += b.exponent[i] * ws[i];
	}
	if (totA != totB)
		return totA > totB;

	for (int i = (int)rank - 1; i >= 0; --i)
		if (a.exponent[i] < b.exponent[i])
			return true;
		else if (a.exponent[i] > b.exponent[i])
			return false;
	return false;
}

} // namespace chalk
