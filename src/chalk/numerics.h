#ifndef CHALK_NUMERICS_H
#define CHALK_NUMERICS_H

#include "chalk/polynomial.h"
#include <vector>

namespace chalk {
std::vector<double> roots(Polynomial<double> const &poly);

/**
 * jacobi polynomial (d/dx)^k P_n^(α,β)(x)
 */
template <typename T, typename T2 = T>
T jacobiPolynomial(T2 alpha, T2 beta, int n, int k, T x)
{
	assert(n >= 0 && k >= 0);

	// derivatives
	if (k > n)
		return T(0);
	if (k != 0)
	{
		T factor = T(1);
		for (int i = 1; i <= k; ++i)
			factor *= (alpha + beta + n + i) / 2;
		return factor * jacobiPolynomial(alpha + k, beta + k, n - k, 0, x);
	}

	// jacobi without derivatives by simple recursion formula
	if (n == 0)
		return T(1);
	auto y0 = T(1);
	auto y1 = (alpha + 1) + (alpha + beta + 2) / 2 * (x - 1);
	for (int m = 2; m <= n; ++m)
	{
		auto y2 = (2 * m + alpha + beta - 1) *
		          ((2 * m + alpha + beta) * (2 * m + alpha + beta - 2) * x +
		           alpha * alpha - beta * beta) *
		          y1;
		y2 -=
		    2 * (m + alpha - 1) * (m + beta - 1) * (2 * m + alpha + beta) * y0;

		y2 /= 2 * m * (m + alpha + beta) * (2 * m + alpha + beta - 2);
		y0 = y1;
		y1 = y2;
	}

	return y1;
}

/** value and derivative of the n'th Legendre polynomial */
template <typename T> std::pair<T, T> legendrePolynomial(int n, T x)
{
	// "Bonnet's recursion formula"
	auto y0 = T(0);
	auto y1 = T(1);
	auto d1 = T(0);
	for (int m = 1; m <= n; ++m)
	{
		d1 = m * y1 + x * d1;
		auto y2 = ((2 * m - 1) * x * y1 - (m - 1) * y0) / m;
		y0 = y1;
		y1 = y2;
	}

	return {y1, d1};
}

struct LegendreRule
{
	// odd N, normalized to integration on [-1, 1]
	// N = 2*xs.size() - 1
	// xs[0] = 0
	std::vector<double> xs;
	std::vector<double> ws;
};

/** compute points/weights for legendre quadrature */
LegendreRule computeLegendreQuadrature(int n);

} // namespace chalk

#endif
