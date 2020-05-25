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
T jacobi(T2 alpha, T2 beta, int n, int k, T x)
{
	assert(n >= 0 && k >= 0);

	// derivatives
	if (k > n)
		return T(0);
	if (k != 0)
	{
		T factor = 1;
		for (int i = 1; i <= k; ++i)
			factor *= (alpha + beta + n + i) / 2;
		return factor * jacobi(alpha + k, beta + k, n - k, 0, x);
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

} // namespace chalk

#endif
