#pragma once

#include "chalk/polynomial.h"
#include "util/complex.h"
#include <vector>

namespace chalk {
std::vector<double> roots(Polynomial<double> const &poly);

// polynomial root finding using Durand–Kerner method
// should work for float/double/Floating/...
template <typename Real>
std::vector<util::complex<Real>>
rootsDurandKerner(Polynomial<util::complex<Real>> poly, Real const &eps)
{
	using Complex = util::complex<Real>;

	if (poly.degree() == 0)
		throw std::runtime_error("tried to compute zeros of a 0-polynomial");

	// normalize the polynomial to avoid some weirdnesses
	poly /= poly.coefficients().back();

	// staring points: distribute over unit circle using golden ratio
	auto phi = Complex(Real(-0.73736887807831990), Real(0.67549029426152364));
	std::vector<Complex> r;
	r.reserve(poly.degree());
	r.push_back(phi);
	for (int i = 1; i < poly.degree(); ++i)
		r.push_back(r.back() * phi);

	// refine using Newton method
	int last_change = 0;
	for (int iter = 0; iter <= last_change + 3; ++iter)
	{
		if (iter > 100)
			throw std::runtime_error(
			    "polynomial root finding did not converge (max iter)");

		for (size_t i = 0; i < r.size(); ++i)
		{
			auto tmp = Complex(Real(1.0));
			for (size_t j = 0; j < r.size(); ++j)
				if (j != i)
					tmp *= r[i] - r[j];
			auto step = poly(r[i]) / tmp;
			if (abs(step) > eps)
				last_change = iter;
			r[i] = r[i] - poly(r[i]) / tmp;
		}
	}

	// another sanity check (also catches NaN)
	for (auto &x : r)
		if (!(abs(poly(x)) < 1.0e-12))
			throw std::runtime_error(
			    "polynomial root finding did not converge (weird result)");
	return r;
}

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
