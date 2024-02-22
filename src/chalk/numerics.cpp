#include "chalk/numerics.h"

#include "chalk/floating.h"
#include "fmt/format.h"
#include <algorithm>
#include <complex>
#include <stdexcept>

namespace chalk {

std::vector<double> roots(Polynomial<double> const &poly)
{
	// TODO: more sophisticated cutoffs. hardcoded 1.0e-12 is kinda stupid...
	// leading coefficient should determine the overall scale

	// analytical formulas for low degrees:
	if (poly.degree() == -1)
		throw std::runtime_error("cant factorize zero-polynomial");
	if (poly.degree() == 1)
		return {-poly[0] / poly[1]};
	if (poly.degree() == 2)
	{
		auto p = poly[1] / poly[2];
		auto q = poly[0] / poly[2];
		auto d = 0.25 * p * p - q;
		if (d < 0)
			return {};
		d = std::sqrt(d);

		// roots are very close together -> just return one
		if (std::abs(d) < 1.0e-12)
			return {-0.5 * p};

		// solutions are now -p/2 +- d. But it is numerically more precise to
		// determine only one root this way, and the other by Vieta lemma
		if (p > 0)
		{
			auto x = -0.5 * p - d;
			return {x, q / x};
		}
		else
		{
			auto x = -0.5 * p + d;
			return {q / x, x};
		}
	}

	std::vector<util::complex<double>> coeffs;
	coeffs.reserve(poly.degree() + 1);
	for (auto &c : poly.coefficients())
		coeffs.push_back(util::complex<double>(c));
	double eps = 1e-12;
	auto rc = rootsDurandKerner(
	    Polynomial<util::complex<double>>(std::move(coeffs)), eps);
	std::vector<double> r;
	for (auto &x : rc)
		if (std::abs(x.imag()) < eps)
			r.push_back(x.real());
	std::sort(r.begin(), r.end());

	return r;
}

LegendreRule computeLegendreQuadrature(int n)
{
	assert(n >= 1 && n % 2 == 1);
	LegendreRule rule;
	rule.xs.resize(n / 2);
	rule.ws.resize(n / 2);

	// compute the roots
	for (int k = 0; k < n / 2; ++k)
	{
		// starting guess (very precise, especially for large n)
		double guess = std::sin(M_PI * 2. * k / (2. * n + 1.));
		guess *= 1. - (n - 1.) / (8. * n * n * n); // + O(n^-4)

		// a few iterations of Newton method
		FloatingOctuple x = FloatingOctuple(guess);
		FloatingOctuple fd;
		for (int iter = 1; iter < 10; ++iter)
		{
			auto [f, fd_] = legendrePolynomial(n, x);
			fd = std::move(fd_);
			x -= f / fd;
		}
		rule.xs[k] = x.to_double();
		rule.ws[k] = (2 / ((1 - x * x) * fd * fd)).to_double();
	}
	return rule;
}

} // namespace chalk
