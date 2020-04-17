#include "chalk/numerics.h"

#include "fmt/format.h"
#include <algorithm>
#include <complex>
#include <stdexcept>

namespace chalk {

/** polynomial root finding using Durand–Kerner method */
std::vector<std::complex<double>> roots(Polynomial<std::complex<double>> poly)
{
	if (poly.degree() == 0)
		throw std::runtime_error("tried to compute zeros of a 0-polynomial");

	// normalize the polynomial to avoid some weirdnesses
	poly /= poly.coefficients().back();

	// staring points: distribute over unit circle using golden ratio
	auto phi = std::complex<double>(-0.73736887807831990, 0.67549029426152364);
	std::vector<std::complex<double>> r;
	r.reserve(poly.degree());
	r.push_back(phi);
	for (int i = 1; i < poly.degree(); ++i)
		r.push_back(r.back() * phi);

	// refine using newton method
	int last_change = 0;
	for (int iter = 0; iter <= last_change + 3; ++iter)
	{
		if (iter > 100)
			throw std::runtime_error(
			    "polynomial root finding did not converge (max iter)");
		/*fmt::print("iter = {}\n", iter);
		for (size_t i = 0; i < r.size(); ++i)
		    fmt::print("  {} + {} i\n", r[i].real(), r[i].imag());*/

		for (size_t i = 0; i < r.size(); ++i)
		{
			std::complex<double> tmp = 1;
			for (size_t j = 0; j < r.size(); ++j)
				if (j != i)
					tmp *= r[i] - r[j];
			auto step = poly(r[i]) / tmp;
			if (std::abs(step) > 1.0e-12)
				last_change = iter;
			r[i] = r[i] - poly(r[i]) / tmp;
		}
	}

	// another sanity check (also catches NaN)
	for (auto &x : r)
		if (!(std::abs(poly(x)) < 1.0e-12))
			throw std::runtime_error(
			    "polynomial root finding did not converge (weird result)");
	return r;
}

std::vector<double> roots(Polynomial<double> const &poly)
{
	std::vector<std::complex<double>> coeffs;
	coeffs.reserve(poly.degree() + 1);
	for (auto &c : poly.coefficients())
		coeffs.push_back(c);
	auto rc = roots(Polynomial<std::complex<double>>(std::move(coeffs)));
	std::vector<double> r;
	for (auto &x : rc)
		if (std::abs(x.imag()) < 1.0e-12)
			r.push_back(x.real());
	std::sort(r.begin(), r.end());
	return r;
}
} // namespace chalk
