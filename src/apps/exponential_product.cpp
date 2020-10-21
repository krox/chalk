#include "chalk/bch.h"
#include "chalk/fraction.h"
#include "chalk/free_group.h"
#include "chalk/group_ring.h"
#include "chalk/ideal.h"
#include "chalk/integer.h"
#include "chalk/polynomial.h"
#include "chalk/sparse_polynomial.h"
#include <fmt/format.h>
using namespace chalk;

using Algebra = FreeAlgebra<SparsePolynomial<Rational, 8>>;

Algebra exp(Algebra const &a)
{
	// assert(a.max_order() < 100);
	auto inc = Algebra(1);
	Algebra r = Algebra(0);
	for (int i = 1; !(inc == 0); ++i)
	{
		r += inc;
		inc *= a;
		inc /= i;
	}
	return r;
}

Algebra getEpsOrder(Algebra const &a, int k)
{
	return map_coefficients(a, [&](SparsePolynomial<Rational, 8> const &poly) {
		return get_coefficient(poly, 0, k);
	});
}

auto rationalToDouble(Rational const &x)
{
	return (double)x.num() / (double)x.denom();
};

int main()
{
	int max_order = 4;
	auto R = PolynomialRing<Rational, 8>();
	auto eps = R.generator("eps", max_order);
	auto alpha = R.generator("alpha");
	auto beta = R.generator("beta");
	auto gamma = R.generator("gamma");
	auto delta = R.generator("delta");
	auto rho = R.generator("rho");
	auto lambda = R.generator("lambda");
	auto sigma = R.generator("sigma");
	delta = R("1/2-beta");
	rho = R("1-2*(alpha+gamma)");
	auto double_ring = PolynomialRing<double, 8>(R.var_names());

	auto a = eps * Algebra::generator(0);
	auto b = eps * Algebra::generator(1);

	auto error = exp(alpha * a);
	error *= exp(beta * b);
	error *= exp(gamma * a);
	error *= exp(delta * b);
	error *= exp(rho * a);
	error *= exp(delta * b);
	error *= exp(gamma * a);
	error *= exp(beta * b);
	error *= exp(alpha * a);

	error -= exp(a + b);
	fmt::print("{}\n", error);

	std::vector<SparsePolynomial<Rational, 8>> cond_list;
	for (int order = 0; order <= max_order; ++order)
	{
		auto cond = getEpsOrder(error, order);
		fmt::print("---- conditions (eps^{}) ----\n", order);
		for (auto &[poly, word] : cond.terms())
		{
			fmt::print("{}:   {}\n", word, poly);
			cond_list.push_back(poly);
		}
	}

	fmt::print("---- resulting ideal of coefficients ----\n");
	// cond_list.push_back(R("alpha-3/10"));
	auto ideal = Ideal(cond_list);
	dump(ideal);

	analyze_variety(ideal.change_ring(&double_ring, rationalToDouble), {},
	                true);
}
