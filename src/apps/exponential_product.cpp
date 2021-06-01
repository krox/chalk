#include "chalk/bch.h"
#include "chalk/fraction.h"
#include "chalk/free_group.h"
#include "chalk/group_ring.h"
#include "chalk/ideal.h"
#include "chalk/integer.h"
#include "chalk/polynomial.h"
#include "chalk/series.h"
#include "chalk/sparse_polynomial.h"
#include <fmt/format.h>
using namespace chalk;

using Poly = SparsePolynomial<Rational, 8>;
using Algebra = FreeAlgebra<Poly>;
using Term = Series<Algebra, 4, 'e'>;

Term exp(Term const &a)
{
	// assert(a.max_order() < 100);
	auto inc = Term(1);
	Term r = Term(0);
	for (int i = 1; !(inc == 0); ++i)
	{
		r += inc;
		inc *= a;
		inc /= i;
	}
	return r;
}

int main()
{
	/**
	  - we only consider symmetric schemes
	  - 3 steps, order 2 -> unique solution, leapfrog
	  - 7 steps, order 4 -> unique solution, negative coeffs
	                        (Eq. 42 in arxiv.org/pdf/math-ph/0506007.pdf)
	  - 11 steps, order 6 -> no solutions
	*/
	auto R = PolynomialRing<Rational, 8>();
	auto eps = Term::generator() * Algebra(R(1));
	auto alpha = Algebra(R.generator("alpha"));
	auto beta = Algebra(R.generator("beta"));
	auto gamma = Algebra(R.generator("gamma"));
	auto delta = Algebra(R.generator("delta"));

	auto a = eps * Algebra::generator(0);
	auto b = eps * Algebra::generator(1);

	auto error = exp(alpha * a);
	error *= exp(beta * b);
	error *= exp(gamma * a);
	error *= exp(delta * b);
	error *= exp(gamma * a);
	error *= exp(beta * b);
	error *= exp(alpha * a);

	error -= exp(a + b);
	// fmt::print("error = {}\n", error);

	std::vector<Poly> ideal;
	for (auto const &cond : error.coefficients())
		for (auto &[poly, word] : cond.terms())
		{
			fmt::print("{}:   {}\n", word, poly);
			ideal.push_back(poly);
		}

	reduce(ideal);

	fmt::print("---- resulting ideal of coefficients ----\n");
	dump(ideal);

	analyze_variety(change_ring<double>(ideal), {}, true);
}
