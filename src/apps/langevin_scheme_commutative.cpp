#include "chalk/bch.h"
#include "chalk/covariant.h"
#include "chalk/fraction.h"
#include "chalk/ideal.h"
#include "chalk/sparse_polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

static constexpr int max_order = 4; // orders of epsilon we want everything in
#define RANK 32 // enough space for all parameters (plus one for 'eps')

using R = SparsePolynomial<Rational, RANK>;
auto ring = PolynomialRing<Rational, RANK>();
R seps = ring.generator("s", max_order * 2);
R eps = seps * seps;
R cA = ring.generator("cA") * 0;
using Cov = Covariant<R>;

/** extract the coefficient of eps^k */
auto getEpsOrder(Cov const &a, int k)
{
	return map_coefficients(
	    a, [&](R const &poly) { return get_coefficient(poly, 0, 2 * k); });
}
auto getSepsOrder(Cov const &a, int k)
{
	return map_coefficients(
	    a, [&](R const &poly) { return get_coefficient(poly, 0, k); });
}

/** diff(A*exp(-S))*exp(S) */
auto myDiff(Cov const &a)
{
	// NOTE: '*' not abelian in this case due to index-naming
	return diff(a) - a * Cov(Indexed("S", 1));
}

auto rationalToDouble(Rational const &x)
{
	return (double)x.num() / (double)x.denom();
};

/** Same as 'a*b', but faster because it knows about finite order of epsilon */
Cov myMul(Cov const &a, Cov const &b)
{
	// split into orders of epsilon
	std::array<Cov, max_order * 2 + 1> A, B, C;
	for (int i = 0; i <= 2 * max_order; ++i)
	{
		A[i] = getSepsOrder(a, i);
		B[i] = getSepsOrder(b, i);
		C[i] = Cov(0);
	}

	// multiply
	for (int i = 0; i <= 2 * max_order; ++i)
		for (int j = 0; j <= i; ++j)
			C[i] += A[j] * B[i - j];

	// collect into one again
	auto c = Cov(0);
	for (int i = 0; i <= 2 * max_order; ++i)
		c += C[i] * pow(seps, i);
	return c;
}

void add_term(Cov &a, std::string const &name, Cov term, int order)
{
	while (order-- > 0)
		term *= seps;
	term = simplify_lie_commutative(term, "S");

	fmt::print("{}      : {}\n", name, term);

	a += ring(name) * term;
}

int main()
{
	/** NOTE:
	 *  use something like 'ring.generator(..., INT_MAX, 20)' to increase
	 *  the weight of a variable (for term-ordering), to force the gröbner-
	 *  algorithm to solve for that variable first.
	 */

	/** the forces */
	auto eta1 = Cov(Indexed("eta1", 1));
	auto eta2 = Cov(Indexed("eta2", 1));
	auto eta3 = Cov(Indexed("eta3", 1));
	auto S0 = Cov(Indexed("S", 1));

	fmt::print("---------- Terms of the scheme ----------\n");

	auto f1 = Cov(0);
	add_term(f1, "k1", S0, 2);
	add_term(f1, "k2", eta1, 1);
	// add_term(f1, "k3", eta2, 1);
	// add_term(f1, "k3b", eta3, 1);
	auto S1 = taylor(S0, -f1, 2 * max_order);

	// even terms
	auto f2 = Cov(0);
	add_term(f2, "k4", S0, 2);
	add_term(f2, "k5", S1, 2);
	add_term(f2, "k6", eta1, 1);
	add_term(f2, "k7", eta2, 1);
	// add_term(f2, "k7b", eta3, 1);
	auto S2 = taylor(S0, -f2, 2 * max_order);

	// even terms
	auto f3 = Cov(0);
	add_term(f3, "k8", S0, 2);
	add_term(f3, "k9", S1, 2); // remove this for "Torrero-style" at order 3
	add_term(f3, "k10", S2, 2);
	add_term(f3, "k11", eta1, 1);
	add_term(f3, "k12", eta2, 1);
	auto S3 = taylor(S0, -f3, 2 * max_order);

	auto f4 = Cov(0);
	add_term(f4, "k13", S0, 2);
	add_term(f4, "k14", S1, 2);
	add_term(f4, "k15", S2, 2);
	add_term(f4, "1-k15-k14-k13", S3, 2);
	add_term(f4, "k17", eta1, 1); // scale setting at order 4
	add_term(f4, "k18", eta2, 1); // scale setting at order 4
	// add_term(f4, "k18", eta2, 1); // scale setting at order 4
	auto f = f4;

	// full condition is
	//  (∇ <f> + 1/2 ∇∇<ff> + ...) e^-S = (T-1)e^-S
	fmt::print("---------- Expectation values ----------\n");
	auto full_condition = Cov(0);
	for (int k = 1; k <= max_order * 2; ++k)
	{
		// build ev = <f^k>
		auto ev = Cov(1);
		for (int i = 0; i < k; ++i)
			ev = myMul(ev, f); // same as 'ev *= f'
		ev = wick_contract(ev, "eta1", R(2));
		ev = wick_contract(ev, "eta2", R(2));
		ev = wick_contract(ev, "eta3", R(2));

		fmt::print("----- <f^{}> -----\n", k);
		fmt::print("{}\n", ev);

		// build ev = 1/k! ∇^k <f^k> and add it to the complete
		// condition
		for (int i = 1; i <= k; ++i)
		{
			ev = myDiff(ev);
			ev = contract(ev, 0, k + 1 - i);
			ev = ev / i;
		}
		full_condition += ev;
	}

	full_condition = simplify_lie_commutative(full_condition, "S");

	// collect order conditions for all epsilon
	// these are 'T^(k) e^-S = 0'
	std::vector<R> cond_list;
	for (int k = 0; k <= max_order; ++k)
	{
		auto tmp = getEpsOrder(full_condition, k);
		if (tmp == 0)
			continue;

		fmt::print("---------- T^({})e^-S ----------\n", k);
		dump_summary(tmp);
		for (auto term : tmp.terms())
			cond_list.push_back(term.coefficient);
	}

	/*std::map<std::string, std::pair<double, double>> constraints = {
	    {"k1", {0.0, 1.0}}};*/
	std::map<std::string, std::pair<double, double>> constraints = {};
	// cond_list.push_back(ring("c11-1/2"));
	fmt::print("\ngeneral form (before reduce)\n");
	dump(cond_list);
	dump_singular(cond_list, "ideal.singular");

	reduce_partial(cond_list);
	/*fmt::print("---------- Dependence of conditions on coefficients "
	           "----------\n");
	for (int i = 2; i < RANK; ++i)
	{
	    fmt::print("----- {} -----\n", ring.var_names()[i]);
	    for (size_t j = 0; j < cond_list.size(); ++j)
	    {
	        auto tmp = diff(cond_list[j], i);
	        if (!(tmp == 0))
	            fmt::print("df{}/d{} = {}\n", j, ring.var_names()[i], tmp);
	    }
	}*/
	fmt::print("\ngeneral form (after reduce)\n");
	dump(cond_list);
	/*fmt::print("\ngeneral form (after gröbner)\n");
	groebner(cond_list);
	dump(cond_list);

	auto double_ring = PolynomialRing<double, RANK>(ring.var_names());
	analyze_variety(change_ring(cond_list, &double_ring, rationalToDouble),
	                constraints, true);*/
}
