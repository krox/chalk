#include "chalk/bch.h"
#include "chalk/covariant.h"
#include "chalk/fraction.h"
#include "chalk/ideal.h"
#include "chalk/sparse_polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

/**
 * For the 1D case, we represent the derivatives as S', S'', ..., i.e. without
 * any indexes (even though we use the 'Indexed' class).
 * taylor() -> taylor_1d()
 * diff() -> diff_1d()
 * noise term(s) 'eta' are represented as polynomial vars, instead of 'indexed'
 */

static constexpr int max_order = 4; // orders of epsilon we want everything in
#define RANK 32 // enough space for all parameters (plus one for 'eps')

using R = SparsePolynomial<Rational, RANK>;
auto ring = PolynomialRing<Rational, RANK>();
R seps = ring.generator("s", max_order * 2);
R eps = seps * seps;
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
	return diff_1d(a) - a * Cov(Indexed("S'"));
}

/** 1D 'wick contract' (just expectation values of normal distributions) */
Cov my_wick_contract(Cov const &a, std::string const &eta)
{
	// NOTE: we assume Variance==2
	return map_coefficients(a, [&eta](R const &poly) {
		R r = ring((int)0);
		std::vector<R> terms = get_coefficients(poly, ring.var_id(eta));
		int coeff = 1;
		for (int i = 0; i < (int)terms.size(); i += 2)
		{
			r += terms[i] * coeff;
			coeff *= 2; // this is the Variance
			coeff *= i + 1;
		}
		return r;
	});
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
	// term = simplify_lie_commutative(term, "S");

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
	/*ring("k1");
	ring("k4");
	ring("k5");
	ring("k8");
	ring("k9");
	ring("k10");
	ring("k13");
	ring("k14");
	ring("k15");
	ring("k16");

	ring("k2");
	ring("k3");
	ring("k6");
	ring("k7");
	ring("k11");
	ring("k12");*/
	// ring("k17");
	/** the forces */
	auto eta1 = ring("eta1");
	auto eta2 = ring("eta2");
	auto eta3 = ring("eta3");
	auto S0 = Cov(Indexed("S'"));

	fmt::print("---------- Terms of the scheme ----------\n");

	auto f1 = Cov(0);
	add_term(f1, "k1", S0, 2);
	add_term(f1, "k2", Cov(eta1), 1);
	add_term(f1, "k3", Cov(eta2), 1);
	// add_term(f1, "k3b", eta3, 1);
	auto S1 = taylor_1d(S0, -f1, 2 * max_order);

	// even terms
	auto f2 = Cov(0);
	add_term(f2, "k4", S0, 2); // this one was negative in early results
	add_term(f2, "k5", S1, 2);
	add_term(f2, "k6", Cov(eta1), 1);
	add_term(f2, "k7", Cov(eta2), 1);
	// add_term(f2, "k7b", eta3, 1);
	auto S2 = taylor_1d(S0, -f2, 2 * max_order);

	// even terms
	auto f3 = Cov(0);
	add_term(f3, "k8", S0, 2);
	add_term(f3, "k9", S1, 2); // remove this for "Torrero-style" at order 3
	add_term(f3, "k10", S2, 2);
	add_term(f3, "k11", Cov(eta1), 1);
	add_term(f3, "k12", Cov(eta2), 1);
	auto S3 = taylor_1d(S0, -f3, 2 * max_order);

	auto f4 = Cov(0);
	add_term(f4, "k13", S0, 2);
	add_term(f4, "k14", S1, 2);
	add_term(f4, "k15", S2, 2);
	add_term(f4, "k16", S3, 2);
	add_term(f4, "1", Cov(eta1), 1); // scale setting at order 4
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
		ev = my_wick_contract(ev, "eta1");
		ev = my_wick_contract(ev, "eta2");
		ev = my_wick_contract(ev, "eta3");

		fmt::print("----- <f^{}> -----\n", k);
		fmt::print("{}\n", ev);

		// build ev = 1/k! ∇^k <f^k> and add it to the complete
		// condition
		for (int i = 1; i <= k; ++i)
		{
			ev = myDiff(ev);
			ev = ev / i;
		}
		full_condition += ev;
	}

	// full_condition = simplify_lie_commutative(full_condition, "S");

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
