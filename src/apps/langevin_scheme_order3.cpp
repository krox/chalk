#include "chalk/bch.h"
#include "chalk/covariant.h"
#include "chalk/fraction.h"
#include "chalk/ideal.h"
#include "chalk/sparse_polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

static constexpr int max_order = 3; // orders of epsilon we want everything in
#define RANK 32 // enough space for all parameters (plus two for 'cA' and 'eps')

using R = SparsePolynomial<Rational, RANK>;
auto ring = PolynomialRing<Rational, RANK>();
R seps = ring.generator("s", max_order * 2);
R eps = seps * seps;
R cA = ring.generator("cA");
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

/** extract coefficient of cA^k */
auto getcAOrder(Cov const &a, int k)
{
	return map_coefficients(
	    a, [&](R const &poly) { return get_coefficient(poly, 1, k); });
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

/** commutator c=[a,b]. all three interpreted as  'iT^a f_a' */
Cov comm(Cov const &a, Cov const &b)
{
	// need exactly one open index
	for (auto &[_, indexed] : a.terms())
		assert(indexed.rank() == 1);
	for (auto &[_, indexed] : b.terms())
		assert(indexed.rank() == 1);
	auto c = -Cov(Indexed("f", 0, 1, 2)) * a * b; // order is important here!
	c = contract(c, 1, 3);
	c = contract(c, 1, 2);
	return c;
}
/** commutator [a,[b,c]] */
Cov comm(Cov const &a, Cov const &b, Cov const &c)
{
	return comm(a, comm(b, c));
}
Cov comm(Cov const &a, Cov const &b, Cov const &c, Cov const &d)
{
	return comm(a, comm(b, comm(c, d)));
}

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
	term = simplify_lie(term, "S", cA);
	a += ring(name) * term;
}

int main()
{
	/** NOTE:
	 *  use something like 'ring.generator(..., INT_MAX, 20)' to increase
	 *  the weight of a variable (for term-ordering), to force the gröbner-
	 *  algorithm to solve for that variable first.
	 */
	ring.generator("k5", INT_MAX, 100);

	/** the forces */
	auto eta = Cov(Indexed("eta", 1));
	auto S0 = Cov(Indexed("S", 1));

	// first step
	auto f1 = Cov(0);
	add_term(f1, "k1", S0, 2);
	add_term(f1, "k2", comm(eta, eta, S0), 4);
	// odd terms
	add_term(f1, "k3", eta, 1);
	add_term(f1, "k4", comm(eta, S0), 3);
	auto S1 = taylor(S0, -f1, 2 * max_order);

	// even terms
	auto f2 = Cov(0);
	add_term(f2, "k5", S0, 2);
	add_term(f2, "k6", S1, 2);
	add_term(f2, "k7", comm(eta, eta, S0), 4);
	// odd terms
	add_term(f2, "k8", eta, 1);
	add_term(f2, "k9", comm(eta, S0), 3);
	add_term(f2, "k10", comm(eta, S1), 3);
	auto S2 = taylor(S0, -f2, 2 * max_order);

	// even terms
	auto f3 = Cov(0);
	add_term(f3, "k11", S0, 2);
	add_term(f3, "k12", S1, 2);
	add_term(f3, "k13", S2, 2);
	add_term(f3, "k14", comm(S0, S1), 4);
	add_term(f3, "k15", comm(S0, S2), 4);
	add_term(f3, "k16", comm(S1, S2), 4);
	add_term(f3, "k17", comm(eta, eta, S0), 4);
	add_term(f3, "k18", comm(eta, eta, S1), 4);
	add_term(f3, "k19", comm(eta, eta, S2), 4);
	add_term(f3, "k20", comm(eta, eta, comm(eta, eta, S0)), 6);
	// odd terms
	add_term(f3, "1", eta, 1); // scale setting
	add_term(f3, "k21", comm(eta, S0), 3);
	add_term(f3, "k22", comm(eta, S1), 3);
	add_term(f3, "k23", comm(eta, S2), 3);
	add_term(f3, "k24", comm(eta, eta, eta, S0), 5); // does not contribute???
	add_term(f3, "k25", comm(eta, eta, eta, S1), 5);
	add_term(f3, "k26", comm(eta, eta, eta, S2), 5);
	add_term(f3, "k27", comm(S0, eta, S1), 5);
	// add_term(f3, "k28", comm(S0, eta, S2), 5); // k27,k28 are equivalent
	add_term(f3, "k29", comm(S1, eta, S2), 5); // but k29 is different
	// TODO: THERE MIGHT BE SOME eps^5/2 TERMS MISSING
	auto f = f3;

	// NOTE:
	//   -

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
		ev = wick_contract(ev, "eta", R(2));

		fmt::print("----- <f^{}> -----\n", k);
		fmt::print("{}\n", ev);

		// build ev = 1/k! ∇^k <f^k> and add it to the complete condition
		for (int i = 1; i <= k; ++i)
		{
			ev = myDiff(ev);
			ev = contract(ev, 0, k + 1 - i);
			ev = ev / i;
		}
		full_condition += ev;
	}

	full_condition = simplify_lie(full_condition, "S", cA);

	// collect order conditions for all epsilon
	// these are 'T^(k) e^-S = 0'
	std::vector<R> cond_list;
	for (int k = 0; k <= max_order; ++k)
	{
		auto tmp = getEpsOrder(full_condition, k);
		if (tmp == 0)
			continue;
		for (int l = 0; l <= max_order + 3; ++l)
		{
			auto tmp2 = getcAOrder(tmp, l);
			if (tmp2 == 0)
				continue;
			fmt::print("---------- T^({},{})e^-S ----------\n", k, l);
			dump_summary(tmp2);
			for (auto term : tmp2.terms())
				cond_list.push_back(term.coefficient);
		}
	}

	/*std::map<std::string, std::pair<double, double>> constraints = {
	    {"k1", {0.0, 1.0}}};*/
	std::map<std::string, std::pair<double, double>> constraints = {};
	// cond_list.push_back(ring("c11-1/2"));
	auto ideal = Ideal(cond_list);
	fmt::print("\ngeneral form (before gröbner)\n");
	dump(ideal);
	fmt::print("\ngeneral form (after gröbner)\n");
	ideal.groebner();
	dump(ideal);

	auto double_ring = PolynomialRing<double, RANK>(ring.var_names());
	analyze_variety(ideal.change_ring(&double_ring, rationalToDouble),
	                constraints, true);
}
