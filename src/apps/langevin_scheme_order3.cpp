#include "chalk/bch.h"
#include "chalk/covariant.h"
#include "chalk/floating.h"
#include "chalk/fraction.h"
#include "chalk/ideal.h"
#include "chalk/series.h"
#include "chalk/sparse_polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

// polynomial ring of coefficients
static constexpr size_t RANK = 20;
using R = SparsePolynomial<Rational, RANK>;
auto ring = PolynomialRing<Rational, RANK>();

// terms of the transition operator */
static constexpr int max_order = 3;
using Term = Series<Covariant<R>, max_order * 2>;
auto seps = Term::generator(); // sqrt(epsilon)
auto eps = seps * seps;
auto cA = Term(Covariant<R>(Indexed("cA")));

/** diff(A*exp(-S))*exp(S) */
Term myDiff(Term const &a)
{
	// NOTE: '*' not abelian in this case due to index-naming
	return mapCoefficients(
	    [](Covariant<R> const &b) {
		    return diff(b) - b * Covariant<R>(Indexed("S", 1));
	    },
	    a);
}
auto rationalToDouble(Rational const &x)
{
	return (double)x.num() / (double)x.denom();
};

auto rationalToOctuple(Rational const &x)
{
	return FloatingOctuple(x.num()) / FloatingOctuple(x.denom());
}

void checkRank(Term const &a)
{
	int rank = -1;
	for (auto const &tmp : a.coefficients())
		for (auto &[_, indexed] : tmp.terms())
		{
			assert(rank == -1 || indexed.rank() == rank);
			rank = indexed.rank();
		}
}

void checkRank(Term const &a, int rank)
{
	for (auto const &tmp : a.coefficients())
		for (auto &[_, indexed] : tmp.terms())
			assert(indexed.rank() == rank);
}

/** commutator c=[a,b]. all three interpreted as  'iT^a f_a' */
Term comm(Term const &a, Term const &b)
{
	// need exactly one open index
	checkRank(a, 1);
	checkRank(b, 1);
	auto c = -Term(Covariant<R>(Indexed("f", 0, 1, 2))) * a * b; // order matter
	// fmt::print("---------\n{}\n---------\n", c);
	checkRank(c, 5);
	c = mapCoefficients([](Covariant<R> const &x) { return contract(x, 1, 3); },
	                    c);
	c = mapCoefficients([](Covariant<R> const &x) { return contract(x, 1, 2); },
	                    c);
	checkRank(c, 1);
	return c;
}
/** commutator [a,[b,c]] */
Term comm(Term const &a, Term const &b, Term const &c)
{
	return comm(a, comm(b, c));
}
Term comm(Term const &a, Term const &b, Term const &c, Term const &d)
{
	return comm(a, comm(b, comm(c, d)));
}

/*
 * S(f) = S + S'f + 1/2 S''ff + ...
 *   - f should have exactly one open index, contracted with derivatives
 *   - S can have arbitrary open indices (which simply stay that way)
 */
Term myTaylor(Term const &S, Term const &force, int degree)
{
	checkRank(S, 1);
	checkRank(force, 1);
	Term result = Term(0);
	for (int d = 0; d <= degree; ++d)
	{
		Term term = S;
		R prefactor = R(1);
		for (int i = 0; i < d; ++i)
			term = mapCoefficients(
			    [](Covariant<R> const &c) { return diff(c); }, term);
		checkRank(term, 1 + d);
		for (int i = 0; i < d; ++i)
		{
			term = term * force;
			term = term.mapCoefficients(CHALK_LIFT(contract), d - i, d + 1 - i);
			prefactor /= (i + 1);
		}
		checkRank(term, 1);
		result += term * Scalar(prefactor);
	}
	return result;
}

void add_term(Term &a, std::string const &name, Term term, int order)
{
	while (order-- > 0)
		term *= seps;
	term = mapCoefficients([](auto const &x) { return simplify_lie(x, "S"); },
	                       term);
	fmt::print("{}      : {}\n", name, term);
	a += term * Scalar(Scalar(ring(name)));
}

int main()
{
	/** NOTE:
	 *  use something like 'ring.generator(..., INT_MAX, 20)' to increase
	 *  the weight of a variable (for term-ordering), to force the gröbner-
	 *  algorithm to solve for that variable first.
	 */

	/** the forces */
	auto eta = Term(Covariant<R>(Indexed("eta", 1)));
	auto S0 = Term(Covariant<R>(Indexed("S", 1)));

	fmt::print("---------- Terms of the scheme ----------\n");
	// first step
	auto f1 = Term(0);
	add_term(f1, "k1", S0, 2);
	add_term(f1, "k2", comm(eta, eta, S0),
	         4); // killed when removing k12 (Torrero term)
	// odd terms
	add_term(f1, "k3", eta, 1);
	// add_term(f1, "k4", comm(eta, S0), 3, 4);
	auto S1 = myTaylor(S0, -f1, 2 * max_order);

	// even terms
	auto f2 = Term(0);
	add_term(f2, "k5", S0, 2); // primitive!
	add_term(f2, "k6", S1, 2);
	// add_term(f2, "k7", comm(eta, eta, S0), 4, 4); // primitive! cA term?
	add_term(f2, "k7", S0 * cA, 4); // primitive! cA term?
	// odd terms
	add_term(f2, "k8", eta, 1);
	// add_term(f2, "k9", comm(eta, S0), 3, 4);
	// add_term(f2, "k10", comm(eta, S1), 3, 4);
	auto S2 = myTaylor(S0, -f2, 2 * max_order);

	// even terms
	auto f3 = Term(0);
	// add_term(f3, "k14", comm(S0, S1), 4);
	// add_term(f3, "k15", comm(S0, S2), 4);
	// add_term(f3, "k16", comm(S1, S2), 4);
	// add_term(f3, "k17", comm(eta, eta, S0), 4);
	// add_term(f3, "k18", comm(eta, eta, S1), 4);
	// add_term(f3, "k19", comm(eta, eta, S2), 4);
	// add_term(f3, "k24", comm(eta, eta, eta, S0), 5);
	// add_term(f3, "k25", comm(eta, eta, eta, S1), 5);
	// add_term(f3, "k26", comm(eta, eta, eta, S2), 5);
	// add_term(f3, "k21", comm(eta, S0), 3); // not needed
	// add_term(f3, "k22", comm(eta, S1), 3); // not needed

	add_term(f3, "k11", S0, 2); // primitive!
	// add_term(f3, "k12", S1, 2); // remove this for "Torrero-style"
	add_term(f3, "k13", S2, 2); // primitive!

	// add_term(f3, "k17", cA * S0, 4); // vanishes on its own
	// add_term(f3, "k18", cA * S1, 4);
	add_term(f3, "k19", S2 * cA, 4);

	add_term(f3, "k20", S0 * (cA * cA), 6); // cA^2 term

	add_term(f3, "1", eta, 1); // scale setting

	// add_term(f3, "k23", comm(eta, S2), 3); // not needed
	add_term(f3, "k30", comm(S0, eta, S0), 5); // primite!
	// add_term(f3, "k27", comm(S0, eta, S1), 5);
	// add_term(f3, "k28", comm(S0, eta, S2), 5);
	// add_term(f3, "k31", comm(S1, eta, S1), 5);
	// add_term(f3, "k29", comm(S1, eta, S2), 5);
	// add_term(f3, "k32", comm(S2, eta, S2), 5);

	// add_term(f3, "k33", comm(eta, S0, S1), 5);
	// add_term(f3, "k34", comm(eta, S0, S2), 5);
	// add_term(f3, "k35", comm(eta, S1, S2), 5);
	// MORE MISSING eps^5/2 TERMS. (probabaly some duplicates/linear
	// combinations)
	auto f = f3;

	// NOTE:
	//   - k27 = k28 = k30 (exactly)
	//   - k29 = k31 = k32 (k31/29 with 'k3', and k32 with 'k8' factor)
	//   - k33 = k34 = k35 (with k3,k8,(-k3-k8) respectively)
	//   - k30 = k31 = k33 (with '1','-k3+1','k3' factors respectively)
	//   - k4 = k9 = k14 (with k12+2*k22, k13+2*k23, -1 respectively)
	//   - k14 = k30 (with k3, 2 as factors)
	//   - In the commutative case, the torrero-style exists uniquely,
	//     so we set here k12=0 as well
	//   - k10 contributes as combination of k7 an k30, so remove it
	//   - k20 is the 'cA^2' term we expect from dimensional analysis.
	//   - k14/15/16 are essentially equivalent (appear with k3, k8, k3+k8
	//     respectively). So only one is enough.
	//   - k24 does not contribute at all, and k25/26 are equivalent, so one is
	//     enough. Furthermore, it only contributes as 'eps^3*cA^2*{k3,k8}*S0'
	//     term in <f> (and not in <f^2> at all). So it can be removed
	//     completely by rescaling rescaling epsilon as s' = s + # cA^2 s^5.
	//   - k31/32 are equivalent (appear with k3/k8 respectively)

	// full condition is
	//  (∇ <f> + 1/2 ∇∇<ff> + ...) e^-S = (T-1)e^-S
	fmt::print("---------- Expectation values ----------\n");
	auto full_condition = Term(0);
	for (int k = 1; k <= max_order * 2; ++k)
	{
		// build ev = <f^k>
		auto ev = Term(1);

		fmt::print("f = {}\n", f);
		for (int i = 0; i < k; ++i)
		{
			ev *= f;
			checkRank(ev);
		}
		checkRank(ev, k);
		ev = ev.mapCoefficients(CHALK_LIFT(wick_contract), "eta", R(2));

		fmt::print("----- <f^{}> -----\n", k);
		fmt::print("{}\n", ev);

		// build ev = 1/k! ∇^k <f^k> and add it to the complete condition
		for (int i = 1; i <= k; ++i)
		{
			ev = myDiff(ev);
			checkRank(ev, k - i + 2);
			ev = ev.mapCoefficients(CHALK_LIFT(contract), 0, k + 1 - i);
			ev = ev / i;
		}
		full_condition += ev;
	}

	full_condition =
	    full_condition.mapCoefficients(CHALK_LIFT(simplify_lie), "S");

	// collect order conditions for all epsilon
	// these are 'T^(k) e^-S = 0'
	std::vector<R> cond_list;
	for (int k = 0; k <= max_order; ++k)
	{
		fmt::print("---------- T^({})e^-S ----------\n", k);
		dump_summary(full_condition[2 * k]);
		for (auto term : full_condition[2 * k].terms())
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
	fmt::print(
	    "---------- Dependence of conditions on coefficients ----------\n");
	for (int i = 2; i < (int)RANK; ++i)
	{
		fmt::print("----- {} -----\n", ring.var_names()[i]);
		for (size_t j = 0; j < cond_list.size(); ++j)
		{
			auto tmp = diff(cond_list[j], i);
			if (!(tmp == 0))
				fmt::print("df{}/d{} = {}\n", j, ring.var_names()[i], tmp);
		}
	}

	fmt::print("\ngeneral form (after reduce)\n");
	dump(cond_list);

	auto oct_ring = PolynomialRing<FloatingOctuple, RANK>(ring.var_names());
	solve_numerical(change_ring(cond_list, &oct_ring, rationalToOctuple));

	fmt::print("\ngeneral form (after gröbner)\n");
	groebner(cond_list);
	dump(cond_list);

	auto double_ring = PolynomialRing<double, RANK>(ring.var_names());
	analyze_variety(change_ring(cond_list, &double_ring, rationalToDouble),
	                constraints, true);
}
