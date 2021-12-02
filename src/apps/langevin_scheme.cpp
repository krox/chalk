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
static constexpr size_t RANK = 36;
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

/** check that all terms have the same number of open indices */
void checkRank(Term const &a, int rank = -1)
{
	for (auto const &tmp : a.coefficients())
		for (auto &[_, indexed] : tmp.terms())
		{
			assert(rank == -1 || indexed.rank() == rank);
			rank = indexed.rank();
		}
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
/** commutator [a,[b,[c,d]]] */
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

std::vector<R> makeOrderConditions(Term const &f, int order = max_order)
{
	assert(order <= max_order);

	// full condition is
	//  (∇ <f> + 1/2 ∇∇<ff> + ...) e^-S = (T-1)e^-S
	fmt::print("---------- Expectation values ----------\n");
	auto full_condition = Term(0);
	for (int k = 1; k <= order * 2; ++k)
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
		ev = ev.mapCoefficients(CHALK_LIFT(wick_contract), "eta1", R(2));
		ev = ev.mapCoefficients(CHALK_LIFT(wick_contract), "eta2", R(2));
		ev = ev.mapCoefficients(CHALK_LIFT(wick_contract), "eta3", R(2));
		ev = ev.mapCoefficients(CHALK_LIFT(wick_contract), "eta4", R(2));

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
	std::vector<R> ideal;
	for (int k = 0; k <= order; ++k)
	{
		fmt::print("---------- T^({})e^-S ----------\n", k);
		dump_summary(full_condition[2 * k]);
		for (auto term : full_condition[2 * k].terms())
			ideal.push_back(term.coefficient);
	}
	return ideal;
}

template <typename T> void append(std::vector<T> &a, std::vector<T> const &b)
{
	a.reserve(a.size() + b.size());
	a.insert(a.end(), b.begin(), b.end());
}

int main()
{
	std::vector<R> ideal = {};
	std::vector<R> nlo = {};

	/** second order schemes */
	if (true)
	{
		// ring.generator("", 100); // cA term first

		/** the forces */
		auto eta1 = Term(Covariant<R>(Indexed("eta1", 1)));
		auto eta2 = Term(Covariant<R>(Indexed("eta2", 1)));
		auto S0 = Term(Covariant<R>(Indexed("S", 1)));

		fmt::print("---------- Terms of the scheme ----------\n");

		auto f1 = Term(0);
		add_term(f1, "k1", S0, 2);
		add_term(f1, "k2", eta1, 1);
		add_term(f1, "k3", eta2, 1);
		add_term(f1, "k6", comm(eta1, S0), 3);
		add_term(f1, "k7", comm(eta2, S0), 3);
		add_term(f1, "k8", eta1 * cA, 3);
		add_term(f1, "k9", eta2 * cA, 3);
		add_term(f1, "k10", S0 * cA, 4);

		auto S1 = myTaylor(S0, -f1, 2 * max_order);

		auto f2 = Term(0);
		add_term(f2, "k4", S0, 2);
		add_term(f2, "k5", S1, 2);
		add_term(f2, "1", eta1, 1); // scale setting
		add_term(f2, "k11", comm(S0, S1), 4);
		add_term(f2, "k12", S0 * cA, 4);
		add_term(f2, "k13", S1 * cA, 4);
		auto f = f2;

		append(ideal, makeOrderConditions(f, 2));
		append(nlo, makeOrderConditions(f, 3));

		// ideal.push_back(ring("k3")); // no secondary noise
		// ideal.push_back(ring("k4")); // Torrero condition
		// ideal.push_back(ring("k1-1")); // BF condition

		ideal.push_back(ring("k1-1"));       // buerger
		ideal.push_back(ring("k2-5/6"));     // buerger
		ideal.push_back(ring("k3^2-11/36")); // buerger
		ideal.push_back(ring("k4-1/4"));     // buerger
		ideal.push_back(ring("k5-3/4"));     // buerger
	}

	/** third order schemes */
	if (false)
	{
		/** the forces */
		auto eta = Term(Covariant<R>(Indexed("eta", 1)));
		auto S0 = Term(Covariant<R>(Indexed("S", 1)));

		fmt::print("---------- Terms of the scheme ----------\n");

		auto f1 = Term(0);
		add_term(f1, "k1", S0, 2);
		// add_term(f1, "k2", comm(eta, eta, S0), 4);
		add_term(f1, "k3", eta, 1);
		// add_term(f1, "k4", comm(eta, S0), 3);
		auto S1 = myTaylor(S0, -f1, 2 * max_order);

		auto f2 = Term(0);
		add_term(f2, "k5", S0, 2); // primitive!
		add_term(f2, "k6", S1, 2);
		add_term(f2, "k7", S0 * cA, 4); // [eta,eta,S0]
		add_term(f2, "k8", eta, 1);
		// add_term(f2, "k9", comm(eta, S0), 3);
		// add_term(f2, "k10", comm(eta, S1), 3);
		auto S2 = myTaylor(S0, -f2, 2 * max_order);

		// even terms
		auto f3 = Term(0);
		add_term(f3, "k11", S0, 2);
		// add_term(f3, "k12", S1, 2); // torrero
		add_term(f3, "k13", S2, 2);
		// add_term(f3, "k14", comm(S0, S1), 4);

		// add_term(f3, "k16", comm(S1, S2), 4);
		add_term(f3, "k17", cA * S0, 4); // [eta,eta,S0]
		// add_term(f3, "k18", cA * S1, 4);        // [eta,eta,S1]
		add_term(f3, "k19", S2 * cA, 4);        // [eta,eta,S2]
		add_term(f3, "k20", S0 * (cA * cA), 6); // [eta,eta,eta,eta,S0]
		add_term(f3, "1", eta, 1);              // scale setting

		// add_term(f3, "k22", comm(eta, S1), 3);

		// add_term(f3, "k24", comm(eta, eta, eta, S0), 5);
		// add_term(f3, "k25", comm(eta, eta, eta, S1), 5);
		// add_term(f3, "k26", comm(eta, eta, eta, S2), 5);
		// add_term(f3, "k27", comm(S0, eta, S1), 5);
		// add_term(f3, "k28", comm(S0, eta, S2), 5);
		// add_term(f3, "k29", comm(S1, eta, S2), 5);

		// add_term(f3, "k31", comm(S1, eta, S1), 5);
		// add_term(f3, "k32", comm(S2, eta, S2), 5);
		// add_term(f3, "k33", comm(eta, S0, S1), 5);
		// add_term(f3, "k34", comm(eta, S0, S2), 5);
		// add_term(f3, "k35", comm(eta, S1, S2), 5);

		// add_term(f3, "k15", comm(S0, S2), 4);
		// add_term(f3, "k21", comm(S0, eta), 3);
		// add_term(f3, "k23", comm(eta, S2), 3);
		add_term(f3, "k30", comm(S0, S0, eta), 5);

		// MORE MISSING eps^5/2 TERMS. (probabaly some duplicates/linear
		// combinations)

		// discussion:
		// we impose the commutative solution (torrero + local consistency)
		// remove k14,16,18,22,25,27,29,31,33,35 as 'generalized torrero'
		// - k7,17,18,19,20 are proportional to cA/cA^2
		// - k24 does not contribute (this is general actually)
		// - k2,4 do not contribute
		// - k28,30 are the same
		// now still 7 dims left
		// - remove commutators from intermediate steps (9,10,32,34)
		// now 3 dims left
		// - 4 terms with commutators left (15,21,23,30), three can be removed
		//    * keeping k23 leads to negative and complex coeffs
		//    * the other three seem to work and even coincide in other coeffs
		//    * keeping k30 is the scheme we chose in the past

		// old notes:
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
		//   - k24 does not contribute at all, and k25/26 are equivalent, so one
		//   is
		//     enough. Furthermore, it only contributes as
		//     'eps^3*cA^2*{k3,k8}*S0' term in <f> (and not in <f^2> at all). So
		//     it can be removed completely by rescaling rescaling epsilon as s'
		//     = s + # cA^2 s^5.
		//   - k31/32 are equivalent (appear with k3/k8 respectively)

		append(ideal, makeOrderConditions(f1, 1));
		append(ideal, makeOrderConditions(f2, 2));
		append(ideal, makeOrderConditions(f3, 3));
	}

	fmt::print("\n===== ideal (saved to ideal.singular) =====\n");
	reduce(ideal);
	groebner(ideal);
	dump(ideal);
	dump_singular(ideal, "ideal.singular");

	fmt::print("\n===== NLO =====\n");
	for (auto &f : nlo)
	{
		reduce(f, ideal);
		if (!(f == 0))
			f /= f.lc();
	}
	dump(nlo);

	/*fmt::print(
	    "---------- Dependence of conditions on coefficients ----------\n");
	for (int i = 0; i < (int)RANK; ++i)
	{
	    fmt::print("----- {} -----\n", ring.var_names()[i]);
	    for (size_t j = 0; j < ideal.size(); ++j)
	    {
	        auto tmp = diff(ideal[j], i);
	        if (!(tmp == 0))
	            fmt::print("df{}/d{} = {}\n", j, ring.var_names()[i], tmp);
	    }
	}*/

	// analyze_variety(change_ring<double>(ideal), {}, true);
}
