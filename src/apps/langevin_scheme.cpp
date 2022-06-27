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

// terms of the transition operator
bool abelian = false; // if true, set cA and structure constant to zero
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
		// dump_summary(full_condition[2 * k]);
		for (auto term : full_condition[2 * k].terms())
		{
			bool skip = false;
			if (abelian)
				for (auto const &atom : term.index.atoms())
					if (atom.symbol == "f" || atom.symbol == "cA")
						skip = true;
			if (skip)
				continue;

			ideal.push_back(term.coefficient);
			fmt::print("{} :: {} + ... ( {} terms total )\n", term.index,
			           term.coefficient.lm(), term.coefficient.terms().size());
		}
	}
	return ideal;
}

std::vector<R> simplifyNLO(std::vector<R> const &nlo,
                           std::vector<R> const &ideal)
{
	std::vector<R> r;

	// NOTE: nlo conditions are overconstrained. Therefore, we only remove
	//       exact duplicates, not linear combinations.
	for (auto f : nlo)
	{
		reduce(f, ideal);
		if (f == 0)
			continue;
		f /= f.lc();
		for (auto const &g : r)
			if (g == f)
				goto next;
		r.push_back(f);
	next:;
	}
	return r;
}

template <typename T> void append(std::vector<T> &a, std::vector<T> const &b)
{
	a.reserve(a.size() + b.size());
	a.insert(a.end(), b.begin(), b.end());
}

int main()
{
	int order = 3;
	abelian = false;

	std::vector<R> ideal = {};
	std::vector<R> nlo = {};
	std::vector<R> errs = {};

	// second order, two noises
	if (order == 2)
	{
		// NOTE: abelian and 1D are only equivalent to second order.
		//       (so it even breaks for two-step NLO)

		// general abalian constraints:
		// f0 = k2^2 + k3^2 - k1         // eliminates k2 (only exists squared)
		// f1 = k2*k5 - 1/2*k5*k1 - 1/4  // eliminates ???
		// f2 = k4 + k5 - 1              // eliminates k4

		// abelian NLO:
		// E0 = k2 - 1/2*k1 - 1/3
		// E1 = k5^2*k1 - 3/4*k2 + 1/4*k5*k1^2 - k5*k1 + 3/8*k1 + 1/4
		// E2 = k5^2*k1 - k5*k1
		// E3 = k5^2*k1 - 1/2*k2 + 1/2*k5*k1^2 - k5*k1 + 1/4*k1 + 1/6

		// rewritten
		// E0 = k2 - 1/2*k1 - 1/3
		// E1 = E2 - 3/4*E0 + 1/4*k5*k1^2
		// E2 = k5*(k5-1)*k1
		// E3 = E2 - 1/2*E0 + 1/2*k5*k1^2

		// cleanest result (especially for NLO)
		ring.generator("k4", 100);
		ring.generator("k3", 90);
		ring.generator("k2", 80);
		ring.generator("k5", 70);
		ring.generator("k1", 10);

		// the basic building blocks
		auto eta1 = Term(Covariant<R>(Indexed("eta1", 1)));
		auto eta2 = Term(Covariant<R>(Indexed("eta2", 1)));
		auto S0 = Term(Covariant<R>(Indexed("S", 1)));

		fmt::print("---------- Terms of the scheme ----------\n");

		// NOTE: * for second order, we dont need explicit commutators,
		//         so we dont even consider them
		//       * also, we dont consider terms that contribute just to NLO
		//       * finally, take the S0/S1 terms in the same proportion as
		//         S0*cA / S1*cA. This makes the non-abelian stuff unique
		//         for all cases

		auto f1 = Term(0);
		add_term(f1, "k1", S0, 2);
		add_term(f1, "k2", eta1, 1);
		add_term(f1, "k3", eta2, 1);

		auto S1 = myTaylor(S0, -f1, 2 * max_order);

		auto f2 = Term(0);
		add_term(f2, "k4", S0, 2);
		add_term(f2, "k5", S1, 2);
		add_term(f2, "1", eta1, 1); // scale setting
		add_term(f2, "k4*k6", S0 * cA, 4);
		add_term(f2, "k5*k6", S1 * cA, 4);

		append(ideal, makeOrderConditions(f2, 2));
		append(nlo, makeOrderConditions(f2, 3));

		// make "k1 != 0" explicit. this removes some spurious imaginary-noise
		// solutions and might make groebner faster...

		ideal.push_back(ring("k1*k1inv-1"));

		errs.push_back(ring("k2 - 1/2*k1 - 1/3"));
		errs.push_back(
		    ring("k5^2*k1 - 3/4*k2 + 1/4*k5*k1^2 - k5*k1 + 3/8*k1 + 1/4"));
		errs.push_back(ring("k5^2*k1 - k5*k1"));
		errs.push_back(
		    ring("k5^2*k1 - 1/2*k2 + 1/2*k5*k1^2 - k5*k1 + 1/4*k1 + 1/6"));

		// abelian schemes:
		//   * set e2=0 -> no other vanishing possible
		//              -> other errors are of the form k1^2 + const
		//              -> minimized by k3=k4=0
		//              -> torrero (0.086 | 0 + 1)
		//              -> (also one weird real but >1 solution...)
		//   * B2a: set e0=0, e1=0 -> one solution (1.0 | 0.25 + 0.75)
		//   * B2b: set e0=0, e3=0 -> one solution (0.5 | 0.25 + 0.75)
		//   * B2c: set e1=0, e3=0 -> one solution (0.183 | 0.183 + 0.817)

		// torrero
		// ideal.push_back(errs[2]);
		// ideal.push_back(ring("k3"));

		// B2a
		// ideal.push_back(errs[0]);
		// ideal.push_back(errs[1]);

		// B2b
		// ideal.push_back(errs[0]);
		// ideal.push_back(errs[3]);

		// B2c
		// ideal.push_back(errs[1]);
		// ideal.push_back(errs[3]);

		reduce(ideal);
		groebner(ideal);
	}

	else if (order == 3)
	{
		// NOTE on abelian case
		// * single noise:
		//     - implies first-order condition for aux steps, 1D left
		// * three noises:
		//     - 5D left, (example independent set k3,5,6,7,8)
		//     - does NOT imply first-order for aux steps
		//     - dimension is weird because it has 3 more vars than single-noise

		auto eta1 = Term(Covariant<R>(Indexed("eta1", 1)));
		auto eta2 = Term(Covariant<R>(Indexed("eta2", 1)));
		auto eta3 = Term(Covariant<R>(Indexed("eta3", 1)));
		auto S0 = Term(Covariant<R>(Indexed("S", 1)));

		fmt::print("---------- Terms of the scheme ----------\n");

		auto f1 = Term(0);
		add_term(f1, "k1", S0, 2);
		add_term(f1, "k1*k12", S0 * cA, 4);
		add_term(f1, "k2", eta1, 1);
		add_term(f1, "k2*k13", eta1 * cA, 3);
		// add_term(f1, "k3", eta2, 1);
		// add_term(f1, "k4", eta3, 1);
		auto S1 = myTaylor(S0, -f1, 2 * max_order);

		auto f2 = Term(0);
		add_term(f2, "k5", S0, 2);
		add_term(f2, "k5*k15", S0 * cA, 4);
		add_term(f2, "k6", S1, 2);
		add_term(f2, "k6*k16", S1 * cA, 4);
		add_term(f2, "k7", eta1, 1);
		add_term(f2, "k7*k17", eta1 * cA, 3);
		// add_term(f2, "k8", eta2, 1);
		auto S2 = myTaylor(S0, -f2, 2 * max_order);

		// even terms
		auto f3 = Term(0);
		add_term(f3, "k9", S0, 2);
		add_term(f3, "k9*k18", S0 * cA, 4);
		add_term(f3, "k9*k19", S0 * cA * cA, 6);
		add_term(f3, "k10", S1, 2);               // torrero-like
		add_term(f3, "k10*k22", S1 * cA, 4);      // torrero-like
		add_term(f3, "k10*k23", S1 * cA * cA, 6); // torrero-like
		add_term(f3, "k11", S2, 2);
		add_term(f3, "k11*k20", S2 * cA, 4);
		add_term(f3, "k11*k21", S2 * cA * cA, 6);
		add_term(f3, "1", eta1, 1); // scale setting

		// ideal.push_back(ring("k11*k11inv-1"));

		append(ideal, makeOrderConditions(f1, 1)); // implied with single noise
		append(ideal, makeOrderConditions(f2, 1)); // implied with single noise
		append(ideal, makeOrderConditions(f3, 3));
		// append(nlo, makeOrderConditions(f2, 4));

		reduce(ideal);
		// groebner(ideal); // too slow, use Singular/Maple instead
	}

	else if (order == 4)
	{
		// NOTE on abelian case:
		// * 4 steps with 4 noises have no solution at all
		//   (the old numerival solution must have been a 1D fluke)

		auto eta1 = Term(Covariant<R>(Indexed("eta1", 1)));
		auto eta2 = Term(Covariant<R>(Indexed("eta2", 1)));
		auto eta3 = Term(Covariant<R>(Indexed("eta3", 1)));
		auto eta4 = Term(Covariant<R>(Indexed("eta3", 1)));
		auto S0 = Term(Covariant<R>(Indexed("S", 1)));

		fmt::print("---------- Terms of the scheme ----------\n");

		auto f1 = Term(0);
		add_term(f1, "k1", S0, 2);
		add_term(f1, "k2", eta1, 1);
		add_term(f1, "k3", eta2, 1);
		add_term(f1, "k4", eta3, 1);
		add_term(f1, "k5", eta4, 1);
		auto S1 = myTaylor(S0, -f1, 2 * max_order);

		auto f2 = Term(0);
		add_term(f2, "k6", S0, 2);
		add_term(f2, "k7", S1, 2);
		add_term(f2, "k8", eta1, 1);
		add_term(f2, "k9", eta2, 1);
		add_term(f2, "k10", eta3, 1);
		auto S2 = myTaylor(S0, -f2, 2 * max_order);

		auto f3 = Term(0);
		add_term(f3, "k11", S0, 2);
		add_term(f3, "k12", S1, 2);
		add_term(f3, "k13", S2, 2);
		add_term(f3, "k14", eta1, 1);
		add_term(f3, "k15", eta2, 1);
		auto S3 = myTaylor(S0, -f3, 2 * max_order);

		auto f4 = Term(0);
		add_term(f4, "k16", S0, 2);
		add_term(f4, "k17", S1, 2);
		add_term(f4, "k18", S2, 2);
		add_term(f4, "k19", S3, 2);
		add_term(f4, "1", eta1, 1);

		append(ideal, makeOrderConditions(f4, 4));

		reduce(ideal);
		groebner(ideal);
	}

	else
	{
		fmt::print("ERROR: invalid order = {}\n", order);
		return -1;
	}

	fmt::print("\n===== ideal (saved to ideal.singular) =====\n");
	dump_singular(ideal, "ideal.singular");
	dump_maple(ideal, "ideal.maple");
	dump(ideal);

	fmt::print("\n===== NLO =====\n");
	nlo = simplifyNLO(nlo, ideal);
	dump(nlo);

	fmt::print("\n===== errors =====\n");
	for (size_t i = 0; i < errs.size(); ++i)
	{
		auto e = errs[i];
		reduce(e, ideal);
		fmt::print("e{} = {}\n", i, e);
	}

	analyze_variety(change_ring<double>(ideal), {}, true);
}
