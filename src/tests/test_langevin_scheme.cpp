
#include "chalk/bch.h"
#include "chalk/covariant.h"
#include "chalk/floating.h"
#include "chalk/fraction.h"
#include "chalk/ideal.h"
#include "chalk/series.h"
#include "chalk/sparse_polynomial.h"
#include <cassert>
#include <fmt/format.h>
#include <gtest/gtest.h>
using namespace chalk;

namespace {

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
	a += term * Scalar(Scalar(ring(name)));
}

std::vector<R> makeOrderConditions(Term const &f, int order = max_order)
{
	assert(order <= max_order);

	// full condition is
	//  (∇ <f> + 1/2 ∇∇<ff> + ...) e^-S = (T-1)e^-S
	auto full_condition = Term(0);
	for (int k = 1; k <= order * 2; ++k)
	{
		// build ev = <f^k>
		auto ev = Term(1);

		for (int i = 0; i < k; ++i)
		{
			ev *= f;
			checkRank(ev);
		}
		checkRank(ev, k);
		ev = ev.mapCoefficients(CHALK_LIFT(wick_contract), "eta", R(2));

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
		for (auto term : full_condition[2 * k].terms())
			ideal.push_back(term.coefficient);

	return ideal;
}

template <typename T> void append(std::vector<T> &a, std::vector<T> const &b)
{
	a.reserve(a.size() + b.size());
	a.insert(a.end(), b.begin(), b.end());
}
} // namespace

TEST(Langevin, Torrero)
{
	std::vector<R> ideal = {};

	/** second order schemes */
	ring.generator("k5", 100);

	/** the forces */
	auto eta = Term(Covariant<R>(Indexed("eta", 1)));
	auto S0 = Term(Covariant<R>(Indexed("S", 1)));

	auto f1 = Term(0);
	add_term(f1, "k1", S0, 2);
	add_term(f1, "k2", eta, 1);
	auto S1 = myTaylor(S0, -f1, 2 * max_order);

	auto f2 = Term(0);
	add_term(f2, "k3", S0, 2);
	add_term(f2, "k4", S1, 2);
	add_term(f2, "k5", S0 * cA, 4);
	add_term(f2, "1", eta, 1);

	append(ideal, makeOrderConditions(f2, 2));

	std::map<std::string, std::pair<double, double>> constraints = {
	    {"k1", {0.0, 1.0}}};

	{
		reduce(ideal);
		groebner(ideal);
		// dump(ideal);
		auto result =
		    analyze_variety(change_ring<double>(ideal), constraints, false);
		ASSERT_EQ(result.size(), 0); // no unique solution (variety dimension=2)
	}

	{
		// BF scheme
		auto ideal2 = ideal;
		ideal2.push_back(ring("k2-1"));
		groebner(ideal2);
		auto result =
		    analyze_variety(change_ring<double>(ideal2), constraints, false);

		ASSERT_EQ(result.size(), 1);
		EXPECT_NEAR(result[0].at("k1"), 1.0, 1.0e-10);
		EXPECT_NEAR(result[0].at("k2"), 1.0, 1.0e-10);
		EXPECT_NEAR(result[0].at("k3"), 0.5, 1.0e-10);
		EXPECT_NEAR(result[0].at("k4"), 0.5, 1.0e-10);
		EXPECT_NEAR(result[0].at("k5"), 0.16666666666666666, 1.0e-10);
	}

	{
		// Torrero / BBPT scheme
		auto ideal2 = ideal;
		ideal2.push_back(ring("k3"));
		groebner(ideal2);
		auto result =
		    analyze_variety(change_ring<double>(ideal2), constraints, false);

		ASSERT_EQ(result.size(), 1);
		EXPECT_NEAR(result[0].at("k1"), 0.08578643762690494, 1.0e-10);
		EXPECT_NEAR(result[0].at("k2"), 0.2928932188134524, 1.0e-10);
		EXPECT_NEAR(result[0].at("k3"), 0.0, 1.0e-10);
		EXPECT_NEAR(result[0].at("k4"), 1.0, 1.0e-10);
		EXPECT_NEAR(result[0].at("k5"), 0.0631132760733929, 1.0e-10);
	}
}
