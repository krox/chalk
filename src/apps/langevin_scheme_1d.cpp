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

using Real = double;
// using Real = FloatingDoubleOctuple;

/**
 * For the 1D case, we represent the derivatives as S', S'', ..., i.e. without
 * any indexes (even though we use the 'Indexed' class).
 * taylor() -> taylor_1d()
 * diff() -> diff_1d()
 * noise term(s) 'eta' are represented as polynomial vars, instead of 'indexed'
 */

// polynomial ring of coefficients
static constexpr size_t RANK = 23; // needs space for eta1,2,3 too
using R = SparsePolynomial<Rational, RANK>;
auto ring = PolynomialRing<Rational, RANK>();

// terms of the transition operator */
static constexpr int max_order = 4;
using Term = Series<Covariant<R>, max_order * 2>;
auto seps = Term::generator(); // sqrt(epsilon)

/** diff(A*exp(-S))*exp(S) */
Term myDiff(Term const &a)
{
	// NOTE: '*' not abelian in this case due to index-naming
	return mapCoefficients(
	    [](Covariant<R> const &b) {
		    return diff_1d(b) - b * Covariant<R>(Indexed("S'"));
	    },
	    a);
}

/*
 * S(f) = S + S'f + 1/2 S''ff + ...
 *   - f may not have any open index
 *   - S can have arbitrary open indices (which simply stay that way)
 */
Term myTaylor(Term const &S, Term const &force, int degree)
{
	Term result = Term(0);
	for (int d = 0; d <= degree; ++d)
	{
		Term term = S;
		R prefactor = R(1);
		for (int i = 0; i < d; ++i)
			term = mapCoefficients([](Covariant<R> c) { return diff_1d(c); },
			                       term);
		for (int i = 0; i < d; ++i)
		{
			term = term * force;
			prefactor /= (i + 1);
		}
		result += term * Scalar(prefactor);
	}
	return result;
}

/** 1D 'wick contract' (just expectation values of normal distributions) */
Term my_wick_contract(Term const &a, std::string const &eta)
{
	// NOTE: we assume Variance==2
	return mapCoefficientsNested(
	    [&eta](R const &poly) {
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
	    },
	    a);
}

void add_term(Term &a, std::string const &name, Term term, int order)
{
	while (order-- > 0)
		term *= seps;
	// term = simplify_lie_commutative(term, "S");

	fmt::print("{}      : {}\n", name, term);

	a += term * Scalar(Scalar(ring(name)));
}

template <typename T> void append(std::vector<T> &a, std::vector<T> const &b)
{
	a.reserve(a.size() + b.size());
	a.insert(a.end(), b.begin(), b.end());
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
		for (int i = 0; i < k; ++i)
			ev *= f;
		ev = my_wick_contract(ev, "eta1");
		ev = my_wick_contract(ev, "eta2");

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

	// collect order conditions for all epsilon
	// these are 'T^(k) e^-S = 0'
	std::vector<R> cond_list;
	for (int k = 0; k <= order; ++k)
	{
		auto tmp = full_condition[2 * k];
		if (tmp == 0)
			continue;

		fmt::print("---------- T^({})e^-S ----------\n", k);
		dump_summary(tmp);
		for (auto term : tmp.terms())
			cond_list.push_back(term.coefficient);
	}
	return cond_list;
}

void analyze(std::vector<R> ideal, bool full = false)
{
	reduce(ideal);
	fmt::print("\ngeneral form\n");
	dump(ideal);
	dump_singular(ideal, "ideal.singular");

	// print dependence of conditions on coefficients, in order to find
	// redunant coefficients by hand
	/*
	fmt::print("-------- Dependence of conditions on coefficients --------\n");
	for (int i = 2; i < RANK; ++i)
	{
	    fmt::print("----- {} -----\n", ring.var_names()[i]);
	    for (size_t j = 0; j < cond_list.size(); ++j)
	    {
	        auto tmp = diff(cond_list[j], i);
	        if (!(tmp == 0))
	            fmt::print("df{}/d{} = {}\n", j, ring.var_names()[i], tmp);
	    }
	}
	*/

	// numerical root finding
	// solve_numerical(change_ring<Real>(cond_list), values, {});

	// hybrid algebraic/numerical root finding
	reduce(ideal);
	dump_singular(ideal, "ideal.singular");
	dump_maple(ideal, "ideal.maple");
	fmt::print("ideal written to 'ideal.singular' and 'ideal.maple'\n");

	if (full)
	{
		fmt::print("\ngeneral form (after gröbner)\n");
		groebner(ideal);
		dump(ideal);
		analyze_variety(change_ring<Real>(ideal), {}, true);
	}
}

int main()
{
	/** NOTE:
	 *  use something like 'ring.generator(..., INT_MAX, 20)' to increase
	 *  the weight of a variable (for term-ordering), to force the gröbner-
	 *  algorithm to solve for that variable first.
	 */

	/** the forces */
	auto eta1 = ring("eta1");
	auto eta2 = ring("eta2");
	auto S0 = Term(Covariant<R>(Indexed("S'")));

	fmt::print("---------- Terms of the scheme ----------\n");

	auto f1 = Term(0);
	add_term(f1, "a1", S0, 2);
	add_term(f1, "b1", Term(Covariant<R>(eta1)), 1);
	add_term(f1, "b2", Term(Covariant<R>(eta2)), 1);
	auto S1 = myTaylor(S0, -f1, 2 * max_order);

	auto f2 = Term(0);
	add_term(f2, "a2", S0, 2);
	add_term(f2, "a3", S1, 2);
	add_term(f2, "b3", Term(Covariant<R>(eta1)), 1);
	add_term(f2, "b4", Term(Covariant<R>(eta2)), 1);
	auto S2 = myTaylor(S0, -f2, 2 * max_order);

	auto f3 = Term(0);
	add_term(f3, "a4", S0, 2);
	add_term(f3, "a5", S1, 2);
	add_term(f3, "a6", S2, 2);
	add_term(f3, "b5", Term(Covariant<R>(eta1)), 1);
	add_term(f3, "b6", Term(Covariant<R>(eta2)), 1);
	auto S3 = myTaylor(S0, -f3, 2 * max_order);

	auto f4 = Term(0);
	add_term(f4, "a7", S0, 2);
	add_term(f4, "a8", S1, 2);
	add_term(f4, "a9", S2, 2);
	add_term(f4, "a10", S3, 2);
	add_term(f4, "1", Term(Covariant<R>(eta1)), 1);
	// add_term(f4, "b7", Term(Covariant<R>(eta2)), 1);

	// order 2 torrero
	if (false)
	{
		std::vector<R> ideal = {};
		append(ideal, makeOrderConditions(f2, 2));
		ideal.push_back(ring("a2+a3-1")); // scale setting
		ideal.push_back(ring("a2"));      // torrero
		ideal.push_back(ring("b2"));
		ideal.push_back(ring("b4"));
		analyze(ideal);
	}

	// order 3 buerger
	if (false)
	{
		// only main condition: 1D variety, does not imply torrero
		// with torrero: still 1D
		// with 1st-order consistent f1: still 1D
		// with both:  0 dimensions, which also imply 2nd order f2
		std::vector<R> ideal = {};
		append(ideal, makeOrderConditions(f1, 1));
		// append(ideal, makeOrderConditions(f2, 2));
		append(ideal, makeOrderConditions(f3, 3));
		ideal.push_back(ring("a4+a5+a6-1")); // scale setting
		ideal.push_back(ring("a5"));         // torrero condition
		ideal.push_back(ring("b2"));         // no secondary noise
		ideal.push_back(ring("b4"));         // ditto
		ideal.push_back(ring("b6"));         // ditto
		analyze(ideal);
	}

	// order 4 buerger
	if (true)
	{
		// * we always want final order 4
		// * explicitly coding a10 != 0 can help performance...
		// * single noise has no solution (thanks to Maple), so we take two
		// * first approach: order 1+2+3+4 for intermediate steps
		//     * contains a1^2*a5 term, so we set a5=0
		//     * additional a9=0 OR b6=0 -> no reasonable solution
		//      (contains a1^2=0 -> b1^2+b2^2=0 -> complex noise)
		//     * kinda think this branch wont work...
		// * second approach: torrero-style a9=0 ?
		//     * actually possible to solve (barely^^)
		//     * maple can solve it modulo prime, but not in fractions
		//       (so there probably is a at least a complex solution. but
		//        not neccessarily a real one)
		//     * lets try with more conditions. maybe?
		// * third approach: order 1+1+1+4
		// * order 0+1+1+4 + first step deterministic :) (b1=b2=0)
		//      -> implies a3=a5=a6=a8=0
		// * classic RK4: (1,2,2,1)/6 weights, two noises -> no solution
		// * 3/8th rule (1,3,3,1)/8 weights, two noises -> no solution
		// * b4=b6=0 ??
		std::vector<R> ideal = {};
		ideal.push_back(ring("a7+a8+a9+a10-1")); // scale setting

		// append(ideal, makeOrderConditions(f1, 1));
		// append(ideal, makeOrderConditions(f2, 1));
		// append(ideal, makeOrderConditions(f3, 1));
		append(ideal, makeOrderConditions(f4, 4));

		// classic RK4 in the final step
		ideal.push_back(ring("a7-1/8"));
		ideal.push_back(ring("a8-3/8"));
		ideal.push_back(ring("a9-3/8"));
		ideal.push_back(ring("a10-1/8"));

		// ideal.push_back(ring("a9"));             // torrero-style
		// ideal.push_back(ring("a5")); // this is implied for order 1+2+3+4

		// no secondary noise -> no solution
		// ideal.push_back(ring("b2"));
		// ideal.push_back(ring("b4"));
		// ideal.push_back(ring("b6"));

		// ideal.push_back(ring("b6*b6inv-1"));

		// some constraints for non-zero coeffs (purely an optimization)
		ideal.push_back(ring("a1*a1inv-1"));
		ideal.push_back(ring("a10*a10inv-1"));

		// ideal.push_back(ring("a3-22/100")); // fixing the known scheme
		// ideal.push_back(ring("a1-39/1000")); // fixing the known scheme

		analyze(ideal);
	}

	// clang-format off
	/*
	// old solution with one negative coeff :/
	std::map<std::string, Real> values;
	values["a1"] = "3.908321893316119026777794566668589568283401149509362110424178005984885432e-02";
	values["b1"] = "1.976947620276298323650439914795655112885521032670732248024095524602647502e-01";
	values["a2"] = "-4.620041923971396347070594856631942093372344970703125e-02";
	values["a3"] = "2.199999999884254941719075611902054922906873347853657130082453179867988394e-01";
	values["b2"] = "3.13006172903248066641609807270210065593728387787400151206723762972436617e-01";
	values["b3"] = "2.753665129843741665736203963900263231728657221307780834856404733628975428e-01";
	values["a4"] = "2.034816600036364722872692491525818820250224115172562465471236127468572616e-02";
	values["a5"] = "1.688988894232784232856296747934205620583317238006163223452880181317366614e-01";
	values["a6"] = "2.950852361472194394726440667876265148841128210422320774218102888531632754e-01";
	values["b4"] = "2.726958234492407830607154581780401402360575799387478118607438633623190581e-01";
	values["b5"] = "6.402884345700787454302048876824305036164568682185923885382016865235213235e-01";
	values["a7"] = "1.673009623048142407352103302552170404135662511363483595150755404410507553e-02";
	values["a8"] = "4.335605412782167732260382266356757442556845133200647524264956687702019805e-01";
	values["a9"] = "2.845766541723607439590701683025452680522902612705103145525209828795567149e-01";
	values["a10"] = "2.651327083189410587413705720362572836506686002957900970694757943061362213e-01";
	values["b6"] = "6.112713096216171893236362634950972685822938167093318175308924026549222241e-01";
	values["b7"] = "7.914211180108053021299056291610290099230055271388432783984667489968495298e-01";
	*/
	// clang-format on
}
