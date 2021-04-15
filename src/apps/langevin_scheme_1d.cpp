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

// using Real = double;
using Real = FloatingDoubleOctuple;

/**
 * For the 1D case, we represent the derivatives as S', S'', ..., i.e. without
 * any indexes (even though we use the 'Indexed' class).
 * taylor() -> taylor_1d()
 * diff() -> diff_1d()
 * noise term(s) 'eta' are represented as polynomial vars, instead of 'indexed'
 */

static constexpr int max_order = 4; // orders of epsilon we want everything in

// polynomial ring of coefficients
#define RANK 20
using R = SparsePolynomial<Rational, RANK>;
auto ring = PolynomialRing<Rational, RANK>();

// R seps = ring.generator("s", max_order * 2);
// R eps = seps * seps;
using Term = Series<Covariant<R>, max_order * 2>;
auto seps = Term::generator();

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

Real toReal(Rational const &x) { return Real(x.num()) / Real(x.denom()); }

void add_term(Term &a, std::string const &name, Term term, int order)
{
	while (order-- > 0)
		term *= seps;
	// term = simplify_lie_commutative(term, "S");

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
	auto eta1 = ring("eta1");
	auto eta2 = ring("eta2");
	auto eta3 = ring("eta3");
	auto S0 = Term(Covariant<R>(Indexed("S'")));

	fmt::print("---------- Terms of the scheme ----------\n");

	auto f1 = Term(0);
	add_term(f1, "a1", S0, 2);
	add_term(f1, "b1", Term(Covariant<R>(eta1)), 1);
	auto S1 = myTaylor(S0, -f1, 2 * max_order);

	auto f2 = Term(0);
	add_term(f2, "a2", S0, 2);
	add_term(f2, "a3", S1, 2);
	add_term(f2, "b2", Term(Covariant<R>(eta1)), 1);
	add_term(f2, "b3", Term(Covariant<R>(eta2)), 1);
	auto S2 = myTaylor(S0, -f2, 2 * max_order);

	auto f3 = Term(0);
	add_term(f3, "a4", S0, 2);
	add_term(f3, "a5", S1, 2);
	add_term(f3, "a6", S2, 2);
	add_term(f3, "b4", Term(Covariant<R>(eta1)), 1);
	add_term(f3, "b5", Term(Covariant<R>(eta2)), 1);
	auto S3 = myTaylor(S0, -f3, 2 * max_order);

	auto f4 = Term(0);
	add_term(f4, "a7", S0, 2);
	add_term(f4, "a8", S1, 2);
	add_term(f4, "a9", S2, 2);
	add_term(f4, "a10", S3, 2);
	add_term(f4, "b6", Term(Covariant<R>(eta1)), 1);
	add_term(f4, "b7", Term(Covariant<R>(eta2)), 1);
	auto f = f4;

	// full condition is
	//  (∇ <f> + 1/2 ∇∇<ff> + ...) e^-S = (T-1)e^-S
	fmt::print("---------- Expectation values ----------\n");
	auto full_condition = Term(0);
	for (int k = 1; k <= max_order * 2; ++k)
	{
		// build ev = <f^k>
		auto ev = Term(1);
		for (int i = 0; i < k; ++i)
			ev *= f;
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
	cond_list.push_back(ring("a7+a8+a9+a10-1")); // overall scale setting !!!
	for (int k = 0; k <= max_order; ++k)
	{
		auto tmp = full_condition[2 * k];
		if (tmp == 0)
			continue;

		fmt::print("---------- T^({})e^-S ----------\n", k);
		dump(tmp);
		for (auto term : tmp.terms())
			cond_list.push_back(term.coefficient);
	}

	/*std::map<std::string, std::pair<double, double>> constraints = {
	    {"k1", {0.0, 1.0}}};*/
	// std::map<std::string, std::pair<double, double>> constraints = {};
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
	/* old solution
	std::map<std::string, Real> values = {
	    {"a1", 0.03908321893732925},
	    {"b1", 0.1976947620381715},
	    {"a2", -0.046200419239713963},
	    {"b2", 0.31300617291850696},
	    {"a3", 0.22},
	    {"b3", 0.2753665129875376},
	    {"a4", 0.020348166006750767},
	    {"b4", 0.2726958234381351},
	    {"a5", 0.16889888941280423},
	    {"b5", 0.6402884345777209},
	    {"a6", 0.2950852361550359},
	    {"a7", 0.016730096231808433},
	    {"b6", 0.6112713096065316},
	    {"a8", 0.433560541281889},
	    {"b7", 0.791421118022457},
	    {"a9", 0.28457665417547295},
	    {"a10", 0.2651327083108296}};*/
	std::map<std::string, Real> values;
	// clang-format off
	/* old solution with one negative coeff :/
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
	// funny new scheme with negative a4, and a9 seems to be zero?
	values["a1"] = "0.079734";
	values["b1"] = "0.282363";
	values["a2"] = "0.716393";
	values["a3"] = "0.083155";
	values["b2"] = "0.649631";
	values["b3"] = "-0.702100";
	values["a4"] = "-0.131367";
	values["a5"] = "0.444071";
	values["a6"] = "0.008058";
	values["b4"] = "0.384715";
	values["b5"] = "0.415647";
	values["a7"] = "0.065171";
	values["a8"] = "0.484072";
	//values["a9"] = "0.000007";
	values["a9"] = 0.0;
	values["a10"] = "0.450749";
	values["b6"] = "0.634689";
	values["b7"] = "0.772768";

	// clang-format on

	fmt::print("\ngeneral form (after reduce)\n");
	dump(cond_list);

	auto real_ring = PolynomialRing<Real, RANK>(ring.var_names());
	solve_numerical(change_ring(cond_list, &real_ring, toReal), values, {});

	/*for (auto const &cond : cond_list)
	    fmt::print("{}\n", cond.change_ring(&double_ring, rationalToDouble)
	                           .substitute(values));*/

	/*
	fmt::print("\ngeneral form (after gröbner)\n");
	groebner(cond_list);
	dump(cond_list);
	analyze_variety(change_ring(cond_list, &double_ring, rationalToDouble),
	                constraints, true);
	*/
}
