#pragma once

/** ideals in multivariate polynomial rings */

#include "chalk/numerics.h"
#include "chalk/sparse_polynomial.h"
#include "fmt/ranges.h"
#include "util/random.h"
#include <Eigen/Dense>
#include <fmt/os.h>
#include <iostream>
#include <map>

namespace chalk {

/**
 * NOTE: previously, there was an actual 'Ideal' class. But after some thought,
 *       thats not too useful because we don't actually want to enforce a
 *       strong normal form all the time. Also it lead to plenty of
 *       code-duplication, because operations on raw vector<Polynomial> were
 *       still neccessary.
 */
template <typename R, size_t rank>
using Ideal = std::vector<SparsePolynomial<R, rank>>;

/** reduce f using g */
template <typename R, size_t rank>
bool reduce(SparsePolynomial<R, rank> &f, SparsePolynomial<R, rank> const &g)
{
	if (g.terms().empty())
		return false;
	assert(f.ring() == g.ring()); // need consistent term-order

	bool change = false;
	for (size_t i = 0; i < f.terms().size(); ++i)
		if (auto d = f.terms()[i] / g.terms()[0]; d)
		{
			f -= g * d.value();
			--i;
			change = true;
		}
	return change;
}

/** reduce f using g, leading term only */
template <typename R, size_t rank>
bool reduce_partial(SparsePolynomial<R, rank> &f,
                    SparsePolynomial<R, rank> const &g)
{
	if (g.terms().empty() || f.terms().empty())
		return false;
	assert(f.ring() == g.ring()); // need consistent term-order

	if (auto d = f.terms()[0] / g.terms()[0]; d)
	{
		f -= g * d.value();
		reduce_partial(f, g);
		return true;
	}
	else
		return false;
}

template <typename R, size_t rank>
bool reduce(SparsePolynomial<R, rank> &f, Ideal<R, rank> const &gs)
{
	bool r = false;
	while (true)
	{
		bool change = false;
		for (auto const &g : gs)
			change |= reduce(f, g);
		if (change)
			r = true;
		else
			break;
	}
	return r;
}

template <typename R, size_t rank>
bool reduce_partial(SparsePolynomial<R, rank> &f, Ideal<R, rank> const &gs)
{
	bool r = false;
	while (true)
	{
		bool change = false;
		for (auto const &g : gs)
			change |= reduce_partial(f, g);
		if (change)
			r = true;
		else
			break;
	}
	return r;
}

template <typename R, size_t rank> void reduce(Ideal<R, rank> &fs)
{
	if (fs.empty())
		return;
	auto ring = fs[0].ring();

	// normalize polynomials
	for (auto &f : fs)
		if (!(f == 0))
			f /= f.terms()[0].coefficient;

	// remove zero entries
	fs.erase(std::remove_if(fs.begin(), fs.end(),
	                        [](SparsePolynomial<R, rank> const &poly) {
		                        return poly == 0;
	                        }),
	         fs.end());

	// sort by leading term
	std::sort(fs.begin(), fs.end(), [ring](auto &a, auto &b) {
		return order_degrevlex(a.terms()[0], b.terms()[0], ring->weights());
	});

	bool change = false;
	for (size_t i = 0; i < fs.size(); ++i)
		for (size_t j = 0; j < i; ++j)
			change |= reduce(fs[j], fs[i]);

	if (change)
		reduce(fs);
}

template <typename R, size_t rank> void reduce_partial(Ideal<R, rank> &fs)
{
	if (fs.empty())
		return;
	auto ring = fs[0].ring();

	// normalize polynomials
	for (auto &f : fs)
		if (!(f == 0))
			f /= f.terms()[0].coefficient;

	// remove zero entries
	fs.erase(std::remove_if(fs.begin(), fs.end(),
	                        [](SparsePolynomial<R, rank> const &poly) {
		                        return poly == 0;
	                        }),
	         fs.end());

	// sort by leading term
	std::sort(fs.begin(), fs.end(), [ring](auto &a, auto &b) {
		return order_degrevlex(a.terms()[0], b.terms()[0], ring->weights());
	});

	bool change = false;
	for (size_t i = 0; i < fs.size(); ++i)
		for (size_t j = 0; j < i; ++j)
			change |= reduce_partial(fs[j], fs[i]);

	if (change)
		reduce_partial(fs);
}

template <typename R, size_t rank>
SparsePolynomial<R, rank> resolvent(SparsePolynomial<R, rank> const &a,
                                    SparsePolynomial<R, rank> const &b)
{
	assert(!a.terms().empty());
	assert(!b.terms().empty());
	std::array<int, rank> ex;
	for (size_t i = 0; i < rank; ++i)
		ex[i] = std::max(a.terms()[0].exponent[i], b.terms()[0].exponent[i]);
	auto mon = Monomial<R, rank>{R(1), ex};

	return a * (mon / a.terms()[0]).value() - b * (mon / b.terms()[0]).value();
}

/** enumerate variables that actually occur in the polynomials */
template <typename R, size_t rank>
std::vector<std::string> active_variables(Ideal<R, rank> const &fs)
{
	if (fs.empty())
		return {};

	std::array<bool, rank> occ;
	occ.fill(false);
	for (auto const &f : fs)
		for (auto &mon : f.terms())
			for (size_t i = 0; i < rank; ++i)
				if (mon.exponent[i] != 0)
					occ[i] = true;
	std::vector<std::string> r;
	for (size_t i = 0; i < rank; ++i)
		if (occ[i])
			r.push_back(fs[0].ring()->var_names()[i]);
	return r;
}

/** change basis ring */
template <typename R2, typename R, size_t rank>
Ideal<R2, rank> change_ring(Ideal<R, rank> const &fs)
{
	if (fs.empty())
		return {};
	auto convert = [](R const &a) -> R2 { return R2(a); };

	// this leaks memory. in the future all rings should be stored globally
	auto ring2 = new PolynomialRing<R2, rank>(fs[0].ring()->var_names());

	Ideal<R2, rank> r;
	r.reserve(fs.size());
	for (auto const &f : fs)
		r.push_back(f.change_ring(ring2, convert));
	return r;
}

template <typename R, size_t rank> inline void dump(Ideal<R, rank> const &ideal)
{
	auto vars = active_variables(ideal);
	fmt::print("polynomial ideal with {} equations for {} variables ({})\n",
	           ideal.size(), vars.size(), vars);
	for (size_t i = 0; i < ideal.size(); ++i)
		fmt::print("    f{} = {}\n", i, ideal[i]);
}

template <typename R, size_t rank>
inline void dump_summary(Ideal<R, rank> const &ideal)
{
	auto vars = active_variables(ideal);
	fmt::print("polynomial ideal with {} equations for {} variables ({})\n",
	           ideal.size(), vars.size(), vars);
	for (auto const &poly : ideal)
	{
		if (poly.terms().size() <= 10)
			fmt::print("    {} ( {} terms total )\n", poly,
			           poly.terms().size());
		else
			fmt::print("    {} + ... + {}   ( {} terms total )\n", poly.lm(),
			           poly.trailing_terms(2), poly.terms().size());
	}
}

/**
 * Write the ideal with the current basis polynomials into a textfile that
 * can be read by Singular.
 */
template <typename R, size_t rank>
inline void dump_singular(Ideal<R, rank> const &ideal,
                          std::string const &filename)
{
	// 1) define the ring of polynomials (using only the used variables)
	auto file = fmt::output_file(filename);
	// "QQ" is the field of coefficients, here its rational numbers.
	// Alternatively, some finite field could be used to speed up computations.
	file.print("ring r = QQ,(");
	auto vars = active_variables(ideal);
	for (size_t i = 0; i < vars.size(); ++i)
		file.print(i == 0 ? "{}" : ",{}", vars[i]);
	// "dp" is the monomial ordering , here its "degree reverse lexicographical"
	// which is the suggested default (usually much better than just
	// lexicographic). The more advanced orderings (using weights and stuff)
	// require more knowlege about the variables.
	file.print("),dp;\n");

	// 2) define the basis polynomials
	for (size_t i = 0; i < ideal.size(); ++i)
		file.print("poly f{} = {};\n", i, ideal[i]);

	// 3) define the ideal
	file.print("ideal I = ");
	for (size_t i = 0; i < ideal.size(); ++i)
		file.print(i == 0 ? "f{}" : ",f{}", i);
	file.print(";\n");

	// 4) compute groebner basis
	// file.print("ideal Is = groebner(I);\n");
	// file.print("print(dim(Is));\n");
}

/** write a file to be used by Maple to find a Gröbner basis. */
template <typename R, size_t rank>
inline void dump_maple(Ideal<R, rank> const &ideal, std::string const &filename)
{
	auto file = fmt::output_file(filename);

	for (size_t i = 0; i < ideal.size(); ++i)
		file.print("f{} := {};\n", i, ideal[i]);

	// 3) define the ideal
	file.print("ideal := [");
	for (size_t i = 0; i < ideal.size(); ++i)
		file.print(i == 0 ? "f{}" : ",f{}", i);
	file.print("];\n");

	file.print("with(Groebner);\n");
	file.print("infolevel[GroebnerBasis] := 5;\n");
	file.print("vars := [");
	auto vars = active_variables(ideal);
	for (size_t i = 0; i < vars.size(); ++i)
		file.print(i == 0 ? "{}" : ",{}", vars[i]);
	file.print("];\n");
	file.print("mybasis := Basis(ideal, tdeg);\n");
}

/** write a file to be used by PHCpack */
template <typename R, size_t rank>
inline void dump_phcpack(Ideal<R, rank> const &ideal,
                         std::string const &filename)
{
	auto file = fmt::output_file(filename);
	file.print("{}\n", ideal.size());
	for (auto const &f : ideal)
		file.print("{};\n", f);
}

template <typename R, size_t rank>
inline void groebner(Ideal<R, rank> &fs, size_t max_polys = 999999999)
{
	int batch_size = 5;
	// really naive algorithm adding resolvents until saturation
	while (true)
	{
		int new_polys = 0;
		int old_count = (int)fs.size();
		for (int i = old_count - 1; i >= 0; --i)
			for (int j = old_count - 1;
			     j > i && fs.size() < max_polys && new_polys < batch_size; --j)
			{
				auto r = resolvent(fs[i], fs[j]);
				reduce_partial(r, fs);
				if (r.terms().empty())
					continue;
				r /= r.lc();

				/*fmt::print("poly({}) = {} + ... ( {} terms total )\n",
				           basis_.size(), r.lm(), r.terms().size());*/

				new_polys++;
				fs.push_back(r);
			}
		if (new_polys)
		{
			// dump_summary(*this);
			// fmt::print("sweep...\n");
			reduce_partial(fs);
		}
		else
			break;
	}
	reduce(fs);
}

/** if a variable occurs in exactly one polynomial, remove that polynomial */
template <typename R, size_t rank>
std::vector<std::string> removeTrivialVariables(Ideal<R, rank> &ideal)
{
	std::vector<std::string> r;
	// build table of which variable occurs in which polynomial
	auto occs = std::vector<std::bitset<rank>>(ideal.size());
	for (size_t i = 0; i < ideal.size(); ++i)
		occs[i] = ideal[i].var_occs();

// step 1: look for variables which only occur in a single polynomial
again:
	for (size_t k = 0; k < rank; ++k)
	{
		size_t pivot = (size_t)-1;
		for (size_t i = 0; i < occs.size(); ++i)
			if (occs[i][k])
			{
				if (pivot == (size_t)-1)
					pivot = i;
				else
					goto next;
			}
		if (pivot != (size_t)-1)
		{
			r.push_back(ideal[pivot].ring()->var_names()[k]);
			ideal[pivot] = SparsePolynomial<R, rank>(0);
			occs[pivot].reset();
			goto again;
		}
	next:;
	}

	if (!r.empty())
		reduce_partial(ideal);
	return r;
}

/* some future project to have some more clever analysis than just gröbner */
template <typename R, size_t rank>
Ideal<R, rank> analyze_ideal(Ideal<R, rank> &ideal)
{
	fmt::print("--------- analyzing polynomial ideal ----------\n");
	assert(!ideal.empty());
	auto ring = ideal[0].ring();
	using Poly = SparsePolynomial<R, rank>;

	reduce(ideal);
	auto removedVars = removeTrivialVariables(ideal);
	fmt::print("trivial variables: {}\n", removedVars);

	auto occs = std::bitset<rank>();
	for (size_t i = 0; i < ideal.size(); ++i)
		occs |= ideal[i].var_occs();

	// look for variables that only occur together with certain others
	for (int v = 0; v < (int)rank; ++v)
	{
		if (!occs[v])
			continue;
		std::array<int, rank> m;
		for (int i = 0; i < (int)rank; ++i)
			m[i] = i == v ? 0 : INT_MAX;
		for (Poly const &f : ideal)
			for (auto const &term : f.terms())
			{
				if (term.exponent[v] == 0)
					continue;

				for (int w = 0; w < (int)rank; ++w)
					if (w != v)
						m[w] =
						    std::min(m[w], term.exponent[w] / term.exponent[v]);
			}

		fmt::print("{} appears with other vars: {}\n", ring->var_names()[v], m);
	}

	// look for polys that can be used directly as definition of a variable
	for (int v = 0; v < (int)rank; ++v)
	{
		for (Poly const &f : ideal)
			if (f.max_order(v) == 1 && get_coefficient(f, v, 1).isConstant())
			{
				auto s = f.solveFor(v);
				fmt::print("variable {} is defined by poly {} as {} = {}\n",
				           ring->var_names()[v], f, ring->var_names()[v], s);
				for (Poly &g : ideal)
					g = g.substitute(v, s);
				break;
			}
	}

	reduce(ideal);

	return {};
}

/**
 * Numerically solve as many basis polynomials as possible.
 * For best result, compute groebner basis first, but even without, it might
 * produce something useful.
 */
template <size_t rank>
inline std::vector<std::map<std::string, double>> analyze_variety(
    Ideal<double, rank> const &ideal,
    std::map<std::string, std::pair<double, double>> const &constraints = {},
    bool verbose = true)
{
	assert(!ideal.empty());
	using poly_list = std::vector<SparsePolynomial<double, rank>>;
	using value_map = std::map<int, double>;
	auto &var_names = ideal[0].ring()->var_names();
	int solution_count = 0;
	std::vector<std::map<std::string, double>> result;

	// stack of partial solutions
	std::vector<std::pair<poly_list, value_map>> stack;
	stack.push_back({ideal, {}});
	while (!stack.empty())
	{
		auto state = std::move(stack.back());
		stack.pop_back();

		// look for a univariate polynomial
		bool found_univariate = false;
		for (auto &f : state.first)
		{
			// univariate found -> solve and branch
			if (f.var_count() == 1)
			{
				auto [pivot, poly] = f.to_univariate();
				auto &var_name = var_names[pivot];
				/*fmt::print("found univariate (var={}): {}\n",
				           ideal.ring()->var_names()[pivot], poly);*/
				auto sols = roots(poly);
				for (double sol : sols)
				{
					if (constraints.count(var_name))
						if (sol < constraints.at(var_name).first ||
						    sol > constraints.at(var_name).second)
							continue;

					stack.push_back(state);
					stack.back().second[pivot] = sol;
					for (auto &g : stack.back().first)
						g = g.substitute(pivot, sol);
				}
				found_univariate = true;
				break;
			}
		}

		// no univariate found -> print solution (and remaining polys)
		if (!found_univariate)
		{
			// check for contradictions
			bool contradiction = false;
			for (auto &f : state.first)
				if (f.terms().size() == 1)
					if (f.var_count() == 0)
						if (std::abs(f.terms()[0].coefficient) > 1.0e-10)
						{
							contradiction = true;
							break;
						}
			if (contradiction)
				continue;

			// print (partial) variable assignment
			if (verbose)
			{
				fmt::print("solution {}:\n", ++solution_count);
				for (auto &[var, val] : state.second)
					fmt::print("    {} = {}\n", var_names[var], val);
			}
			bool has_conditions = false;
			// print remaining conditions
			for (auto &f : state.first)
			{
				if (f.terms().empty())
					continue;
				if (f.var_count() == 0)
					if (std::abs(f.terms()[0].coefficient) < 1.0e-14)
						continue;
				has_conditions = true;
				if (verbose)
					fmt::print("    0 = {}\n", f);
			}

			if (!has_conditions)
			{
				result.push_back({});
				for (auto &[var, val] : state.second)
					result.back()[ideal[0].ring()->var_names()[var]] = val;
			}
		}
	}
	if (solution_count == 0 && verbose)
		fmt::print("no solutions found\n");
	return result;
}

template <typename R, size_t rank>
void solve_numerical(Ideal<R, rank> const &ideal,
                     std::map<std::string, R> startingValues = {},
                     std::vector<std::string> fixedVars = {})
{
	using Poly = SparsePolynomial<R, rank>;
	using Matrix = Eigen::Matrix<R, Eigen::Dynamic, Eigen::Dynamic>;
	using Vector = Eigen::Matrix<R, Eigen::Dynamic, 1>;
	auto ring = ideal[0].ring();
	auto rng = util::xoshiro256{std::random_device()()};

	// compute derivative polynomials
	auto derivs = std::vector<std::vector<Poly>>(ideal.size());
	for (size_t i = 0; i < ideal.size(); ++i)
	{
		derivs[i].resize(rank);
		for (size_t j = 0; j < rank; ++j)
			derivs[i][j] = diff(ideal[i], j);
	}
start:
	// random starting points
	auto x = Vector(rank);
	auto xs = util::span<const R>(x.data(), (int)x.size());
	for (size_t i = 0; i < rank; ++i)
	{
		if (startingValues.count(ring->var_names()[i]))
			x[i] = startingValues[ring->var_names()[i]];
		else
			x[i] = R(std::uniform_real_distribution<double>()(rng));
	}

	// (squared) error
	double err = 0;
	for (auto const &f : ideal)
	{
		auto e = (double)f.substitute(xs);
		err += e * e;
	}
	fmt::print("error^2 = {}\n", err);
	double last_err = err;

	for (int iter = 0;; ++iter)
	{
		// compute (pseudo-)inverse and update x
		auto f = Vector(ideal.size());
		auto fd = Matrix(ideal.size(), rank);
		for (size_t i = 0; i < ideal.size(); ++i)
		{
			f(i) = ideal[i].substitute(xs);
			for (size_t j = 0; j < rank; ++j)
				fd(i, j) = derivs[i][j].substitute(xs);
		}
		for (auto const &var : fixedVars)
			fd.col(ring->var_id(var)).setZero();
		for (size_t i = 0; i < ideal.size(); ++i)
			fmt::print("f[{}] = {}\n", i, (double)f(i));
		auto tmp = fd.completeOrthogonalDecomposition();
		std::cout << fd << std::endl;
		assert(tmp.info() == Eigen::Success);
		fmt::print("matrix is {} x {}, rank = {}\n", fd.rows(), fd.cols(),
		           tmp.rank());
		Vector step = tmp.solve(f);
		auto scale = step.template lpNorm<Eigen::Infinity>();
		for (size_t i = 0; i < rank; ++i)
			fmt::print("step[{}] = {}\n", i, (double)step(i));
		fmt::print("steps size = {}\n", (double)scale);
		if (scale > 0.05)
			x -= (0.05 / scale) * step;
		else
			x -= step;

		// (squared) error
		err = 0;
		for (auto const &f : ideal)
		{
			auto e = (double)f.substitute(xs);
			err += e * e;
		}

		if (iter % 50 == 0)
		{
			fmt::print("iter = {}, error^2 = {}\n", iter, err);
			for (size_t i = 0; i < rank; ++i)
				fmt::print("values[\"{}\"] = \"{}\";\n",
				           ideal[0].ring()->var_names()[i], x(i));
			if (x.template lpNorm<Eigen::Infinity>() > 10.0)
				goto start;
		}
		if (last_err < 1e-100)
			break;
		last_err = err;
	}
	fmt::print("--------- end result ----------\n");
	for (size_t i = 0; i < rank; ++i)
		fmt::print("values[\"{}\"] = \"{}\";\n",
		           ideal[0].ring()->var_names()[i], x(i));
}

} // namespace chalk
