#ifndef CHALK_IDEAL_H
#define CHALK_IDEAL_H

/** ideals in multivariate polynomial rings */

#include "chalk/numerics.h"
#include "chalk/sparse_polynomial.h"
#include <map>

namespace chalk {

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

template <typename R, size_t rank> class Ideal
{
  public:
	using polynomial_ring_t = PolynomialRing<R, rank>;
	using polynomial_t = SparsePolynomial<R, rank>;

  private:
	std::vector<polynomial_t> basis_;

	void cleanup()
	{
		// pivotize
		bool change = false;
		do
		{
			change = false;
			for (size_t i = 0; i < basis_.size(); ++i)
			{
				if (basis_[i] == 0)
					continue;
				basis_[i] /= basis_[i].terms()[0].coefficient;
				for (size_t j = 0; j < basis_.size(); ++j)
					if (j != i)
						change |= chalk::reduce(basis_[j], basis_[i]);
			}
		} while (change == true);

		// remove zero entries
		basis_.erase(
		    std::remove_if(basis_.begin(), basis_.end(),
		                   [](polynomial_t const &poly) { return poly == 0; }),
		    basis_.end());

		// sort by leading term
		std::sort(basis_.begin(), basis_.end(), [](auto &a, auto &b) {
			return order_degrevlex(a.terms()[0], b.terms()[0],
			                       a.ring()->weights());
		});
	}

  public:
	Ideal() = default;
	explicit Ideal(std::vector<polynomial_t> polys) : basis_{polys}
	{
		cleanup();
	}

	/** reduce f using all polynomials in this basis */
	void reduce(polynomial_t &f) const
	{
		while (true)
		{
			bool change = false;
			for (auto &g : basis_)
			{
				change |= chalk::reduce(f, g);
				if (f.terms().empty())
					return;
			}
			if (!change)
				return;
		}
	}

	/** partially reduce f using all polynomials in this basis */
	void reduce_partial(polynomial_t &f) const
	{
		while (true)
		{
			bool change = false;
			for (auto &g : basis_)
			{
				change |= chalk::reduce_partial(f, g);
				if (f.terms().empty())
					return;
			}
			if (!change)
				return;
		}
	}

	std::vector<polynomial_t> const &basis() const { return basis_; }
	PolynomialRing<R, rank> const *ring() const
	{
		if (basis_.empty())
			return SparsePolynomial<R, rank>::default_ring();
		else
			return basis_[0].ring();
	}

	/** transform into groebner basis */
	void groebner(size_t max_polys = SIZE_MAX);

	/** change basis ring */
	template <typename R2, size_t rank2, typename F>
	Ideal<R2, rank2> change_ring(PolynomialRing<R2, rank2> const *ring2,
	                             F &&convert) const
	{
		std::vector<SparsePolynomial<R2, rank2>> r;
		r.reserve(basis_.size());
		for (auto &f : basis_)
			r.push_back(f.change_ring(ring2, convert));
		return Ideal<R2, rank2>(std::move(r));
	}
};

template <typename R, size_t rank> inline void dump(Ideal<R, rank> const &ideal)
{
	fmt::print("polynomial ideal with {} variables and {} equations:\n", rank,
	           ideal.basis().size());
	for (auto &poly : ideal.basis())
		fmt::print("    {}\n", poly);
}

template <typename R, size_t rank>
inline void dump_summary(Ideal<R, rank> const &ideal)
{
	fmt::print("polynomial ideal with {} variables and {} equations:\n", rank,
	           ideal.basis().size());
	for (auto &poly : ideal.basis())
		fmt::print("    {} + ... ( {} terms total )\n", poly.lm(),
		           poly.terms().size());
}

template <typename R, size_t rank>
inline void Ideal<R, rank>::groebner(size_t max_polys)
{
	int batch_size = 5;
	// really naive algorithm adding resolvents until saturation
	while (true)
	{
		int new_polys = 0;
		int old_count = (int)basis_.size();
		for (int i = old_count - 1; i >= 0; --i)
			for (int j = old_count - 1;
			     j >= 0 && basis_.size() < max_polys && new_polys < batch_size;
			     --j)
			{
				auto r = resolvent(basis_[i], basis_[j]);
				reduce_partial(r);
				if (r.terms().empty())
					continue;
				r /= r.lc();

				/*fmt::print("poly({}) = {} + ... ( {} terms total )\n",
				           basis_.size(), r.lm(), r.terms().size());*/

				new_polys++;
				basis_.push_back(r);
			}
		if (new_polys)
		{
			// dump_summary(*this);
			// fmt::print("sweep...\n");
			cleanup();
		}
		else
			break;
	}
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
	using poly_list = std::vector<SparsePolynomial<double, rank>>;
	using value_map = std::map<int, double>;
	auto &var_names = ideal.ring()->var_names();
	int solution_count = 0;
	std::vector<std::map<std::string, double>> result;

	// stack of partial solutions
	std::vector<std::pair<poly_list, value_map>> stack;
	stack.push_back({ideal.basis(), {}});
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
					result.back()[ideal.ring()->var_names()[var]] = val;
			}
		}
	}
	if (solution_count == 0 && verbose)
		fmt::print("no solutions found\n");
	return result;
}

} // namespace chalk

#endif
