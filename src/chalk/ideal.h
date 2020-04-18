#ifndef CHALK_IDEAL_H
#define CHALK_IDEAL_H

/** ideals in multivariate polynomial rings */

#include "chalk/sparse_polynomial.h"

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
		if (auto d = f.terms()[i] / g.terms()[0])
		{
			f -= g * d.value();
			--i;
			change = true;
		}
	return change;
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
						change |= reduce(basis_[j], basis_[i]);
			}
		} while (change == true);

		// remove zero entries
		basis_.erase(
		    std::remove_if(basis_.begin(), basis_.end(),
		                   [](polynomial_t const &poly) { return poly == 0; }),
		    basis_.end());

		// sort by leading term
		std::sort(basis_.begin(), basis_.end(), [](auto &a, auto &b) {
			return order_degrevlex(b.terms()[0], a.terms()[0]);
		});
	}

  public:
	Ideal() = default;
	explicit Ideal(std::vector<polynomial_t> polys) : basis_{polys}
	{
		cleanup();
	}

	std::vector<polynomial_t> const &basis() const { return basis_; }
};

template <typename R, size_t rank> inline void dump(Ideal<R, rank> const &ideal)
{
	fmt::print("polynomial ideal with {} variables and {} equations:\n", rank,
	           ideal.basis().size());
	for (auto &poly : ideal.basis())
		fmt::print("{}\n", poly);
}

} // namespace chalk

#endif
