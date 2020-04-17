#ifndef CHALK_IDEAL_H
#define CHALK_IDEAL_H

/** ideals in multivariate polynomial rings */

#include "chalk/polynomial.h"

namespace chalk {

template <typename R, size_t rank> class Ideal
{
  public:
	using polynomial_ring_t = PolynomialRing<R, rank>;
	using polynomial_t = SparsePolynomial<R, rank>;

  private:
	std::vector<polynomial_t> basis_;

	static bool divides(std::array<int, rank> const &a,
	                    std::array<int, rank> const &b)
	{
		for (size_t i = 0; i < rank; ++i)
			if (a[i] > b[i])
				return false;
		return true;
	}

	static std::array<int, rank> divide(std::array<int, rank> const &a,
	                                    std::array<int, rank> const &b)
	{
		std::array<int, rank> c;
		for (size_t i = 0; i < rank; ++i)
		{
			assert(a[i] <= b[i]);
			c[i] = b[i] - a[i];
		}
		return c;
	}

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
				{
					if (i == j)
						continue;

					for (auto &t : basis_[j].terms())
						if (divides(basis_[i].terms()[0].exponent, t.exponent))
						{
							auto tmp = basis_[i] * polynomial_t(t.coefficient);
							basis_[j] -=
							    basis_[i] * polynomial_t(t.coefficient) *
							    polynomial_t(divide(
							        basis_[i].terms()[0].exponent, t.exponent));
							change = true;

							// 't' is dangling now, so get out asap
							break;
						}
				}
			}
		} while (change == true);

		// remove zero entries
		basis_.erase(
		    std::remove_if(basis_.begin(), basis_.end(),
		                   [](polynomial_t const &poly) { return poly == 0; }),
		    basis_.end());
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
