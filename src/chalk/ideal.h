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
		if (auto d = f.terms()[i] / g.terms()[0]; d)
		{
			f -= g * d.value();
			--i;
			change = true;
		}
	return change;
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
			return order_degrevlex(a.terms()[0], b.terms()[0]);
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

	std::vector<polynomial_t> const &basis() const { return basis_; }

	/** transform into groebner basis */
	void groebner(size_t max_polys = SIZE_MAX);
};

template <typename R, size_t rank> inline void dump(Ideal<R, rank> const &ideal)
{
	fmt::print("polynomial ideal with {} variables and {} equations:\n", rank,
	           ideal.basis().size());
	for (auto &poly : ideal.basis())
		fmt::print("{}\n", poly);
}

template <typename R, size_t rank>
inline void Ideal<R, rank>::groebner(size_t max_polys)
{
	// really naive algorithm adding resolvents until saturation
	while (true)
	{
		bool change = false;
		for (size_t i = 0; i < basis_.size(); ++i)
			for (size_t j = i + 1;
			     j < basis_.size() && basis_.size() < max_polys; ++j)
			{
				auto r = resolvent(basis_[i], basis_[j]);
				reduce(r);
				if (r.terms().empty())
					continue;

				change = true;
				basis_.push_back(r);
			}
		if (change)
			cleanup();
		else
			break;
	}
}

} // namespace chalk

#endif
