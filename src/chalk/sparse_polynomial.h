#ifndef CHALK_SPARSE_POLYNOMIAL_H
#define CHALK_SPARSE_POLYNOMIAL_H

#include "chalk/polynomial.h"
#include "fmt/format.h"
#include "util/span.h"
#include <array>
#include <cassert>
#include <climits>
#include <utility>
#include <vector>

namespace chalk {

/** a single term of a polynomial: a*x1^e1*...*xn^en */
template <typename R, size_t rank> struct Monomial
{
	R coefficient = {};
	std::array<int, rank> exponent = {};

	Monomial operator*(Monomial const &b) const
	{
		std::array<int, rank> ex;
		for (size_t i = 0; i < rank; ++i)
			ex[i] = exponent[i] + b.exponent[i];
		return Monomial{coefficient * b.coefficient, ex};
	}
	std::optional<Monomial> operator/(Monomial const &b) const
	{
		std::array<int, rank> ex;
		for (size_t i = 0; i < rank; ++i)
			if (exponent[i] >= b.exponent[i])
				ex[i] = exponent[i] - b.exponent[i];
			else
				return std::nullopt;
		return Monomial{coefficient / b.coefficient, ex};
	}
};

template <typename R, size_t rank>
bool order_degrevlex(Monomial<R, rank> const &a, Monomial<R, rank> const &b)
{
	int totA = 0, totB = 0;
	for (auto e : a.exponent)
		totA += e;
	for (auto e : b.exponent)
		totB += e;
	if (totA != totB)
		return totA > totB;

	for (size_t i = 0; i < rank; ++i)
		if (a.exponent[i] < b.exponent[i])
			return true;
		else if (a.exponent[i] > b.exponent[i])
			return false;
	return false;
}

template <typename R, size_t rank> class PolynomialRing;
template <typename R, size_t rank> class SparsePolynomial;

template <typename R, size_t rank> class PolynomialRing
{
	std::array<std::string, rank> varNames_;
	std::array<int, rank> max_order_;
	size_t namedVars_ = 0;

  public:
	constexpr PolynomialRing()
	{
		namedVars_ = 0;
		for (size_t i = 0; i < rank; ++i)
		{
			varNames_[i] = fmt::format("x_{}", i);
			max_order_[i] = INT_MAX;
		}
	}
	explicit constexpr PolynomialRing(std::array<std::string, rank> varNames)
	    : varNames_(std::move(varNames))
	{
		namedVars_ = rank;
		for (size_t i = 0; i < rank; ++i)
			max_order_[i] = INT_MAX;
	}

	// pointers to the ring are stored inside the polynomial, so moving the
	// ring is usually bad, therefore we disable copy/move
	PolynomialRing(PolynomialRing const &) = delete;
	PolynomialRing &operator=(PolynomialRing const &) = delete;

	/** names of the variables */
	auto const &varNames() const { return varNames_; }
	auto const &max_order() const { return max_order_; }

	/** constant polynomial */
	SparsePolynomial<R, rank> operator()(R const &value) const;
	SparsePolynomial<R, rank> operator()(int value) const;

	/** generator x_k */
	SparsePolynomial<R, rank> generator(int k);
	SparsePolynomial<R, rank> generator(std::string const &varName,
	                                    int maxOrder = INT_MAX);
};

/**
 * multivariate sparse polynomial.
 *   - term-order = lexicographic
 */
template <typename R, size_t rank> class SparsePolynomial
{
  public:
	inline static const PolynomialRing<R, rank> trivialRing = {};

  private:
	PolynomialRing<R, rank> const *ring_ = &trivialRing;
	std::vector<Monomial<R, rank>> terms_;

	/** sort terms and collect common terms */
	void cleanup()
	{
		std::sort(terms_.begin(), terms_.end(), order_degrevlex<R, rank>);
		size_t j = 0;
		for (size_t i = 0; i < terms_.size(); ++i)
		{
			// effectively remove terms of too-high order
			for (size_t k = 0; k < rank; ++k)
				if (terms_[i].exponent[k] > ring_->max_order()[k])
					goto next_term;

			// if exponent == previous exponent -> join
			if (j != 0 && terms_[j - 1].exponent == terms_[i].exponent)
				terms_[j - 1].coefficient += terms_[i].coefficient;

			// else new term
			else if (i != j)
				terms_[j++] = std::move(terms_[i]);
			else
				j++;

			// if coeff = 0 -> remove term
			assert(j);
			if (terms_[j - 1].coefficient == 0)
				--j;

		next_term:;
		}
		terms_.resize(j);
	}

  public:
	/** zero polynomial */
	SparsePolynomial() = default;
	explicit SparsePolynomial(PolynomialRing<R, rank> const *r) : ring_(r) {}

	/** constant polynomial */
	explicit SparsePolynomial(int value) : terms_{{R(value), {}}} { cleanup(); }
	explicit SparsePolynomial(R const &value) : terms_{{value, {}}}
	{
		cleanup();
	};

	/** monomial */
	explicit SparsePolynomial(std::array<int, rank> const &exponent)
	    : terms_{{R(1), exponent}}
	{}

	/** general constructor */
	explicit SparsePolynomial(PolynomialRing<R, rank> const *r,
	                          std::vector<Monomial<R, rank>> ts)
	    : ring_(r), terms_(std::move(ts))
	{
		cleanup();
	}

	/** (read-only) access to the ring and the terms */
	PolynomialRing<R, rank> const *ring() const { return ring_; }
	std::vector<Monomial<R, rank>> const &terms() const { return terms_; }

	/** lead coefficient and exponent */
	R const &lc() const
	{
		assert(!terms_.empty());
		return terms_[0].coefficient;
	}
	std::array<int, rank> const &le() const
	{
		assert(!terms_.empty());
		return terms_[0].exponent;
	}

	/** "optimized" inplace operations */
	void operator+=(SparsePolynomial<R, rank> const &b)
	{
		terms_.reserve(terms_.size() + b.terms().size());
		for (auto &term : b.terms())
			terms_.push_back(term);
		cleanup();
	}
	void operator-=(SparsePolynomial<R, rank> const &b)
	{
		terms_.reserve(terms_.size() + b.terms().size());
		for (auto &term : b.terms())
			terms_.push_back({-term.coefficient, term.exponent});
		cleanup();
	}
	void operator*=(Monomial<R, rank> const &b)
	{
		for (auto &term : terms_)
			term *= b;

		cleanup(); // this does nothing if R and the term-order are nice
	}
	void operator*=(R b) // dont take b by reference due to aliasing problems
	{
		for (auto &term : terms_)
			term.coefficient *= b;
		cleanup();
	}
	void operator/=(R b) // dont take b by reference due to aliasing problems
	{
		for (auto &term : terms_)
			term.coefficient /= b;
		cleanup();
	}
	void operator*=(int b)
	{
		for (auto &term : terms_)
			term.coefficient *= b;
		cleanup();
	}
	void operator/=(int b)
	{
		for (auto &term : terms_)
			term.coefficient /= b;
		cleanup();
	}
};

/** helper function to determine ring of result of binary operations */
template <typename R, size_t rank>
PolynomialRing<R, rank> const *unify(PolynomialRing<R, rank> const *a,
                                     PolynomialRing<R, rank> const *b)
{
	if (a == b)
		return a;
	if (a == &SparsePolynomial<R, rank>::trivialRing)
		return b;
	if (b == &SparsePolynomial<R, rank>::trivialRing)
		return a;
	throw std::runtime_error("incompatible polynomial rings");
}

/** unary operations */
template <typename R, size_t rank>
SparsePolynomial<R, rank> operator-(SparsePolynomial<R, rank> const &a)
{
	std::vector<Monomial<R, rank>> terms = a.terms();
	for (auto &t : terms)
		t.coefficient = -t.coefficient;
	return SparsePolynomial(a.ring(), std::move(terms));
}

/** binary operations (Poly <-> int) */
template <typename R, size_t rank>
SparsePolynomial<R, rank> operator+(SparsePolynomial<R, rank> const &a, int b)
{
	auto terms = a.terms();
	terms.push_back({R(b), {}});
	return SparsePolynomial(a.ring(), std::move(terms));
}
template <typename R, size_t rank>
SparsePolynomial<R, rank> operator-(SparsePolynomial<R, rank> const &a, int b)
{
	auto terms = a.terms();
	terms.push_back({R(-b), {}});
	return SparsePolynomial(a.ring(), std::move(terms));
}
template <typename R, size_t rank>
SparsePolynomial<R, rank> operator*(SparsePolynomial<R, rank> const &a, int b)
{
	std::vector<Monomial<R, rank>> terms;
	terms.reserve(a.terms().size());
	for (auto &t : a.terms())
		terms.push_back({t.coefficient * b, t.exponent});
	return SparsePolynomial(a.ring(), std::move(terms));
}
template <typename R, size_t rank>
SparsePolynomial<R, rank> operator/(SparsePolynomial<R, rank> const &a, int b)
{
	std::vector<Monomial<R, rank>> terms;
	terms.reserve(a.terms().size());
	for (auto &t : a.terms())
		terms.push_back({t.coefficient / b, t.exponent});
	return SparsePolynomial(a.ring(), std::move(terms));
}

/** binary operations (Poly <-> R) */
template <typename R, size_t rank>
SparsePolynomial<R, rank> operator*(SparsePolynomial<R, rank> const &a,
                                    R const &b)
{
	std::vector<Monomial<R, rank>> terms;
	terms.reserve(a.terms().size());
	for (auto &t : a.terms())
		terms.push_back({t.coefficient * b, t.exponent});
	return SparsePolynomial(a.ring(), std::move(terms));
}
template <typename R, size_t rank>
SparsePolynomial<R, rank> operator/(SparsePolynomial<R, rank> const &a,
                                    R const &b)
{
	std::vector<Monomial<R, rank>> terms;
	terms.reserve(a.terms().size());
	for (auto &t : a.terms())
		terms.push_back({t.coefficient / b, t.exponent});
	return SparsePolynomial(a.ring(), std::move(terms));
}

/** binary operations (Poly <-> Monomial) */
template <typename R, size_t rank>
SparsePolynomial<R, rank> operator*(SparsePolynomial<R, rank> const &a,
                                    Monomial<R, rank> const &b)
{
	std::vector<Monomial<R, rank>> terms;
	terms.reserve(a.terms().size());
	for (auto &t : a.terms())
		terms.push_back(t * b);
	return SparsePolynomial(a.ring(), std::move(terms));
}

/** binary operations (Poly <-> Poly) */
template <typename R, size_t rank>
SparsePolynomial<R, rank> operator+(SparsePolynomial<R, rank> const &a,
                                    SparsePolynomial<R, rank> const &b)
{
	std::vector<Monomial<R, rank>> terms;
	terms.reserve(a.terms().size() + b.terms().size());
	for (auto &t : a.terms())
		terms.push_back(t);
	for (auto &t : b.terms())
		terms.push_back(t);
	return SparsePolynomial(unify(a.ring(), b.ring()), std::move(terms));
}

template <typename R, size_t rank>
SparsePolynomial<R, rank> operator-(SparsePolynomial<R, rank> const &a,
                                    SparsePolynomial<R, rank> const &b)
{
	std::vector<Monomial<R, rank>> terms;
	terms.reserve(a.terms().size() + b.terms().size());
	for (auto &t : a.terms())
		terms.push_back(t);
	for (auto &t : b.terms())
		terms.push_back({-t.coefficient, t.exponent});
	return SparsePolynomial(unify(a.ring(), b.ring()), std::move(terms));
}

template <typename R, size_t rank>
SparsePolynomial<R, rank> operator*(SparsePolynomial<R, rank> const &a,
                                    SparsePolynomial<R, rank> const &b)
{
	std::vector<Monomial<R, rank>> terms;
	terms.reserve(a.terms().size() * b.terms().size());
	for (auto &t : a.terms())
		for (auto &u : b.terms())
			terms.push_back(t * u);
	return SparsePolynomial(unify(a.ring(), b.ring()), std::move(terms));
}

/** comparisons */
template <typename R, size_t rank>
bool operator==(SparsePolynomial<R, rank> const &a,
                SparsePolynomial<R, rank> const &b)
{
	unify(a.ring(), b.ring());
	return a.terms() == b.terms();
}
template <typename R, size_t rank>
bool operator==(SparsePolynomial<R, rank> const &a, int b)
{
	if (b == 0)
		return a.terms().empty();
	if (a.terms().size() != 1)
		return false;
	if (!(a.terms()[0].coefficient == b))
		return false;
	for (auto k : a.terms()[0].exponent)
		if (k != 0)
			return false;
	return true;
}

/** op-assign for convenience */
template <typename R, size_t rank>
void operator*=(SparsePolynomial<R, rank> &a,
                SparsePolynomial<R, rank> const &b)
{
	a = a * b;
}

/** constant polynomial */
template <typename R, size_t rank>
SparsePolynomial<R, rank> PolynomialRing<R, rank>::
operator()(R const &value) const
{
	std::vector<Monomial<R, rank>> terms;
	terms.emplace_back(value, {});
	return SparsePolynomial<R, rank>(this, std::move(terms));
}

template <typename R, size_t rank>
SparsePolynomial<R, rank> PolynomialRing<R, rank>::operator()(int value) const
{
	std::vector<Monomial<R, rank>> terms;
	terms.push_back({R(value), std::array<int, rank>{}});
	return SparsePolynomial<R, rank>(this, std::move(terms));
}

/** generator x_k */
template <typename R, size_t rank>
SparsePolynomial<R, rank> PolynomialRing<R, rank>::generator(int k)
{
	std::array<int, rank> exp = {};
	exp[k] = 1;
	std::vector<Monomial<R, rank>> terms;
	terms.push_back({R(1), exp});
	return SparsePolynomial<R, rank>(this, std::move(terms));
}
template <typename R, size_t rank>
SparsePolynomial<R, rank>
PolynomialRing<R, rank>::generator(std::string const &varName, int max_order)
{
	assert(max_order >= 0);
	for (size_t i = 0; i < namedVars_; ++i)
		if (varNames_[i] == varName)
		{
			assert(max_order_[i] == max_order);
			return generator(i);
		}
	if (namedVars_ >= rank)
		throw std::runtime_error("too many generators in polynomial ring");
	varNames_[namedVars_] = varName;
	max_order_[namedVars_] = max_order;
	return generator(namedVars_++);
}

/** remove all odd powers of x_k */
template <typename R, size_t rank>
SparsePolynomial<R, rank> remove_odd_powers(SparsePolynomial<R, rank> const &a,
                                            size_t k)
{
	assert(0 <= k && k < rank);
	auto terms = a.terms();
	for (auto &t : terms)
		if (t.exponent[k] % 2)
			t.coefficient = R(0);
	return SparsePolynomial<R, rank>(a.ring(), std::move(terms));
}

/** remove all powers above x_k^n */
template <typename R, size_t rank>
SparsePolynomial<R, rank> truncate(SparsePolynomial<R, rank> const &a, size_t k,
                                   int n)
{
	assert(0 <= k && k < rank);
	auto terms = a.terms();
	for (auto &t : terms)
		if (t.exponent[k] > n)
			t.coefficient = R(0);
	return SparsePolynomial<R, rank>(a.ring(), std::move(terms));
}

/** get coefficient of x_k^n */
template <typename R, size_t rank>
SparsePolynomial<R, rank> get_coefficient(SparsePolynomial<R, rank> const &a,
                                          size_t k, int n)
{
	assert(0 <= k && k < rank);
	auto terms = a.terms();
	for (auto &t : terms)
		if (t.exponent[k] == n)
			t.exponent[k] = 0;
		else
			t.coefficient = R(0);
	return SparsePolynomial<R, rank>(a.ring(), std::move(terms));
}

} // namespace chalk

template <typename R> struct fmt::formatter<chalk::Polynomial<R>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Polynomial<R> &poly, FormatContext &ctx)
	    -> decltype(ctx.out())
	{
		// no terms -> output "0"
		if (poly.coefficients().empty())
			return format_to(ctx.out(), "0");

		// otherwise list the terms with " + " inbetween
		auto it = ctx.out();
		for (int i = poly.degree(); i >= 0; --i)
		{
			auto c = poly[i];
			if (c == 0)
				continue;

			if (is_negative(c))
			{
				if (i == poly.degree())
					it = format_to(it, "-");
				else
					it = format_to(it, " - ");
				c = -c;
			}
			else if (i != poly.degree())
				it = format_to(it, " + ");

			if (i == 0)
				it = format_to(it, "{}", c);
			else if (!(c == 1))
				it = format_to(it, "{}*", c);

			if (i == 1)
				it = format_to(it, "x");
			else if (i >= 2)
				it = format_to(it, "x^{}", i);
		}
		return it;
	}
};

template <typename R, size_t rank>
struct fmt::formatter<chalk::SparsePolynomial<R, rank>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::SparsePolynomial<R, rank> &poly,
	            FormatContext &ctx) -> decltype(ctx.out())
	{
		auto const &terms = poly.terms();

		// no terms -> output "0"
		if (terms.empty())
			return format_to(ctx.out(), "0");

		// otherwise list the terms with " + " inbetween
		auto it = ctx.out();
		for (size_t i = 0; i < terms.size(); ++i)
		{
			R coeff = terms[i].coefficient;
			auto &exp = terms[i].exponent;

			if (is_negative(coeff))
			{
				if (i == 0)
					it = format_to(it, "-");
				else
					it = format_to(it, " - ");
				coeff = -coeff;
			}
			else if (i != 0)
				it = format_to(it, " + ");

			bool isOne = coeff == 1;
			if (!isOne || exp == std::array<int, rank>{})
				it = format_to(it, "{}", coeff);
			bool first = true;
			for (size_t k = 0; k < rank; ++k)
			{
				if (exp[k] == 0)
					continue;
				if (!first || !isOne)
					it++ = '*';
				first = false;
				if (exp[k] == 1)
					it = format_to(it, "{}", poly.ring()->varNames()[k]);
				else
					it = format_to(it, "{}^{}", poly.ring()->varNames()[k],
					               exp[k]);
			}
		}
		return it;
	}
};

#endif
