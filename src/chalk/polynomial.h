#ifndef CHALK_POLYNOMIAL_H
#define CHALK_POLYNOMIAL_H

#include "fmt/format.h"
#include "util/span.h"
#include "util/vector2d.h"
#include <array>
#include <cassert>
#include <utility>
#include <vector>

namespace chalk {

template <typename R, size_t rank> class PolynomialRing;
template <typename R, size_t rank> class Monomial;
template <typename R, size_t rank> class SparsePolynomial;

template <typename R, size_t rank> class PolynomialRing
{
	std::array<std::string, rank> varNames_;

  public:
	constexpr PolynomialRing()
	{
		for (size_t i = 0; i < rank; ++i)
			varNames_[i] = fmt::format("x_{}", i);
	}
	explicit constexpr PolynomialRing(std::array<std::string, rank> varNames)
	    : varNames_(std::move(varNames))
	{}

	// pointers to the ring are stored inside the polynomial, so moving the
	// ring is usually bad, therefore we disable copy/move
	PolynomialRing(PolynomialRing const &) = delete;
	PolynomialRing &operator=(PolynomialRing const &) = delete;

	/** names of the variables */
	auto const &varNames() const { return varNames_; }

	/** constant polynomial */
	SparsePolynomial<R, rank> operator()(R const &value) const;
	SparsePolynomial<R, rank> operator()(int value) const;

	/** generator x_k */
	SparsePolynomial<R, rank> generator(int k);
};

template <typename R, size_t rank> struct Monomial
{
	R coefficient = {};
	std::array<int, rank> exponent = {};
};

template <typename R, size_t rank>
bool order_lex(Monomial<R, rank> const &a, Monomial<R, rank> const &b)
{
	for (size_t i = 0; i < rank; ++i)
		if (a.exponent[i] > b.exponent[i])
			return true;
		else if (a.exponent[i] < b.exponent[i])
			return false;
	return false;
}

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
		std::sort(terms_.begin(), terms_.end(), order_lex<R, rank>);
		size_t j = 0;
		for (size_t i = 0; i < terms_.size(); ++i)
		{
			if (j != 0 && terms_[j - 1].exponent == terms_[i].exponent)
				terms_[j - 1].coefficient += terms_[i].coefficient;
			else if (i != j)
				terms_[j++] = std::move(terms_[i]);
			else
				j++;

			if (terms_[j - 1].coefficient == 0)
				--j;
		}
		terms_.resize(j);
	}

  public:
	SparsePolynomial() = default;
	explicit SparsePolynomial(int value) : terms_{{R(value), {}}} { cleanup(); }
	explicit SparsePolynomial(R const &value) : terms_{{value, {}}}
	{
		cleanup();
	};
	explicit SparsePolynomial(PolynomialRing<R, rank> const *r) : ring_(r) {}
	explicit SparsePolynomial(PolynomialRing<R, rank> const *r,
	                          std::vector<Monomial<R, rank>> ts)
	    : ring_(r), terms_(std::move(ts))
	{
		cleanup();
	}

	/** (read-only) access to the ring and the terms */
	PolynomialRing<R, rank> const *ring() const { return ring_; }
	std::vector<Monomial<R, rank>> const &terms() const { return terms_; }
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
		{
			std::array<int, rank> mon;
			for (size_t k = 0; k < rank; ++k)
				mon[k] = t.exponent[k] + u.exponent[k];
			terms.push_back({t.coefficient * u.coefficient, mon});
		}
	return SparsePolynomial(unify(a.ring(), b.ring()), std::move(terms));
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

} // namespace chalk

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
