#pragma once

#include "chalk/monomial.h"
#include "chalk/parser.h"
#include "chalk/polynomial.h"
#include "chalk/rings.h"
#include "fmt/format.h"
#include <array>
#include <bitset>
#include <cassert>
#include <climits>
#include <map>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace chalk {

template <typename R, size_t rank> class SparsePolynomial;

template <typename R, size_t rank> class PolynomialRing
{
	std::array<std::string, rank> var_names_;
	std::array<int, rank> weights_;
	size_t namedVars_ = 0;

  public:
	constexpr PolynomialRing();
	explicit constexpr PolynomialRing(std::array<std::string, rank> var_names);

	// pointers to the ring are stored inside the polynomial, so moving the
	// ring is usually bad, therefore we disable copy/move
	PolynomialRing(PolynomialRing const &) = delete;
	PolynomialRing &operator=(PolynomialRing const &) = delete;

	/** names of the variables */
	auto const &var_names() const { return var_names_; }
	auto const &weights() const { return weights_; }

	/** look up variable name */
	int var_id(std::string const &name) const;

	/** constant polynomial */
	SparsePolynomial<R, rank> operator()(R const &value) const;
	SparsePolynomial<R, rank> operator()(int value) const;

	/** generator x_k */
	SparsePolynomial<R, rank> generator(int k);
	SparsePolynomial<R, rank> generator(std::string_view varName,
	                                    int weight = 1);

	/** parse polynomial from string */
	SparsePolynomial<R, rank> operator()(std::string const &str);
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
	std::vector<Monomial<R, rank>> terms_ = {};

	/** sort terms and collect common terms */
	void cleanup();

  public:
	static PolynomialRing<R, rank> const *default_ring()
	{
		return &trivialRing;
	}

	/** zero polynomial */
	SparsePolynomial() = default;
	explicit SparsePolynomial(PolynomialRing<R, rank> const *r) : ring_(r) {}

	/** constant polynomial */
	explicit SparsePolynomial(int value) : terms_{{R(value), {}}} { cleanup(); }
	explicit SparsePolynomial(R const &value) : terms_{{value, {}}}
	{
		cleanup();
	};
	explicit SparsePolynomial(int value, PolynomialRing<R, rank> const *ring)
	    : ring_(ring), terms_{{R(value), {}}}
	{
		cleanup();
	}
	explicit SparsePolynomial(R const &value,
	                          PolynomialRing<R, rank> const *ring)
	    : ring_(ring), terms_{{value, {}}}
	{
		cleanup();
	};

	/** single monomial */
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

	/** bitset = true for the variables that occur in this polynomial */
	std::bitset<rank> var_occs() const
	{
		std::bitset<rank> seen = {};
		for (auto &term : terms_)
			for (size_t i = 0; i < rank; ++i)
				if (term.exponent[i] != 0)
					seen[i] = true;
		return seen;
	}

	/** number of variables occuring in this polynomials */
	size_t var_count() const { return var_occs().count(); }

	/**
	 * map to a new polynomial ring.
	 *  - coefficients are converted using a callback
	 *  - variables are matched by names
	 *  - reducing rank is okay if some variables dont actually occur
	 */
	template <typename R2, size_t rank2, typename F>
	SparsePolynomial<R2, rank2>
	change_ring(PolynomialRing<R2, rank2> const *ring2, F &&convert) const
	{
		// translation of variable names
		auto trans = std::vector<int>(rank, -1);
		for (int i = 0; i < (int)rank; ++i)
			for (int j = 0; j < (int)rank2; ++j)
				if (ring()->var_names()[i] == ring2->var_names()[j])
				{
					assert(trans[i] == -1);
					trans[i] = j;
				}

		// translate all the terms
		std::vector<Monomial<R2, rank2>> r;
		r.reserve(terms().size());
		for (auto &term : terms())
		{
			r.push_back({convert(term.coefficient), {}});
			for (size_t i = 0; i < rank; ++i)
			{
				if (term.exponent[i] == 0)
					continue;
				if (trans[i] == -1)
					throw std::runtime_error(
					    "SparsePolynomial::change_ring failed");
				assert(r.back().exponent[trans[i]] == 0);
				r.back().exponent[trans[i]] = term.exponent[i];
			}
		}
		return SparsePolynomial<R2, rank2>(ring2, std::move(r));
	}

	bool isConstant() const
	{
		if (terms_.empty())
			return true;
		if (terms_.size() > 1)
			return false;
		for (size_t i = 0; i < rank; ++i)
			if (terms_[0].exponent[i] != 0)
				return false;
		return true;
	}

	/** convert to univariate polynomial, assuming only one variable occurs */
	std::pair<int, Polynomial<R>> to_univariate() const
	{
		std::vector<R> coeffs;
		int pivot = -1;
		for (auto &term : terms_)
		{
			int ex = 0;
			for (int i = 0; i < (int)rank; ++i)
				if (term.exponent[i] != 0)
				{
					if (pivot == -1)
						pivot = i;
					else if (pivot != i)
						throw std::runtime_error(
						    "polynomial with more than one variable");
					assert(ex == 0);
					ex = term.exponent[i];
				}

			if ((int)coeffs.size() < ex + 1)
				coeffs.resize(ex + 1, R(0));
			assert(coeffs[ex] == 0);
			coeffs[ex] = term.coefficient;
		}
		if (pivot == -1)
			pivot = 0;
		return {pivot, Polynomial<R>(std::move(coeffs))};
	}

	/** substitue x_i = val */
	SparsePolynomial substitute(int i, R const &val) const
	{
		auto r = terms();
		for (auto &term : r)
			for (; term.exponent[i] > 0; --term.exponent[i])
				term.coefficient *= val;
		return SparsePolynomial(ring(), std::move(r));
	}
	SparsePolynomial substitute(std::string const &var, R const &val) const
	{
		for (size_t i = 0; i < ring_->var_names().size(); ++i)
			if (ring_->var_names()[i] == var)
				return substitute(i, val);

		// maybe unknown variables should be silently ignored?
		throw std::runtime_error(fmt::format("unknown variable '{}'", var));
	}
	SparsePolynomial substitute(int var, SparsePolynomial const &val) const
	{
		auto c = get_coefficients(*this, var);
		SparsePolynomial r = SparsePolynomial(0, ring());
		for (int i = 0; i < (int)c.size(); ++i)
			r += c[i] * pow(val, i);
		return r;
	}
	R substitute(std::map<std::string, R> const &vals) const
	{
		SparsePolynomial r = *this;
		for (auto const &[var, val] : vals)
			r = r.substitute(var, val);
		assert(r.var_count() == 0);
		assert(r.terms().size() <= 1);
		if (r.terms().size() == 0)
			return R(0);
		else
			return r.terms()[0].coefficient;
	}
	R substitute(std::span<const R> vals) const
	{
		assert(vals.size() == rank);
		auto r = R(0);
		for (auto &term : terms())
		{
			R tmp = term.coefficient;
			for (size_t i = 0; i < rank; ++i)
				if (term.exponent[i] != 0)
					tmp *= pow(vals[i], term.exponent[i]);
			r += tmp;
		}
		return r;
	}

	SparsePolynomial solveFor(int var) const
	{
		auto c = get_coefficients(*this, var);
		assert(c.size() == 2);
		assert(c[1].isConstant());
		return -c[0] / c[1].terms()[0].coefficient;
	}

	/** largest power of x_k */
	int max_order(size_t k) const
	{
		assert(k < rank);
		int r = 0;
		for (auto const &term : terms())
			r = std::max(r, term.exponent[k]);
		return r;
	}

	/** lead coefficient and exponent and monomial*/
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
	SparsePolynomial lm() const
	{
		if (terms_.empty())
			return SparsePolynomial(ring_);
		return SparsePolynomial(ring_, {terms_[0]});
	}
	SparsePolynomial trailing_terms(size_t count) const
	{
		std::vector<Monomial<R, rank>> r;
		for (size_t i = terms_.size() >= count ? terms_.size() - count : 0;
		     i < terms_.size(); ++i)
			r.push_back(terms_[i]);
		return SparsePolynomial(ring_, r);
	}
	R coefficient(std::array<int, rank> ex) const
	{
		for (auto &term : terms_)
			if (term.exponent == ex)
				return term.coefficient;
		return R(0);
	}

	/** "optimized" inplace operations */
	void operator+=(SparsePolynomial<R, rank> const &b)
	{
		assert(this != &b);
		terms_.reserve(terms_.size() + b.terms().size());
		for (auto &term : b.terms())
			terms_.push_back(term);
		cleanup();
	}
	void operator-=(SparsePolynomial<R, rank> const &b)
	{
		assert(this != &b);
		terms_.reserve(terms_.size() + b.terms().size());
		for (auto &term : b.terms())
			terms_.push_back({-term.coefficient, term.exponent});
		cleanup();
	}
	void operator*=(Monomial<R, rank> const b)
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

	/** only division by R is supported (TODO: polynomial divison?) */
	void operator/=(SparsePolynomial<R, rank> const &b)
	{
		assert(this != &b);
		if (b.terms().size() == 0)
			throw std::runtime_error("tried to divide polynomial by zero");
		if (b.terms().size() != 1 ||
		    b.terms()[0].exponent != std::array<int, rank>{})
			throw std::runtime_error(
			    "tried to divide polynomial by non-constant polynomial");
		(*this) /= b.terms()[0].coefficient;
	}
};

///////////////////////////////////////////////////////////////////////////////
// definitions for PolynomialRing
///////////////////////////////////////////////////////////////////////////////

template <typename R, size_t rank>
constexpr PolynomialRing<R, rank>::PolynomialRing()
{
	namedVars_ = 0;
	for (size_t i = 0; i < rank; ++i)
	{
		var_names_[i] = fmt::format("x_{}", i);
		weights_[i] = 1;
	}
}

template <typename R, size_t rank>
constexpr PolynomialRing<R, rank>::PolynomialRing(
    std::array<std::string, rank> var_names)
    : var_names_(std::move(var_names))
{
	namedVars_ = rank;
	for (size_t i = 0; i < rank; ++i)
		weights_[i] = 1;
}

template <typename R, size_t rank>
int PolynomialRing<R, rank>::var_id(std::string const &name) const
{
	for (size_t i = 0; i < rank; ++i)
		if (var_names_[i] == name)
			return (int)i;
	throw std::runtime_error(
	    fmt::format("Variable name '{}' not found in Poly-Ring", name));
}

template <typename R, size_t rank>
SparsePolynomial<R, rank>
PolynomialRing<R, rank>::operator()(std::string const &str)
{
	auto f = [&](std::string_view token) -> SparsePolynomial<R, rank> {
		if (std::isdigit(token[0]))
			return (*this)(util::parse_int(token));
		else
			return generator(token);
	};

	return parse<SparsePolynomial<R, rank>>(str, f);
}

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

///////////////////////////////////////////////////////////////////////////////
// definitions for SparsePolynomial
///////////////////////////////////////////////////////////////////////////////

template <typename R, size_t rank> void SparsePolynomial<R, rank>::cleanup()
{
	std::sort(terms_.begin(), terms_.end(),
	          [&ws = ring_->weights()](auto &a, auto &b) {
		          return order_degrevlex<R, rank>(a, b, ws);
	          });
	size_t j = 0;
	for (size_t i = 0; i < terms_.size(); ++i)
	{
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
	}
	terms_.resize(j);
}

template <typename R, size_t rank>
bool need_parens(SparsePolynomial<R, rank> const &poly)
{
	if (poly.terms().size() == 0)
		return false;
	if (poly.terms().size() == 1)
		return need_parens(poly.terms()[0].coefficient);
	else
		return true;
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

template <typename R, size_t rank>
SparsePolynomial<R, rank> pow(SparsePolynomial<R, rank> a, int b)
{
	if (b < 0)
		throw std::runtime_error("cannot raise polynomial to negative power");

	auto r = SparsePolynomial<R, rank>(1);
	for (; b != 0; b >>= 1, a *= a)
		if (b & 1)
			r *= a;
	return r;
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
SparsePolynomial<R, rank>
PolynomialRing<R, rank>::operator()(R const &value) const
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
PolynomialRing<R, rank>::generator(std::string_view varName, int weight)
{
	assert(weight >= 1);
	for (size_t i = 0; i < namedVars_; ++i)
		if (var_names_[i] == varName)
		{
			// assert(weights_[i] == weight);
			return generator(i);
		}
	if (namedVars_ >= rank)
		throw std::runtime_error("too many generators in polynomial ring");
	var_names_[namedVars_] = varName;
	weights_[namedVars_] = weight;
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

template <typename R, size_t rank, typename F>
SparsePolynomial<R, rank> mapCoefficients(F f,
                                          SparsePolynomial<R, rank> const &a)
{
	// NOTE: the ring knows the coefficient type. therefore we cant change
	//       the coefficient type.

	std::vector<Monomial<R, rank>> terms;
	terms.reserve(a.terms().size());
	for (auto &term : a.terms())
		terms.emplace_back(f(term.coefficient), term.exponent);
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

/** get coefficients of x_k */
template <typename R, size_t rank>
std::vector<SparsePolynomial<R, rank>>
get_coefficients(SparsePolynomial<R, rank> const &a, size_t k)
{
	auto m = a.max_order(k);
	std::vector<SparsePolynomial<R, rank>> r;
	r.reserve(m + 1);
	for (int i = 0; i <= m; ++i)
		r.push_back(get_coefficient(a, k, i));
	return r;
}

/** differentiate w.r.t. x_k */
template <typename R, size_t rank>
SparsePolynomial<R, rank> diff(SparsePolynomial<R, rank> const &a, size_t k)
{
	auto terms = a.terms();
	for (auto &t : terms)
	{
		t.coefficient *= t.exponent[k];
		t.exponent[k] -= 1;
	}
	return SparsePolynomial<R, rank>(a.ring(), std::move(terms));
}

template <typename R, size_t rank>
SparsePolynomial<R, rank> diff(SparsePolynomial<R, rank> const &a,
                               std::string const &var)
{
	return diff(a, a.ring()->var_id(var));
}

template <typename R, size_t rank> struct RingTraits<SparsePolynomial<R, rank>>
{
	static bool is_zero(SparsePolynomial<R, rank> const &a)
	{
		return a.terms().empty();
	}
	static bool is_one(SparsePolynomial<R, rank> const &a) { return a == 1; }
	static bool is_negative(SparsePolynomial<R, rank> const &) { return false; }

	/** these two are not precisely correct, but at least not wrong */
	static bool need_parens_product(SparsePolynomial<R, rank> const &)
	{
		return true;
	}
	static bool need_parens_power(SparsePolynomial<R, rank> const &)
	{
		return true;
	}
};

} // namespace chalk

template <typename R, size_t rank>
struct fmt::formatter<chalk::SparsePolynomial<R, rank>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::SparsePolynomial<R, rank> &poly,
	            FormatContext &ctx) const -> decltype(ctx.out())
	{
		auto const &terms = poly.terms();

		// no terms -> output "0"
		if (terms.empty())
			return fmt::format_to(ctx.out(), "0");

		// otherwise list the terms with " + " inbetween
		auto it = ctx.out();
		for (size_t i = 0; i < terms.size(); ++i)
		{
			R coeff = terms[i].coefficient;
			auto &exp = terms[i].exponent;

			if (chalk::is_negative(coeff))
			{
				if (i == 0)
					it = fmt::format_to(it, "-");
				else
					it = fmt::format_to(it, " - ");
				coeff = -coeff;
			}
			else if (i != 0)
				it = fmt::format_to(it, " + ");

			bool isOne = coeff == 1;
			if (!isOne || exp == std::array<int, rank>{})
				it = fmt::format_to(it, "{}", coeff);
			bool first = true;
			for (size_t k = 0; k < rank; ++k)
			{
				if (exp[k] == 0)
					continue;
				if (!first || !isOne)
					*it++ = '*';
				first = false;
				if (exp[k] == 1)
					it = fmt::format_to(it, "{}", poly.ring()->var_names()[k]);
				else
					it = fmt::format_to(it, "{}^{}",
					                    poly.ring()->var_names()[k], exp[k]);
			}
		}
		return it;
	}
};
