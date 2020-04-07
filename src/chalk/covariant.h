#ifndef CHALK_COVARIANT
#define CHALK_COVARIANT

#include "fmt/format.h"
#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace chalk {

/** a (named) object with indices */
struct CovariantAtom
{
	std::string symbol;
	std::vector<int> indices; // first internal, then derivatives
};

inline bool operator==(CovariantAtom const &a, CovariantAtom const &b)
{
	return a.symbol == b.symbol && a.indices == b.indices;
}

inline bool operator<(CovariantAtom const &a, CovariantAtom const &b)
{
	// ascending name
	if (a.symbol < b.symbol)
		return true;
	if (a.symbol > b.symbol)
		return false;

	// descending number of indices
	return a.indices.size() > b.indices.size();
}

/** product of atoms and a pre-factor */
template <typename R> struct CovariantTerm
{
	R coefficient;
	std::vector<CovariantAtom> atoms;

	/**
	 * Brings atoms into partial normal form:
	 *   - sort atoms
	 *   - rename summed indices
	 *   - return list of outer indices
	 */
	std::vector<int> cleanup()
	{
		// sort atoms by symbol and number of indices
		// (this assumes all atoms commute)
		std::sort(atoms.begin(), atoms.end());

		// count occurences of indices (should be 1 or 2)
		std::map<int, int> count;
		for (auto &a : atoms)
			for (int k : a.indices)
				count[k] += 1;
		std::vector<int> open;
		for (auto &[k, c] : count)
		{
			if (c == 1)
				open.push_back(k);
			else if (c == 2)
				continue;
			else
				throw std::runtime_error("inconsistent indices");
		}
		std::sort(open.begin(), open.end());

		// renmae inner indices
		int z = open.empty() ? 1 : open.back() + 1; // first inner index

		std::map<int, int> trans;
		for (auto &a : atoms)
			for (int &k : a.indices)
			{
				if (count[k] != 2)
					continue;
				if (trans.count(k) == 0)
					trans[k] = z++;
				k = trans[k];
			}

		return open;
	}
};

template <typename R>
inline bool operator<(CovariantTerm<R> const &a, CovariantTerm<R> const &b)
{
	return a.atoms < b.atoms;
}

/**
 * A covariant expression is a sum of covariant terms.
 * These do not quite form a ring, due to the indices
 */
template <typename R> class Covariant
{
  public:
	using Atom = CovariantAtom;
	using Term = CovariantTerm<R>;

  private:
	std::vector<Term> terms_;
	int first_inner_index_ = 1;

	void cleanup()
	{
		if (terms_.size() == 0)
			return;

		// cleanup individual terms
		std::vector<int> open = terms_[0].cleanup();
		first_inner_index_ = open.empty() ? 1 : open.back() + 1;
		for (size_t i = 1; i < terms_.size(); ++i)
			if (open != terms_[i].cleanup())
				throw std::runtime_error("inconsistent indices");

		// sort and collect terms
		std::sort(terms_.begin(), terms_.end());
		size_t j = 0;
		for (size_t i = 0; i < terms_.size(); ++i)
		{
			// same atoms as previous -> just add coefficients
			if (i != 0 && terms_[i].atoms == terms_[j - 1].atoms)
				terms_[j - 1].coefficient += terms_[i].coefficient;

			// otherwise -> term
			else
			{
				if (j != 0 && terms_[j - 1].coefficient == 0)
					--j;

				terms_[j++] = terms_[i];
			}
		}
		terms_.resize(j);
	}

  public:
	Covariant() = default;
	explicit Covariant(std::vector<Term> terms) : terms_(std::move(terms))
	{
		cleanup();
	}

	/** constructor from base ring */
	explicit Covariant(int a) : terms_{Term{R(a), {}}} { cleanup(); }
	explicit Covariant(R const &a) : terms_{Term{a, {}}} { cleanup(); }

	/** constructor from single symbol with indices */
	explicit Covariant(std::string const &symbol)
	    : terms_{Term{R(1), {Atom{symbol, {}}}}}
	{
		cleanup();
	}
	explicit Covariant(std::string const &symbol, int a)
	    : terms_{Term{R(1), {Atom{symbol, {a}}}}}
	{
		cleanup();
	}
	explicit Covariant(std::string const &symbol, int a, int b)
	    : terms_{Term{R(1), {Atom{symbol, {a, b}}}}}
	{
		cleanup();
	}
	explicit Covariant(std::string const &symbol, int a, int b, int c)
	    : terms_{Term{R(1), {Atom{symbol, {a, b, c}}}}}
	{
		cleanup();
	}
	explicit Covariant(std::string const &symbol, int a, int b, int c, int d)
	    : terms_{Term{R(1), {Atom{symbol, {a, b, c, d}}}}}
	{
		cleanup();
	}

	std::vector<Term> const &terms() const { return terms_; }
	int first_inner_index() const { return first_inner_index_; }
};

template <typename R>
Covariant<R> inline operator+(Covariant<R> const &a, Covariant<R> const &b)
{
	std::vector<CovariantTerm<R>> terms = a.terms();
	for (auto &t : b.terms())
		terms.push_back(t);
	return Covariant<R>(std::move(terms));
}
template <typename R>
Covariant<R> inline operator-(Covariant<R> const &a, Covariant<R> const &b)
{
	std::vector<CovariantTerm<R>> terms = a.terms();
	for (auto &t : b.terms())
		terms.push_back({-t.coefficient, t.atoms});
	return Covariant<R>(std::move(terms));
}
template <typename R>
Covariant<R> inline operator*(Covariant<R> const &a, Covariant<R> const &b)
{
	std::vector<CovariantTerm<R>> terms;
	terms.reserve(a.terms().size() * b.terms().size());
	for (auto &t : a.terms())
		for (auto &s : b.terms())
		{
			CovariantTerm<R> term;
			term.coefficient = t.coefficient * s.coefficient;

			// append the two atom lists while renaming all inner indices to
			// avoid collisions. Will be renamed again in the constructor
			term.atoms.reserve(t.atoms.size() + s.atoms.size());
			for (auto &atom : t.atoms)
			{
				term.atoms.push_back(atom);
				for (int &k : term.atoms.back().indices)
					if (k >= a.first_inner_index())
						k += 1000000;
			}
			for (auto &atom : s.atoms)
			{
				term.atoms.push_back(atom);
				for (int &k : term.atoms.back().indices)
					if (k >= b.first_inner_index())
						k += 2000000;
			}

			terms.push_back(std::move(term));
		}
	return Covariant<R>(std::move(terms));
}

template <typename R>
inline Covariant<R> operator*(Covariant<R> const &a, int b)
{
	std::vector<CovariantTerm<R>> terms = a.terms();
	for (auto &t : terms)
		t.coefficient *= b;
	return Covariant<R>(std::move(terms));
}
template <typename R>
inline Covariant<R> operator/(Covariant<R> const &a, int b)
{
	std::vector<CovariantTerm<R>> terms = a.terms();
	for (auto &t : terms)
		t.coefficient /= b;
	return Covariant<R>(std::move(terms));
}

} // namespace chalk

template <typename R> struct fmt::formatter<chalk::Covariant<R>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Covariant<R> &cov, FormatContext &ctx)
	    -> decltype(ctx.out())
	{
		auto const &terms = cov.terms();

		// no terms -> output "0"
		if (terms.empty())
			return format_to(ctx.out(), "0");

		// otherwise list the terms with " + " inbetween
		auto it = ctx.out();
		for (size_t i = 0; i < terms.size(); ++i)
		{
			if (i != 0)
				it = format_to(it, " + ");

			it = format_to(it, "{}", terms[i].coefficient);
			for (auto &atom : terms[i].atoms)
			{
				it = format_to(it, "*{}", atom.symbol);
				if (!atom.indices.empty())
					it++ = '_';
				for (int k : atom.indices)
					it = format_to(it, "{}", k);
			}
		}
		return it;
	}
};

#endif
