#ifndef CHALK_COVARIANT
#define CHALK_COVARIANT

#include "absl/container/flat_hash_map.h"
#include "absl/container/inlined_vector.h"
#include "fmt/format.h"
#include <algorithm>
#include <cassert>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace chalk {

/** a (named) object with indices */
struct CovariantAtom
{
	std::string symbol;
	absl::InlinedVector<int, 4> indices; // first internal, then derivatives

	bool operator==(CovariantAtom const &other) const
	{
		return symbol == other.symbol && indices == other.indices;
	}

	bool operator<(CovariantAtom const &other) const
	{
		if (symbol != other.symbol)
			return symbol < other.symbol;
		if (indices.size() != other.indices.size())
			return indices.size() > other.indices.size();
		return indices < other.indices;
	}
};

/**
 * Sort atoms and rename indices.
 * returns list of open indices, for convenience
 */
std::vector<int> normalize_indices(std::vector<CovariantAtom> &atoms);

/** contract all deltas involving inner indices */
bool simplify_delta(std::vector<CovariantAtom> &atoms);

/** moves double-indices to the front */
void move_double_indices(absl::InlinedVector<int, 4> &inds);

/** simplify antisymmetric symbol by sorting its indices. returns sign */
int simplify_antisymmetric(std::vector<CovariantAtom> &atoms,
                           std::string const &symbol);

/** product of atoms and a pre-factor */
template <typename R> struct CovariantTerm
{
	R coefficient;
	std::vector<CovariantAtom> atoms;

	bool operator<(CovariantTerm const &other) const
	{
		if (atoms.size() != other.atoms.size())
			return atoms.size() < other.atoms.size();
		return atoms < other.atoms;
	}
};

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
	std::vector<int> open_indices_;

	void cleanup()
	{
		if (terms_.size() == 0)
			return;

		// prune zero-elements
		// (has to be done early because 0 is allowed to have wrong indices)
		for (size_t i = 0; i < terms_.size(); ++i)
			if (terms_[i].coefficient == 0)
			{
				std::swap(terms_[i], terms_.back());
				terms_.pop_back();
				--i;
			}
		if (terms_.size() == 0)
			return;

		// cleanup individual terms
		open_indices_ = normalize_indices(terms_[0].atoms);
		first_inner_index_ =
		    open_indices_.empty() ? 1 : open_indices_.back() + 1;
		for (size_t i = 1; i < terms_.size(); ++i)
			if (open_indices_ != normalize_indices(terms_[i].atoms))
				throw std::runtime_error(
				    fmt::format("inconsistent indices in '{}'", *this));
		for (auto &term : terms_)
			if (simplify_delta(term.atoms))
				normalize_indices(term.atoms);

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
	std::vector<int> const &open_indices() const { return open_indices_; }
};

/** unary operations */
template <typename R> Covariant<R> inline operator-(Covariant<R> const &a)
{
	auto terms = a.terms();
	for (auto &t : terms)
		t.coefficient = -t.coefficient;
	return Covariant<R>(std::move(terms));
}

/** binary operations (Covariant <-> Covariant) */
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

/** binary operations (Covariant <-> R) */
template <typename R>
inline Covariant<R> operator*(Covariant<R> const &a, R const &b)
{
	std::vector<CovariantTerm<R>> terms = a.terms();
	for (auto &t : terms)
		t.coefficient *= b;
	return Covariant<R>(std::move(terms));
}
template <typename R>
inline Covariant<R> operator/(Covariant<R> const &a, R const &b)
{
	std::vector<CovariantTerm<R>> terms = a.terms();
	for (auto &t : terms)
		t.coefficient /= b;
	return Covariant<R>(std::move(terms));
}

/** binary operations (Covariant <-> int) */
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

/** convenience operators */
template <typename R>
inline void operator+=(Covariant<R> &a, Covariant<R> const &b)
{
	a = a + b;
}
template <typename R>
inline void operator-=(Covariant<R> &a, Covariant<R> const &b)
{
	a = a - b;
}
template <typename R>
inline void operator*=(Covariant<R> &a, Covariant<R> const &b)
{
	a = a * b;
}
template <typename R> inline void operator+=(Covariant<R> &a, R const &b)
{
	a = a + b;
}
template <typename R> inline void operator-=(Covariant<R> &a, R const &b)
{
	a = a - b;
}
template <typename R> inline void operator*=(Covariant<R> &a, R const &b)
{
	a = a * b;
}
template <typename R> inline void operator+=(Covariant<R> &a, int b)
{
	a = a + b;
}
template <typename R> inline void operator-=(Covariant<R> &a, int b)
{
	a = a - b;
}
template <typename R> inline void operator*=(Covariant<R> &a, int b)
{
	a = a * b;
}
template <typename R> inline void operator/=(Covariant<R> &a, int b)
{
	a = a / b;
}

/**
 * replace open index 'from' to 'to'
 */
template <typename R>
inline Covariant<R> rename_index(Covariant<R> const &a, int from, int to)
{
	if (from == to)
		return a;
	auto terms = a.terms();
	for (auto &t : terms)
	{
		int found = 0;
		for (auto &atom : t.atoms)
			for (int &k : atom.indices)
			{
				// get old inner index 'to' out of the way
				if (k == to)
					k = 1000000;

				// rename 'from' to 'to'
				if (k == from)
				{
					k = to;
					++found;
				}
			}
		if (found != 1)
			throw std::runtime_error(
			    "tried renaming an index that wasn't there");
	}
	return Covariant<R>(std::move(terms));
}

/**
 * apply a function to all coefficients
 */
template <typename R, typename F>
inline Covariant<R> map_coefficients(Covariant<R> const &a, F &&f)
{
	auto terms = a.terms();
	for (auto &t : terms)
		t.coefficient = f(t.coefficient);
	return Covariant<R>(std::move(terms));
}

/**
 * differentiate with index k
 * lower-case symbols are assumed to be constants
 * upper-case symbols simply obtain new indices
 */
template <typename R> inline Covariant<R> diff(Covariant<R> const &a, int k)
{
	std::vector<CovariantTerm<R>> terms;
	for (auto &t : a.terms())
		for (size_t i = 0; i < t.atoms.size(); ++i)
		{
			if (!std::isupper(t.atoms[i].symbol[0]))
				continue;
			terms.push_back(t);
			terms.back().atoms[i].indices.push_back(k);
		}
	return Covariant<R>(std::move(terms));
}

/**
 * build a taylor like expression  S(f) = S + S'f + 1/2 S''ff + ...
 *   - f needs exactly one open index
 *   - S can have arbitrary open indices (which simply stay that way)
 */
template <typename R>
inline Covariant<R> taylor(Covariant<R> const &S, Covariant<R> const &force,
                           int degree)
{
	// the force should have exactly one open index
	assert(force.open_indices().size() == 1);
	int k = force.open_indices()[0];

	fmt::print("open index = {}\n", k);

	Covariant<R> result = Covariant<R>(0);
	for (int d = 0; d <= degree; ++d)
	{
		Covariant<R> term = S;
		R prefactor = R(1);
		for (int i = 0; i < d; ++i)
			term = diff(term, i - 100);
		for (int i = 0; i < d; ++i)
		{
			term *= rename_index(force, k, i - 100);
			prefactor /= (i + 1);
		}
		result += term * prefactor;
	}
	return result;
}

std::vector<std::vector<std::pair<int, int>>> wick_list(int n);

template <typename R>
inline Covariant<R> wick_contract(Covariant<R> const &a, std::string const &eta,
                                  R const &variance)
{
	std::vector<CovariantTerm<R>> result;
	for (auto const &term : a.terms())
	{
		// split term into eta / non-eta
		std::vector<int> indices;
		std::vector<CovariantAtom> other;
		for (auto const &atom : term.atoms)
		{
			if (atom.symbol == eta)
			{
				assert(atom.indices.size() == 1);
				indices.push_back(atom.indices[0]);
			}
			else
				other.push_back(atom);
		}
		// odd powers vanish
		if (indices.size() % 2 != 0)
			continue;

		// otherwise, remainder * wick-contraction terms
		for (auto &l : wick_list((int)indices.size()))
		{
			assert(l.size() == indices.size() / 2);
			result.push_back({term.coefficient, other});
			for (size_t i = 0; i < indices.size() / 2; ++i)
				result.back().coefficient *= variance;
			for (auto &p : l)
				result.back().atoms.push_back(
				    {"delta", {indices[p.first], indices[p.second]}});
		}
	}
	return Covariant<R>(std::move(result));
}

// "private" building-blocks for Lie-simplification
namespace {

/**
 * modifies term inplace, returns new extra term if one was generated
 * (needs to be repeated in that case)
 */
template <typename R>
std::optional<CovariantTerm<R>>
simplify_lie_term_sandwich(CovariantTerm<R> &term, std::string const &symbol,
                           R const &cA)
{
	for (size_t atom_i = 0; atom_i < term.atoms.size(); ++atom_i)
	{
		auto &atom = term.atoms[atom_i];
		if (atom.symbol != symbol)
			continue;
		auto &inds = atom.indices;

		// aba -> baa + cA/2 b
		for (int i = 0; i < (int)inds.size() - 2; ++i)
			if (inds[i] == inds[i + 2])
			{
				std::swap(inds[i], inds[i + 1]);

				CovariantTerm<R> new_term = term;
				auto &new_inds = new_term.atoms[atom_i].indices;
				new_inds.erase(new_inds.begin() + i + 2);
				new_inds.erase(new_inds.begin() + i + 1);
				new_term.coefficient *= cA / 2;
				return new_term;
			}

		// abcab -> aabcb + cA bcb
		for (int i = 0; i < (int)inds.size() - 4; ++i)
			if (inds[i] == inds[i + 3] && inds[i + 1] == inds[i + 4])
			{
				inds[i + 3] = inds[i + 2];
				inds[i + 2] = inds[i + 4];
				inds[i + 1] = inds[i];
				assert(inds[i] == inds[i + 1] && inds[i + 2] == inds[i + 4]);

				CovariantTerm<R> new_term = term;
				auto &new_inds = new_term.atoms[atom_i].indices;
				new_inds.erase(new_inds.begin() + i + 1);
				new_inds.erase(new_inds.begin() + i);
				new_term.coefficient *= cA;
				return new_term;
			}
	}
	return std::nullopt;
}

template <typename R>
std::optional<CovariantTerm<R>>
simplify_lie_term_order(CovariantTerm<R> &term, std::string const &symbol,
                        R const &cA)
{
	// S(..ab..)S(..ba..) = S(..ab..)S(..ab..) - cA/2 S(..c..) S(..c..)
	for (size_t ai = 0; ai < term.atoms.size(); ++ai)
		for (size_t ai2 = ai + 1; ai2 < term.atoms.size(); ++ai2)
		{
			if (term.atoms[ai].symbol != symbol ||
			    term.atoms[ai2].symbol != symbol)
				continue;

			auto &inds = term.atoms[ai].indices;
			auto &inds2 = term.atoms[ai2].indices;

			for (int i = 0; i < (int)inds.size() - 1; ++i)
				for (int i2 = 0; i2 < (int)inds2.size() - 1; ++i2)
					if (inds[i] == inds2[i2 + 1] && inds[i + 1] == inds2[i2])
					{
						std::swap(inds2[i2], inds2[i2 + 1]);
						CovariantTerm<R> new_term = term;
						new_term.atoms[ai].indices.erase(
						    new_term.atoms[ai].indices.begin() + i);
						new_term.atoms[ai2].indices.erase(
						    new_term.atoms[ai2].indices.begin() + i2);
						new_term.coefficient *= -cA / 2;
						return new_term;
					}
		}
	return std::nullopt;
}

/** commute the indices simplfy by explicitly using structure constants */
template <typename R>
std::optional<CovariantTerm<R>>
simplify_lie_term_order_full(CovariantTerm<R> &term, std::string const &symbol)
{
	for (size_t ai = 0; ai < term.atoms.size(); ++ai)
	{
		if (term.atoms[ai].symbol != symbol)
			continue;
		auto &inds = term.atoms[ai].indices;

		for (int i = 0; i < (int)inds.size() - 1; ++i)
		{
			// good order -> nothing to do
			if (inds[i] <= inds[i + 1])
				continue;

			// one index is part of a pair -> leave it alone
			if (i >= 1 && inds[i - 1] == inds[i])
				continue;
			if (i < (int)inds.size() - 2 && inds[i + 1] == inds[i + 2])
				continue;

			// otherwise -> transpose indices by adding extra term involving an
			//              explicit structure constant f_{ijk}
			std::swap(inds[i], inds[i + 1]);
			CovariantTerm<R> new_term = term;
			new_term.atoms[ai].indices.erase(
			    new_term.atoms[ai].indices.begin() + i + 1);
			new_term.atoms[ai].indices[i] = 500;
			new_term.atoms.push_back({"f", {inds[i + 1], inds[i], 500}});
			return new_term;
		}
	}
	return std::nullopt;
}

template <typename R>
void simplify_lie_term_structure_constants(CovariantTerm<R> &term,
                                           std::string const &symbol,
                                           R const &cA)
{
// {xab}{ybc}{zca} = 1/2 cA {xyz}
// {xab}{yac}{zbc} = 1/2 cA {xyz} (this form appears in normalized terms)
foo:
	for (size_t i = 0; i < term.atoms.size(); ++i)
		for (size_t j = i + 1; j < term.atoms.size(); ++j)
			for (size_t k = j + 1; k < term.atoms.size(); ++k)
			{
				if (term.atoms[i].symbol != "f")
					continue;
				if (term.atoms[j].symbol != "f")
					continue;
				if (term.atoms[k].symbol != "f")
					continue;

				auto &a = term.atoms[i].indices;
				auto &b = term.atoms[j].indices;
				auto &c = term.atoms[k].indices;
				assert(a.size() == 3 && b.size() == 3 && c.size() == 3);
				if (a[1] == b[1] && a[2] == c[1] && b[2] == c[2])
				{
					a[1] = b[0];
					a[2] = c[0];
					term.coefficient *= cA / 2;

					term.atoms.erase(term.atoms.begin() + k);
					term.atoms.erase(term.atoms.begin() + j);
					goto foo;
				}
			}

// f_abc S_ab = 1/2 C_A S_c
again:
	for (auto &atom_f : term.atoms)
		for (auto &atom_s : term.atoms)
		{
			if (atom_s.symbol != symbol)
				continue;
			if (atom_f.symbol != "f")
				continue;
			for (int i = 0; i < (int)atom_s.indices.size() - 1; ++i)
			{
				if (atom_s.indices[i] == atom_f.indices[0] &&
				    atom_s.indices[i + 1] == atom_f.indices[1])
				{
					atom_s.indices.erase(atom_s.indices.begin() + i + 1);
					atom_s.indices[i] = atom_f.indices[2];
					term.coefficient *= cA / 2;
					std::swap(atom_f, term.atoms.back());
					term.atoms.pop_back(); // invalidated refs, so get out
					goto again;
				}
				if (atom_s.indices[i] == atom_f.indices[1] &&
				    atom_s.indices[i + 1] == atom_f.indices[2])
				{
					atom_s.indices.erase(atom_s.indices.begin() + i + 1);
					atom_s.indices[i] = atom_f.indices[0];
					term.coefficient *= cA / 2;
					std::swap(atom_f, term.atoms.back());
					term.atoms.pop_back(); // invalidated refs, so get out
					goto again;
				}
			}
		}
} // namespace

} // namespace

/** reorder lie-derivative indices using relations with the casimir operator */
template <typename R>
Covariant<R> simplify_lie(Covariant<R> const &a, std::string const &symbol,
                          R const &cA)
{
	std::vector<CovariantTerm<R>> terms = a.terms();

	// step 1: simple rules that do NOT depend on index naming
	for (size_t term_i = 0; term_i < terms.size(); ++term_i)
	{
		auto &term = terms[term_i];

		// move all double-indices to the front
		for (auto &atom : term.atoms)
			if (atom.symbol == symbol)
				move_double_indices(atom.indices);

		// simple rules (that do NOT depend on index naming)
		auto new_term = simplify_lie_term_sandwich(term, symbol, cA);
		if (!new_term)
			new_term = simplify_lie_term_order(term, symbol, cA);

		// if something was found -> new term and repeat on this term
		if (new_term)
		{
			// NOTE: this can make 'term' a dangling reference
			terms.push_back(new_term.value());
			--term_i;
			continue;
		}

		// brute-force rule (that does depends on index naming)
		normalize_indices(term.atoms);
		term.coefficient *= simplify_antisymmetric(term.atoms, "f");
		simplify_lie_term_structure_constants(term, symbol, cA);
		new_term = simplify_lie_term_order_full(term, symbol);

		if (new_term)
		{
			// NOTE: this can make 'term' a dangling reference
			terms.push_back(new_term.value());
			--term_i;
			continue;
		}
	}

	// step 2: rename indices into standard form
	for (auto &term : terms)
		normalize_indices(term.atoms);

	return Covariant<R>(std::move(terms));
}

/**
 *  Write to stdout using multiple lines.
 *  More readable than fmt::print("{}",a) for large expressions
 */
template <typename R> void dump(Covariant<R> const &a)
{
	for (auto &t : a.terms())
	{
		if (t.atoms.empty())
			fmt::print("1");
		for (auto &atom : t.atoms)
		{
			fmt::print(" {}", atom.symbol);
			if (!atom.indices.empty())
				fmt::print("_");
			for (int k : atom.indices)
				fmt::print("{}", k);
		}

		fmt::print(" : {}\n", t.coefficient);
	}
}
template <typename R> void dump_summary(Covariant<R> const &a)
{
	for (auto &t : a.terms())
	{
		if (t.atoms.empty())
			fmt::print("1");
		for (auto &atom : t.atoms)
		{
			fmt::print(" {}", atom.symbol);
			if (!atom.indices.empty())
				fmt::print("_");
			for (int k : atom.indices)
				fmt::print("{}", k);
		}

		fmt::print(" : ... polynomial with {} terms ...\n",
		           t.coefficient.terms().size());
	}
}

} // namespace chalk

template <> struct fmt::formatter<chalk::CovariantAtom>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::CovariantAtom &atom, FormatContext &ctx)
	    -> decltype(ctx.out())
	{
		auto it = ctx.out();
		it = format_to(it, "{}", atom.symbol);
		if (!atom.indices.empty())
			it++ = '_';
		for (int k : atom.indices)
			it = format_to(it, "{}", k);
		return it;
	}
};

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

			// TODO: skip the parens if the coeffs are already atomic
			// (also skip fully if coeff == 1)
			it = format_to(it, "({})", terms[i].coefficient);
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
