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
	if (a.indices.size() != b.indices.size())
		return a.indices.size() > b.indices.size();

	// otherwise arbitrary
	return a.indices < b.indices;
}

/** product of atoms and a pre-factor */
template <typename R> struct CovariantTerm
{
	R coefficient;
	std::vector<CovariantAtom> atoms;

	void rename_index(int from, int to)
	{
		assert(from != to);
		int count = 0;
		for (auto &atom : atoms)
			for (int &k : atom.indices)
				if (k == from)
				{
					k = to;
					count += 1;
				}
				else if (k == to)
					count += 1;
		assert(count == 1 || count == 2);
	}

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
		absl::flat_hash_map<int, int> count;
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
				throw std::runtime_error(
				    "inconsistent indices in expression (index-count != 1,2)");
		}
		std::sort(open.begin(), open.end());

		// handle delta
		for (size_t i = 0; i < atoms.size(); ++i)
		{
			auto &a = atoms[i];
			if (a.symbol != "delta")
				continue;

			if (a.indices.size() != 2)
				throw std::runtime_error("delta needs exactly two indices");
			if (a.indices[0] == a.indices[1])
				throw std::runtime_error(
				    "'delta_ii' not implemented (unknown dimension)");
			auto ind1 = a.indices[0];
			auto ind2 = a.indices[1];

			if (count[ind1] == 1)
				std::swap(ind1, ind2);
			if (count[ind2] != 2)
				continue;

			std::swap(a, atoms.back());
			atoms.pop_back();
			--i;

			rename_index(ind1, ind2);
		}

		// rename inner indices
		{
			int z = open.empty() ? 1 : open.back() + 1; // first inner index

			absl::flat_hash_map<int, int> trans;
			for (auto &a : atoms)
				for (int &k : a.indices)
				{
					if (count[k] != 2)
						continue;
					if (trans.count(k) == 0)
						trans[k] = z++;
					k = trans[k];
				}
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
		open_indices_ = terms_[0].cleanup();
		first_inner_index_ =
		    open_indices_.empty() ? 1 : open_indices_.back() + 1;
		for (size_t i = 1; i < terms_.size(); ++i)
			if (open_indices_ != terms_[i].cleanup())
				throw std::runtime_error(
				    fmt::format("inconsistent indices in '{}'", *this));

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

	Covariant<R> result = Covariant<R>(0);
	for (int d = 0; d <= degree; ++d)
	{
		Covariant<R> term = S;
		R prefactor = R(1);
		for (int i = 0; i < d; ++i)
			term = diff(term, i + 100);
		for (int i = 0; i < d; ++i)
		{
			term *= rename_index(force, k, i + 100);
			prefactor /= (i + 1);
		}
		result += term * prefactor;
	}
	return result;
}

std::vector<std::vector<std::pair<int, int>>> wick_list(int n)
{
	assert(n >= 0);
	std::vector<std::vector<std::pair<int, int>>> result;
	if (n % 2 != 0)
		return result;

	if (n == 0)
	{
		result.push_back({});
		return result;
	}

	for (auto &partial : wick_list(n - 2))
	{
		result.push_back(partial);
		result.back().push_back({n - 2, n - 1});
		for (int k = 0; k < n / 2 - 1; ++k)
		{
			result.push_back(partial);
			result.back().push_back({n - 2, n - 1});
			std::swap(result.back().back().first, result.back()[k].first);
			result.push_back(partial);
			result.back().push_back({n - 2, n - 1});
			std::swap(result.back().back().first, result.back()[k].second);
		}
	}
	return result;
}

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

/** reorder indices assuming complete symmetry */
template <typename R>
Covariant<R> simplify_commutative_indices(Covariant<R> const &a,
                                          std::string const &symbol)
{
	std::vector<CovariantTerm<R>> terms = a.terms();
	for (auto &term : terms)
		for (auto &atom : term.atoms)
			if (atom.symbol == symbol)
				std::sort(atom.indices.begin(), atom.indices.end());
	for (size_t term_i = 0; term_i < terms.size(); ++term_i)
	{
		auto &term = terms[term_i];
		for (size_t atom_i = 0; atom_i < term.atoms.size(); ++atom_i)
		{
			auto &atom = term.atoms[atom_i];
			if (atom.symbol != symbol)
				continue;
			auto &inds = atom.indices;

			// double-indices commute completey -> move them to the front
			for (int i = 1; i < (int)inds.size() - 1; ++i)
				// is double index?
				if (inds[i] == inds[i + 1])
					// preceding is not double index ?
					if (i == 1 || inds[i - 1] != inds[i - 2])
					{
						std::swap(inds[i + 1], inds[i - 1]);
					}
		}
	}
	return Covariant<R>(std::move(terms));
}

// "private" building-blocks for Lie-simplification
namespace {

/** moves double-indices to the front */
inline void move_double_indices(CovariantAtom &atom)
{
	auto &inds = atom.indices;
	while (true)
	{
		bool change = false;
		for (int i = 1; i < (int)inds.size() - 1; ++i)
			// is double index?
			if (inds[i] == inds[i + 1])
				// preceding is not double index ?
				if (i == 1 || inds[i - 1] != inds[i - 2])
				{
					std::swap(inds[i + 1], inds[i - 1]);
					change = true;
				}
		if (!change)
			break;
	}
}

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
} // namespace

/** reorder lie-derivative indices using relations with the casimir operator */
template <typename R>
Covariant<R> simplify_lie(Covariant<R> const &a, std::string const &symbol,
                          R const &cA)
{
	std::vector<CovariantTerm<R>> terms = a.terms();
	for (size_t term_i = 0; term_i < terms.size(); ++term_i)
	{
		while (true)
		{
			auto &term = terms[term_i];
			for (auto &atom : term.atoms)
				if (atom.symbol == symbol)
					move_double_indices(atom);

			auto new_term = simplify_lie_term_sandwich(term, symbol, cA);
			if (!new_term)
				new_term = simplify_lie_term_order(term, symbol, cA);

			if (new_term)
			{
				// NOTE: this potentially makes 'term' a dangling reference!
				terms.push_back(std::move(new_term.value()));
			}
			else
				break;
		}
	}

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
