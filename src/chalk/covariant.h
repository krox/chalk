#ifndef CHALK_COVARIANT
#define CHALK_COVARIANT

/**
 * A 'Covariant' is a linear combination of indexed expressions, i.e. a
 * monoid ring over indexed expressions. This module provides
 *   - wrappers around operations on 'Indexed' to work with whole expressions
 *   - simplifications based on algebra of Lie-derivatives
 */

#include "absl/container/flat_hash_map.h"
#include "absl/container/inlined_vector.h"
#include "chalk/group_ring.h"
#include "chalk/indexed.h"
#include "fmt/format.h"
#include <algorithm>
#include <cassert>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace chalk {

template <typename R> using Covariant = Multinomial<R, Indexed>;
template <typename R> using CovariantTerm = MultinomialTerm<R, Indexed>;

template <typename R>
Covariant<R> inner_product(Covariant<R> const &a, Covariant<R> const &b)
{
	std::vector<CovariantTerm<R>> terms;
	terms.reserve(a.terms().size() * b.terms().size());
	for (auto &t : a.terms())
		for (auto &s : b.terms())
			terms.push_back({t.coefficient * s.coefficient,
			                 inner_product(t.index, s.index)});
	return Covariant<R>(std::move(terms));
}

template <typename R> Covariant<R> contract(Covariant<R> const &a, int i, int j)
{
	auto terms = a.terms();
	for (auto &term : terms)
		term.index = contract(term.index, i, j);
	return Covariant<R>(std::move(terms));
}

/**
 * derivative adding new index to the right
 * lower-case symbols are assumed to be constants
 * upper-case symbols simply obtain new indices
 */
template <typename R> inline Covariant<R> diff(Covariant<R> const &a)
{
	std::vector<CovariantTerm<R>> terms;
	for (auto &t : a.terms())
		for (size_t i = 0; i < t.index.size(); ++i)
		{
			if (!std::isupper(t.index[i].symbol[0]))
				continue;
			terms.push_back(t);
			std::vector<IndexedAtom> tmp = terms.back().index.release();
			tmp[i].indices.push_back(10000);
			terms.back().index = Indexed(std::move(tmp));
		}
	return Covariant<R>(std::move(terms));
}

/**
 * derivative adding a "'" to symbols
 * lower-case symbols are assumed to be constants
 * upper-case symbols optain new prime.
 */
template <typename R> inline Covariant<R> diff_1d(Covariant<R> const &a)
{
	std::vector<CovariantTerm<R>> terms;
	for (auto &t : a.terms())
		for (size_t i = 0; i < t.index.size(); ++i)
		{
			if (!std::isupper(t.index[i].symbol[0]))
				continue;
			terms.push_back(t);
			std::vector<IndexedAtom> tmp = terms.back().index.release();
			tmp[i].symbol = tmp[i].symbol + "'";
			terms.back().index = Indexed(std::move(tmp));
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
	for (auto &t : force.terms())
		assert(t.index.rank() == 1);

	Covariant<R> result = Covariant<R>(0);
	for (int d = 0; d <= degree; ++d)
	{
		Covariant<R> term = S;
		R prefactor = R(1);
		for (int i = 0; i < d; ++i)
			term = diff(term);
		for (int i = 0; i < d; ++i)
		{
			term = inner_product(term, force);
			prefactor /= (i + 1);
		}
		result += term * prefactor;
	}
	return result;
}

/**
 * build a taylor like expression  S(f) = S + S'f + 1/2 S''ff + ...
 *   - f may not have any open index
 *   - S can have arbitrary open indices (which simply stay that way)
 */
template <typename R>
inline Covariant<R> taylor_1d(Covariant<R> const &S, Covariant<R> const &force,
                              int degree)
{
	// the force should have exactly one open index
	for (auto &t : force.terms())
		assert(t.index.rank() == 0);

	Covariant<R> result = Covariant<R>(0);
	for (int d = 0; d <= degree; ++d)
	{
		Covariant<R> term = S;
		R prefactor = R(1);
		for (int i = 0; i < d; ++i)
			term = diff_1d(term);
		for (int i = 0; i < d; ++i)
		{
			term = term * force;
			prefactor /= (i + 1);
		}
		result += term * prefactor;
	}
	return result;
}

template <typename R>
Covariant<R> wick_contract(Covariant<R> const &a, std::string const &eta,
                           R const &variance)
{
	std::vector<CovariantTerm<R>> result;
	for (auto const &term : a.terms())
	{
		auto [new_terms, count] = wick_contract(term.index, eta);
		if (count % 2 != 0)
			continue;
		for (auto &new_term : new_terms)
		{
			result.push_back({term.coefficient, std::move(new_term)});
			for (int i = 0; i < count / 2; ++i)
				result.back().coefficient *= variance;
		}
	}
	return Covariant<R>(std::move(result));
}

// private building-blocks for Lie-simplification
namespace {

/** move double-indices to the front */
inline void move_double_indices(absl::InlinedVector<int, 4> &inds)
{
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

/** sort indices, returns signum (1,-1,0) of permutation */
inline int sort_antisym(absl::InlinedVector<int, 4> &inds)
{
	int r = 1;
	for (int k = (int)inds.size(); k > 1; --k)
		for (int i = 0; i < k - 1; ++i)
			if (inds[i] > inds[i + 1])
			{
				std::swap(inds[i], inds[i + 1]);
				r = -r;
			}
	for (int i = 0; i < (int)inds.size() - 1; ++i)
		if (inds[i] == inds[i + 1])
			return 0;
	return r;
}

/** returns true if a<->b can be proved to be symmetric, excluding symbol 'f' */
inline bool is_symmetric(std::vector<IndexedAtom> const &atoms, int a, int b)
{
	std::vector<IndexedAtom> atoms2 = atoms;
	int count_a = 0, count_b = 0; // for sanity check
	for (auto &atom : atoms2)
		if (atom.symbol != "f")
			for (auto &k : atom.indices)
			{
				if (k == a)
				{
					count_a++;
					k = b;
				}
				else if (k == b)
				{
					count_b++;
					k = a;
				}
			}
	assert(count_a <= 1);
	assert(count_b <= 1);
	if (count_a == 0 || count_b == 0)
		return false;
	std::sort(atoms2.begin(), atoms2.end());
	return atoms2 == atoms;
};

/**
 * modifies term inplace, returns new extra term if one was generated
 * (needs to be repeated in that case)
 */
template <typename R>
std::optional<CovariantTerm<R>>
sort_indices_brute_force(CovariantTerm<R> &term, std::string const &symbol)
{
	auto atoms = term.index.release();

	/** obsolete rules:
	 * 'sandwich rules':
	     - S(aba) -> S(baa) + cA/2 S(b)
	     - S(abcab) -> S(aabcb) + cA S(bcb)
	 * direct ordering rules:
	     - S(..ab..)S(..ba..) -> S(..ab..)S(..ab..) - cA/2 S(..c..) S(..c..)
	 These are now handled by explicitly ordering indices of S (introducing
	 structure constants f), and than contracting f with f and f with S
	 which produces some cA's.
	 (The explicit index-ordering only makes sense because the index-naming
	 quite strongly follows pure graph-structure.)
	 */

	// sort indices by adding explicit structure constants
	for (size_t ai = 0; ai < atoms.size(); ++ai)
	{
		if (atoms[ai].symbol != symbol)
			continue;
		auto &inds = atoms[ai].indices;

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
			std::vector<IndexedAtom> new_atoms = atoms;
			new_atoms[ai].indices.erase(new_atoms[ai].indices.begin() + i + 1);
			new_atoms[ai].indices[i] = 500;
			new_atoms.push_back({"f", {inds[i + 1], inds[i], 500}});
			term.index = Indexed(std::move(atoms));
			return CovariantTerm<R>{term.coefficient,
			                        Indexed(std::move(new_atoms))};
		}
	}

	// nothing found -> restore original
	term.index = Indexed(std::move(atoms));
	return std::nullopt;
}

} // namespace

/**
 * Simplifications that dont create new terms:
 *  - move double-indices to the front
 *  - sort indices of antisymmetric f
 *  - contract combinations of f using cA
 */
template <typename R>
void simplify_lie_term_basic(CovariantTerm<R> &term, std::string const &symbol,
                             R const &cA)
{
	auto atoms = term.index.release();
	for (auto &atom : atoms)
	{
		if (atom.symbol == symbol)
			move_double_indices(atom.indices);
		if (atom.symbol == "f")
		{
			// sort indices by anti-symmetry and look for sym/anti-sym -> 0
			assert(atom.indices.size() == 3);
			term.coefficient *= sort_antisym(atom.indices);

			if (is_symmetric(atoms, atom.indices[0], atom.indices[1]))
				term.coefficient = R(0);
			else if (is_symmetric(atoms, atom.indices[0], atom.indices[2]))
				term.coefficient = R(0);
			else if (is_symmetric(atoms, atom.indices[1], atom.indices[2]))
				term.coefficient = R(0);
		}
	}

again:

	for (size_t i = 0; i < atoms.size(); ++i)
		for (size_t j = 0; j < atoms.size(); ++j)
			if (i != j && atoms[i].symbol == "f" && atoms[j].symbol == "f")
			{
				// f_xab f_yab -> C_A delta_xy
				if (atoms[i].indices[1] == atoms[j].indices[1] &&
				    atoms[i].indices[2] == atoms[j].indices[2])
				{
					term.coefficient *= cA;
					atoms[i] = IndexedAtom{
					    "delta", {atoms[i].indices[0], atoms[j].indices[0]}};
					atoms.erase(atoms.begin() + j);
					goto again;
				}

				// f_abx f_yab -> C_A delta_xy
				if (atoms[i].indices[0] == atoms[j].indices[1] &&
				    atoms[i].indices[1] == atoms[j].indices[2])
				{
					term.coefficient *= cA;
					atoms[i] = IndexedAtom{
					    "delta", {atoms[i].indices[2], atoms[j].indices[0]}};
					atoms.erase(atoms.begin() + j);
					goto again;
				}
			}

	// f_xab f_yac f_zbc = 1/2 cA f_xyz
	for (size_t ai1 = 0; ai1 < atoms.size(); ++ai1)
		for (size_t ai2 = ai1 + 1; ai2 < atoms.size(); ++ai2)
			for (size_t ai3 = ai2 + 1; ai3 < atoms.size(); ++ai3)
				if (atoms[ai1].symbol == "f" && atoms[ai2].symbol == "f" &&
				    atoms[ai3].symbol == "f")
					if (atoms[ai1].indices[1] == atoms[ai2].indices[1] &&
					    atoms[ai1].indices[2] == atoms[ai3].indices[1] &&
					    atoms[ai2].indices[2] == atoms[ai3].indices[2])
					{
						term.coefficient *= cA / 2;
						atoms[ai1].indices[1] = atoms[ai2].indices[0];
						atoms[ai1].indices[2] = atoms[ai3].indices[0];
						atoms.erase(atoms.begin() + ai3);
						atoms.erase(atoms.begin() + ai2);
						goto again;
					}

	// f_abc S_ab = 1/2 C_A S_c
	for (size_t ai = 0; ai < atoms.size(); ++ai)
		for (size_t ai2 = 0; ai2 < atoms.size(); ++ai2)
			if (atoms[ai].symbol == symbol && atoms[ai2].symbol == "f")
			{
				for (int i = 0; i < (int)atoms[ai].indices.size() - 1; ++i)
				{
					if (atoms[ai].indices[i] == atoms[ai2].indices[0] &&
					    atoms[ai].indices[i + 1] == atoms[ai2].indices[1])
					{
						term.coefficient *= cA / 2;
						atoms[ai].indices.erase(atoms[ai].indices.begin() + i +
						                        1);
						atoms[ai].indices[i] = atoms[ai2].indices[2];
						atoms.erase(atoms.begin() + ai2);
						goto again;
					}
				}
			}
	// f_cab S_ab = 1/2 C_A S_c
	for (size_t ai = 0; ai < atoms.size(); ++ai)
		for (size_t ai2 = 0; ai2 < atoms.size(); ++ai2)
			if (atoms[ai].symbol == symbol && atoms[ai2].symbol == "f")
			{
				for (int i = 0; i < (int)atoms[ai].indices.size() - 1; ++i)
				{
					if (atoms[ai].indices[i] == atoms[ai2].indices[1] &&
					    atoms[ai].indices[i + 1] == atoms[ai2].indices[2])
					{
						term.coefficient *= cA / 2;
						atoms[ai].indices.erase(atoms[ai].indices.begin() + i +
						                        1);
						atoms[ai].indices[i] = atoms[ai2].indices[0];
						atoms.erase(atoms.begin() + ai2);
						goto again;
					}
				}
			}
	term.index = Indexed(std::move(atoms));
}

/** reorder lie-derivative indices using relations with the casimir operator */
template <typename R>
Covariant<R> simplify_lie(Covariant<R> const &a, std::string const &symbol,
                          R const &cA)
{
	(void)symbol;
	(void)cA;
	std::vector<CovariantTerm<R>> terms = a.terms();

	for (int iter = 0; iter < 5; ++iter)
	{
		// Step 1) Basic simplifications that dont create new terms

		for (auto &term : terms)
		{
			simplify_lie_term_basic(term, symbol, cA);
		}

		// more serious simplifications that make one term into multiple
		for (int term_i = 0; term_i < (int)terms.size(); ++term_i)
		{
			auto &term = terms[term_i];

			// if something was found -> new term and repeat on this term
			if (auto new_term = sort_indices_brute_force(term, symbol);
			    new_term)
			{
				// NOTE: the push_back can make 'term' a dangling reference
				terms.push_back(new_term.value());
				--term_i;
				continue;
			}
		}
	}

	return Covariant<R>(std::move(terms));
}

template <typename R>
Covariant<R> simplify_lie_commutative(Covariant<R> const &a,
                                      std::string const &symbol)
{
	std::vector<CovariantTerm<R>> terms = a.terms();
	for (auto &term : terms)
	{
		auto atoms = term.index.release();
		for (auto &atom : atoms)
		{
			if (atom.symbol == symbol)
				std::sort(atom.indices.begin(), atom.indices.end());
			if (atom.symbol == "f")
				term.coefficient = R(0);
		}
		term.index = Indexed(std::move(atoms));
	}
	return Covariant<R>(std::move(terms));
}

} // namespace chalk

#endif
