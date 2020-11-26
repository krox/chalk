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

template <typename R>
std::pair<R, absl::InlinedVector<int, 4>>
simplify_lie_derivatives(absl::InlinedVector<int, 4> &inds, R const &cA)
{
	// aba -> baa + cA/2 b
	for (int i = 0; i < (int)inds.size() - 2; ++i)
		if (inds[i] == inds[i + 2])
		{
			std::swap(inds[i], inds[i + 1]);

			auto new_inds = inds;
			new_inds.erase(new_inds.begin() + i + 2);
			new_inds.erase(new_inds.begin() + i + 1);
			return {cA / 2, new_inds};
		}

	// abcab -> aabcb + cA bcb
	for (int i = 0; i < (int)inds.size() - 4; ++i)
		if (inds[i] == inds[i + 3] && inds[i + 1] == inds[i + 4])
		{
			inds[i + 3] = inds[i + 2];
			inds[i + 2] = inds[i + 4];
			inds[i + 1] = inds[i];
			assert(inds[i] == inds[i + 1] && inds[i + 2] == inds[i + 4]);

			auto &new_inds = inds;
			new_inds.erase(new_inds.begin() + i + 1);
			new_inds.erase(new_inds.begin() + i);
			return {cA, new_inds};
		}

	// nothing found
	return {R(0), {}};
}

/**
 * modifies term inplace, returns new extra term if one was generated
 * (needs to be repeated in that case)
 */
template <typename R>
std::optional<CovariantTerm<R>>
simplify_lie_derivatives(CovariantTerm<R> &term, std::string const &symbol,
                         R const &cA)
{
	auto atoms = term.index.release();

	for (size_t atom_i = 0; atom_i < atoms.size(); ++atom_i)
	{
		if (atoms[atom_i].symbol != symbol)
			continue;

		if (auto [factor, new_inds] =
		        simplify_lie_derivatives(atoms[atom_i].indices, cA);
		    !(factor == 0))
		{
			std::vector<IndexedAtom> new_atoms = atoms;
			new_atoms[atom_i].indices = std::move(new_inds);
			term.index = Indexed(std::move(atoms));
			return CovariantTerm<R>{term.coefficient * factor,
			                        Indexed(std::move(new_atoms))};
		}
	}

	// S(..ab..)S(..ba..) = S(..ab..)S(..ab..) - cA/2 S(..c..) S(..c..)
	for (size_t ai = 0; ai < atoms.size(); ++ai)
		for (size_t ai2 = ai + 1; ai2 < atoms.size(); ++ai2)
		{
			if (atoms[ai].symbol != symbol || atoms[ai2].symbol != symbol)
				continue;

			auto &inds = atoms[ai].indices;
			auto &inds2 = atoms[ai2].indices;

			for (int i = 0; i < (int)inds.size() - 1; ++i)
				for (int i2 = 0; i2 < (int)inds2.size() - 1; ++i2)
					if (inds[i] == inds2[i2 + 1] && inds[i + 1] == inds2[i2])
					{
						std::swap(inds2.at(i2), inds2.at(i2 + 1));
						std::vector<IndexedAtom> new_atoms = atoms;
						new_atoms[ai].indices.erase(
						    new_atoms[ai].indices.begin() + i);
						new_atoms[ai2].indices.erase(
						    new_atoms[ai2].indices.begin() + i2);
						term.index = Indexed(std::move(atoms));
						return CovariantTerm<R>{term.coefficient * (-cA / 2),
						                        Indexed(std::move(new_atoms))};
					}
		}

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

	// sort indices of structure constants using anti-symmetry
	for (auto &atom : atoms)
		if (atom.symbol == "f")
		{
			assert(atom.indices.size() == 3);
			if (is_symmetric(atoms, atom.indices[0], atom.indices[1]))
				term.coefficient = R(0);
			else if (is_symmetric(atoms, atom.indices[0], atom.indices[2]))
				term.coefficient = R(0);
			else if (is_symmetric(atoms, atom.indices[1], atom.indices[2]))
				term.coefficient = R(0);
			else
				term.coefficient *= sort_antisym(atom.indices);
		}

foo:
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
						goto foo;
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
						goto foo;
					}
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
						goto foo;
					}

	// f_xab f_yab = C_A delta_xy
	for (size_t ai1 = 0; ai1 < atoms.size(); ++ai1)
		for (size_t ai2 = ai1 + 1; ai2 < atoms.size(); ++ai2)
			if (atoms[ai1].symbol == "f" && atoms[ai2].symbol == "f")
				if (atoms[ai1].indices[1] == atoms[ai2].indices[1] &&
				    atoms[ai1].indices[2] == atoms[ai2].indices[2])
				{
					term.coefficient *= cA;
					atoms[ai1] = IndexedAtom{
					    "delta",
					    {atoms[ai1].indices[0], atoms[ai2].indices[0]}};
					atoms.erase(atoms.begin() + ai2);
					goto foo;
				}
	// f_xab f_aby = C_A delta_xy
	for (size_t ai1 = 0; ai1 < atoms.size(); ++ai1)
		for (size_t ai2 = ai1 + 1; ai2 < atoms.size(); ++ai2)
			if (atoms[ai1].symbol == "f" && atoms[ai2].symbol == "f")
				if (atoms[ai1].indices[1] == atoms[ai2].indices[0] &&
				    atoms[ai1].indices[2] == atoms[ai2].indices[1])
				{
					term.coefficient *= cA;
					atoms[ai1] = IndexedAtom{
					    "delta",
					    {atoms[ai1].indices[0], atoms[ai2].indices[2]}};
					atoms.erase(atoms.begin() + ai2);
					goto foo;
				}

	// nothing found -> restore original
	term.index = Indexed(std::move(atoms));
	return std::nullopt;
}

} // namespace

/** reorder lie-derivative indices using relations with the casimir operator */
template <typename R>
Covariant<R> simplify_lie(Covariant<R> const &a, std::string const &symbol,
                          R const &cA)
{
	std::vector<CovariantTerm<R>> terms = a.terms();

	for (int iter = 0; iter < 5; ++iter)
		for (size_t term_i = 0; term_i < terms.size(); ++term_i)
		{
			auto &term = terms[term_i];

			// move all double-indices to the front
			{
				auto tmp = term.index.release();
				for (auto &atom : tmp)
					if (atom.symbol == symbol)
						move_double_indices(atom.indices);
				term.index = Indexed(std::move(tmp));
			}

			// if something was found -> new term and repeat on this term
			if (auto new_term = simplify_lie_derivatives(term, symbol, cA);
			    new_term)
			{
				// NOTE: the push_back can make 'term' a dangling reference
				terms.push_back(new_term.value());
				--term_i;
				continue;
			}
		}

	return Covariant<R>(std::move(terms));
}

} // namespace chalk

#endif
