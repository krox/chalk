#ifndef CHALK_INDEXED_H
#define CHALK_INDEXED_H

#include "fmt/format.h"
#include "util/vector.h"
#include <algorithm>
#include <cassert>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace chalk {

/**
 * A named object with indices
 *   - TODO: symbols should have some flags (constant/(anti-)symmetric)
 *           (right now its hardcoded to the symbol name)
 */
struct IndexedAtom
{
	std::string symbol;
	util::small_vector<int, 4> indices = {};

	// temporary helper during sort_atoms()
	int64_t prio_ = {};
};

/** arbitrary order ( not the one used in Indexed::cleanup() ) */
bool operator==(IndexedAtom const &a, IndexedAtom const &b);
bool operator<(IndexedAtom const &a, IndexedAtom const &b);

/**
 * Product of indexed symbols.
 *   - repeated indices are summed over
 *   - multiplication '*' is outer product, so that this is a monoid
 *   - not quite commutative, due to numbering of open indices
 *   - does not distinguish upper/lower indices
 *   - does not know what the indices are for (tensor, derivatives)
 */
class Indexed
{
	std::vector<IndexedAtom> atoms_;
	int rank_ = 0; // number of open indices, which are 0,...,n-1

	// orders atoms and renames inner indices (also sets rank_)
	void cleanup();

  public:
	/** constructors */
	Indexed() = default;
	explicit Indexed(int val)
	{
		if (val != 1)
			throw std::runtime_error("invalid use of Indexed(int)");
	}
	explicit Indexed(std::vector<IndexedAtom> atoms) : atoms_(atoms)
	{
		cleanup();
	}
	explicit Indexed(std::string const &symbol,
	                 util::small_vector<int, 4> indices = {})
	    : atoms_{IndexedAtom{symbol, std::move(indices)}}
	{
		cleanup();
	}
	explicit Indexed(std::string const &symbol, int i)
	    : atoms_{IndexedAtom{symbol, {i}}}
	{
		cleanup();
	}
	explicit Indexed(std::string const &symbol, int i, int j)
	    : atoms_{IndexedAtom{symbol, {i, j}}}
	{
		cleanup();
	}
	explicit Indexed(std::string const &symbol, int i, int j, int k)
	    : atoms_{IndexedAtom{symbol, {i, j, k}}}
	{
		cleanup();
	}

	/** array-like access to atoms */
	std::vector<IndexedAtom> const &atoms() const { return atoms_; }
	IndexedAtom const &operator[](size_t i) const { return atoms_.at(i); }
	size_t size() const { return atoms_.size(); }
	bool empty() const { return atoms_.empty(); }

	int rank() const { return rank_; }

	/** move out internal storage */
	std::vector<IndexedAtom> release() { return std::move(atoms_); }

	void operator*=(Indexed const &b);
};

/** outer product */
Indexed operator*(Indexed const &a, Indexed const &b);

/** inner product contracting last index of a with first index of b */
Indexed inner_product(Indexed const &a, Indexed const &b);

/** contract indices i and j */
Indexed contract(Indexed const &a, int i, int j);

/** returns list of new term, and number of eta's found */
std::pair<std::vector<Indexed>, int> wick_contract(Indexed const &a,
                                                   std::string const &eta);

/** arbitrary order */
bool operator==(Indexed const &a, Indexed const &b);
bool operator<(Indexed const &a, Indexed const &b);

/** only for checking " == 1" */
bool operator==(Indexed const &a, int b);

} // namespace chalk

template <> struct fmt::formatter<chalk::IndexedAtom>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::IndexedAtom &atom, FormatContext &ctx)
	    -> decltype(ctx.out())
	{
		auto it = ctx.out();
		it = format_to(it, "{}", atom.symbol);
		if (!atom.indices.empty())
			*it++ = '_';
		for (int k : atom.indices)
			it = format_to(it, "{}", k);
		return it;
	}
};

template <> struct fmt::formatter<chalk::Indexed>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Indexed &indexed, FormatContext &ctx)
	    -> decltype(ctx.out())
	{

		// empty product -> '1'
		if (indexed.atoms().empty())
			return format_to(ctx.out(), "1");

		// otherwise list the terms with " + " inbetween
		auto it = ctx.out();
		for (size_t i = 0; i < indexed.atoms().size(); ++i)
		{
			if (i != 0)
				*it++ = ' ';
			it = format_to(it, "{}", indexed.atoms()[i]);
		}
		return it;
	}
};

#endif
