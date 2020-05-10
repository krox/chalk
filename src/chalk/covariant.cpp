#include "chalk/covariant.h"

#include <functional>

namespace chalk {

int count_repeated(CovariantAtom const &a)
{
	int count = 0;
	for (size_t i = 0; i < a.indices.size(); ++i)
		for (size_t j = i + 1; j < a.indices.size(); ++j)
			if (a.indices[i] == a.indices[j])
				++count;
	return count;
}

void move_double_indices(absl::InlinedVector<int, 4> &inds)
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

void sort_atoms(std::vector<CovariantAtom> &atoms)
{
	// index -> first atom id (already locked atoms only)
	absl::flat_hash_map<int, int> index_prio;

	auto comp = [&index_prio](CovariantAtom const &a, CovariantAtom const &b) {
		// ascending name
		if (a.symbol < b.symbol)
			return true;
		if (a.symbol > b.symbol)
			return false;

		// descending number of indices
		if (a.indices.size() != b.indices.size())
			return a.indices.size() > b.indices.size();

		// descending number of repeated indices
		int count_a = count_repeated(a);
		int count_b = count_repeated(b);
		if (count_a != count_b)
			return count_a > count_b;

		// connections to already locked atoms
		int firstA = 999999, firstB = 999999;
		for (int k : a.indices)
			if (index_prio.count(k) && index_prio[k] < firstA)
				firstA = index_prio[k];
		for (int k : b.indices)
			if (index_prio.count(k) && index_prio[k] < firstB)
				firstB = index_prio[k];

		return firstA < firstB;
	};

	std::stable_sort(atoms.begin(), atoms.end(), comp);
	for (int start = 1; start < (int)atoms.size() - 1; ++start)
	{
		for (int k : atoms[start - 1].indices)
			if (index_prio.count(k) == 0)
				index_prio[k] = start;
		std::stable_sort(atoms.begin() + start, atoms.end(), comp);
	}
}

void rename_indices(std::vector<CovariantAtom> &atoms)
{
	// in which atoms (1 or 2) does each index occur?
	absl::flat_hash_map<int, absl::InlinedVector<int, 2>> occs;
	absl::flat_hash_map<int, int> global_pos;
	int pos = 0;
	for (size_t i = 0; i < atoms.size(); ++i)
		for (int k : atoms[i].indices)
		{
			occs[k].push_back(i);
			if (global_pos.count(k) == 0)
				global_pos[k] = pos++;
		}

	// list variables by (first atom, second atom, global pos, old name)
	std::vector<std::tuple<int, int, int, int>> index_list;
	int max_open = 0;
	for (auto &[k, occ] : occs)
	{
		// ignore open indices
		if (occ.size() == 1)
		{
			max_open = std::max(max_open, k);
			continue;
		}
		assert(occ.size() == 2);

		index_list.push_back({occ[0], occ[1], global_pos[k], k});
	}
	std::sort(index_list.begin(), index_list.end());

	// translate all inner indices
	absl::flat_hash_map<int, int> trans;
	for (auto &[first, second, global, old] : index_list)
	{
		(void)first;
		(void)second;
		(void)global;
		trans[old] = ++max_open;
	}
	for (auto &atom : atoms)
		for (auto &k : atom.indices)
			if (trans.count(k))
				k = trans[k];
}

std::vector<int> normalize_indices(std::vector<CovariantAtom> &atoms)
{
	sort_atoms(atoms);
	rename_indices(atoms);

	// count index-occurences to determine which indices are open
	absl::flat_hash_map<int, int> count;
	for (auto &atom : atoms)
		for (int k : atom.indices)
			count[k] += 1;
	std::vector<int> open;
	for (auto &[k, c] : count)
		if (c == 1)
			open.push_back(k);
		else
			assert(c == 2);
	std::sort(open.begin(), open.end());

	return open;
}

bool simplify_delta(std::vector<CovariantAtom> &atoms)
{
	bool change = false;
	for (size_t ai = 0; ai < atoms.size(); ++ai)
	{
		auto &inds = atoms[ai].indices;
		if (atoms[ai].symbol != "delta")
			continue;
		if (inds.size() != 2)
			throw std::runtime_error("delta needs exactly two indices");
		if (inds[0] == inds[1])
			throw std::runtime_error(
			    "'delta_ii' not implemented (unknown dimension)");

		// sort by symmetry
		if (inds[0] > inds[1])
			std::swap(inds[0], inds[1]);

		// try replacing inds[1] by inds[0] or vice versa
		bool found = false;
		for (size_t aj = 0; aj < atoms.size() && !found; ++aj)
			if (aj != ai)
				for (int &k : atoms[aj].indices)
				{
					if (k == inds[0])
					{
						k = inds[1];
						found = true;
						break;
					}
					if (k == inds[1])
					{
						k = inds[0];
						found = true;
						break;
					}
				}

		// if something was replaced, remove the delta atom
		if (found)
		{
			std::swap(atoms[ai], atoms.back());
			atoms.pop_back();
			--ai;
			change = true;
		}
	}
	return change;
}

/** returns true if a<->b can be proved to be symmetric, excluding symbol */
bool is_symmetric(std::vector<CovariantAtom> const &atoms, int a, int b,
                  std::string const &symbol)
{
	std::vector<CovariantAtom> atoms2 = atoms;
	int count_a = 0, count_b = 0; // for sanity check
	for (auto &atom : atoms2)
	{
		if (atom.symbol == symbol)
			continue;
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
	}
	assert(count_a <= 1);
	assert(count_b <= 1);
	if (count_a == 0 || count_b == 0)
		return false;
	std::sort(atoms2.begin(), atoms2.end());
	return atoms2 == atoms;
};

int simplify_antisymmetric(std::vector<CovariantAtom> &atoms,
                           std::string const &symbol)
{
	int sign = 1;
	for (auto &atom : atoms)
	{
		if (atom.symbol != symbol)
			continue;
		auto &inds = atom.indices;

		// simple insertion sort with counting swaps
		for (size_t i = 0; i < inds.size(); ++i)
			for (size_t j = i + 1; j < inds.size(); ++j)
				if (inds[j] < inds[i])
				{
					sign = -sign;
					std::swap(inds[i], inds[j]);
				}

		for (size_t i = 0; i < inds.size(); ++i)
			for (size_t j = i + 1; j < inds.size(); ++j)
				if (is_symmetric(atoms, inds[i], inds[j], symbol))
					return 0;
	}
	return sign;
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

} // namespace chalk
