#include "chalk/indexed.h"

namespace chalk {

bool operator==(IndexedAtom const &a, IndexedAtom const &b)
{
	return a.symbol == b.symbol && a.indices == b.indices;
}

bool operator<(IndexedAtom const &a, IndexedAtom const &b)
{
	if (a.symbol != b.symbol)
		return a.symbol < b.symbol;
	if (a.indices.size() != b.indices.size())
		return a.indices.size() > b.indices.size();
	return a.indices < b.indices;
}

bool operator==(Indexed const &a, Indexed const &b)
{
	return a.atoms() == b.atoms();
}
bool operator<(Indexed const &a, Indexed const &b)
{
	if (a.atoms().size() != b.atoms().size())
		return a.atoms().size() < b.atoms().size();
	return a.atoms() < b.atoms();
}

bool operator==(Indexed const &a, int b)
{
	if (b != 1)
		throw std::runtime_error("invalid use of 'Indexed == int'");
	return a.atoms().empty();
}

namespace {
int count_repeated(IndexedAtom const &a)
{
	int count = 0;
	for (size_t i = 0; i < a.indices.size(); ++i)
		for (size_t j = i + 1; j < a.indices.size(); ++j)
			if (a.indices[i] == a.indices[j])
				++count;
	return count;
}

/** does not require normalized index-names */
bool simplify_deltas(std::vector<IndexedAtom> &atoms)
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

/**
 * This kinda has to solve graph-isomorphism, so its a very crude approximation.
 */
void sort_atoms(std::vector<IndexedAtom> &atoms)
{
	auto comp = [](IndexedAtom const &a, IndexedAtom const &b) {
		// ascending name
		if (a.symbol != b.symbol)
			return a.symbol < b.symbol;

		// descending number of indices
		if (a.indices.size() != b.indices.size())
			return a.indices.size() > b.indices.size();

		// descending number of repeated indices
		int count_a = count_repeated(a);
		int count_b = count_repeated(b);
		if (count_a != count_b)
			return count_a > count_b;

		return a.prio_ > b.prio_;
	};

	assert(atoms.size() <= 30);
	for (auto &atom : atoms)
		atom.prio_ = 0;

	for (size_t start = 0; start < atoms.size() - 1; ++start)
	{
		std::stable_sort(atoms.begin() + start, atoms.end(), comp);
		for (int k : atoms[start].indices)
			for (size_t i = start + 1; i < atoms.size(); ++i)
				for (int l : atoms[i].indices)
					if (k == l)
						atoms[i].prio_ |= 1 << (30 - start);
	}
}

int rename_indices(std::vector<IndexedAtom> &atoms)
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
	absl::InlinedVector<int, 4> open_indices;
	absl::InlinedVector<std::tuple<int, int, int, int>, 4> inner_indices;
	for (auto &[k, occ] : occs)
	{
		if (occ.size() == 1)
			open_indices.push_back(k);
		else if (occ.size() == 2)
			inner_indices.push_back({occ[0], occ[1], global_pos[k], k});
		else
			assert(false);
	}
	std::sort(open_indices.begin(), open_indices.end());
	std::sort(inner_indices.begin(), inner_indices.end());

	// translate all indices
	int counter = 0;
	absl::flat_hash_map<int, int> trans;
	for (int k : open_indices)
		trans[k] = counter++;
	for (auto &[first, second, global, old] : inner_indices)
	{
		(void)first;
		(void)second;
		(void)global;
		trans[old] = counter++;
	}
	for (auto &atom : atoms)
		for (auto &k : atom.indices)
			k = trans.at(k);
	return (int)open_indices.size();
}

} // namespace

void Indexed::cleanup()
{
	simplify_deltas(atoms_);
	sort_atoms(atoms_);
	rank_ = rename_indices(atoms_);
}

Indexed operator*(Indexed const &a, Indexed const &b)
{
	std::vector<IndexedAtom> atoms;
	atoms.reserve(a.atoms().size() + b.atoms().size());
	atoms = a.atoms();

	// append the two atom lists while renaming the indices of b
	for (auto &atom : b.atoms())
	{
		atoms.push_back(atom);
		for (int &k : atoms.back().indices)
			k += 1000000;
	}

	return Indexed(std::move(atoms));
}

void Indexed::operator*=(Indexed const &b)
{
	atoms_.reserve(atoms_.size() + b.atoms().size());

	for (auto &atom : b.atoms())
	{
		atoms_.push_back(atom);
		for (int &k : atoms_.back().indices)
			k += 1000000;
	}
	cleanup();
}

Indexed inner_product(Indexed const &a, Indexed const &b)
{
	assert(a.rank() >= 1 && b.rank() >= 1);
	std::vector<IndexedAtom> atoms;
	atoms.reserve(a.atoms().size() + b.atoms().size());
	atoms = a.atoms();

	for (auto &atom : b.atoms())
	{
		atoms.push_back(atom);
		for (int &k : atoms.back().indices)
		{
			if (k == 0)
				k = a.rank() - 1;
			else
				k += 1000000;
		}
	}

	return Indexed(std::move(atoms));
}

Indexed contract(Indexed const &a, int i, int j)
{
	assert(0 <= i && i < a.rank());
	assert(0 <= j && j < a.rank());
	assert(i != j);
	auto atoms = a.atoms();
	int found = 0;
	for (auto &atom : atoms)
		for (int &k : atom.indices)
			if (k == i)
			{
				found += 1;
				k = j;
			}
	assert(found == 1);
	return Indexed(std::move(atoms));
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

std::pair<std::vector<Indexed>, int> wick_contract(Indexed const &a,
                                                   std::string const &eta)
{
	// split term into eta / non-eta
	std::vector<int> indices;
	std::vector<IndexedAtom> other;
	for (auto const &atom : a.atoms())
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
		return {{}, (int)indices.size()};

	// otherwise, remainder * wick-contraction terms
	std::vector<Indexed> result;
	for (auto &l : wick_list((int)indices.size()))
	{
		assert(l.size() == indices.size() / 2);
		std::vector<IndexedAtom> tmp;
		tmp.reserve(other.size() + l.size());
		tmp = other;
		for (auto &p : l)
			tmp.push_back({"delta", {indices[p.first], indices[p.second]}});
		result.push_back(Indexed(std::move(tmp)));
	}

	return {std::move(result), (int)indices.size()};
}

} // namespace chalk
