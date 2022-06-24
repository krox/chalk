#include "chalk/numtheory.h"
#include "chalk/parser.h"
#include "util/lexer.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

int main(int argc, char *argv[])
{
	int64_t limit = 1100000000;

	if (argc <= 1)
	{
		fmt::print("usage: sprp <base(s)>\n");
		return -1;
	}

	std::vector<int64_t> bases;
	for (int i = 1; i < argc; ++i)
		bases.push_back(util::parse_int64(argv[i]));

	auto ps = primes(limit);
	size_t k = 0;
	int nFound = 0;
	for (int64_t n = 3; n <= limit; n += 2)
	{

		while (k < ps.size() - 1 && ps[k] < n)
			++k;
		if (ps[k] == n)
			continue;

		bool composite = false;
		for (auto base : bases)
			if (!is_sprp(base, n))
			{
				composite = true;
				break;
			}

		if (!composite)
		{
			fmt::print("{}\n", n);
			nFound += 1;
		}
	}
}
