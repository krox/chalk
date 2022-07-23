#include "CLI/CLI.hpp"
#include "chalk/numtheory.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

int main(int argc, char *argv[])
{
	// numbers of primes below n (for testing, verified by wolfram alpha)
	// 10    ->           4
	// 10^2  ->          25
	// 10^3  ->         168
	// 10^4  ->       1'229
	// 10^5  ->       9'592
	// 10^6  ->      78'498
	// 10^7  ->     664'579
	// 10^8  ->   5'761'455
	// 10^9  ->  50'847'534
	// 10^10 -> 455'052'511

	int64_t max = 1'000'000'000;
	CLI::App app{"list all prime numbers in [2, max]"};
	app.add_option("--max", max);
	CLI11_PARSE(app, argc, argv);
	if (max < 2)
		return 0;

	// NOTE:
	//     * the bit_vectors returned by prime_table take a lot less space than
	//       a std::vector<int64_t> would, so the latter one is never created
	// TODO:
	//     - segmented sieve (and set default max = INT64_MAX)
	//     - '--min' parameter

	auto table = chalk::prime_table(max, 30);
	for (auto p : {2, 3, 5})
		if (p <= max)
			fmt::print("{}\n", p);

	for (int64_t k = 0; k < (int64_t)table[1].size(); ++k)
		for (int64_t x : {1, 7, 11, 13, 17, 19, 23, 29})
			if (table[x][k])
			{
				auto p = k * 30 + x;
				if (p > max)
					return 0;
				fmt::print("{}\n", p);
			}
}
