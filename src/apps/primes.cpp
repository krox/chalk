#include "CLI/CLI.hpp"
#include "chalk/numtheory.h"
#include <algorithm>
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

	uint64_t min = 2;
	uint64_t max = uint64_t(1) << 60;
	CLI::App app{"list all prime numbers in [max, max]"};
	app.add_option("--min", min, "lower bound (default=2)");
	app.add_option("--max", max, "upper bound (default=2^60)");
	CLI11_PARSE(app, argc, argv);
	if (min > max)
		return 0;

	if (min <= 2 && 2 <= max)
		fmt::print("2\n");

	for (uint64_t n = min / 2; n != (max + 1) / 2;)
	{
		// chunk_size must be at least ~sqrt(n) in order for the prime sieve to
		// run at full speed
		uint64_t chunk_size = 100 * isqrt(n);
		chunk_size = std::max(chunk_size, uint64_t(1) << 16);
		chunk_size = std::min(chunk_size, (max + 1) / 2 - n);
		auto table = odd_prime_table(n, chunk_size);

		for (uint64_t k = 0; k < chunk_size; ++k)
			if (table[k])
				fmt::print("{}\n", 2 * (n + k) + 1);
		n += chunk_size;
	}
}
