#include "CLI/CLI.hpp"
#include "chalk/numtheory.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

bool is_wieferich(int64_t p, int64_t a)
{
	a %= p * p;
	return powmod(a, p - 1, p * p) == 1;
}

int main(int argc, char *argv[])
{
	// complete list for p <= 10^7 (for testing purposes):
	// a=2: 1093, 3511 (most studied case. no other examples up to 2^64
	//                  according to literature)
	// a=3: 11, 1006003
	// a=4: 1093, 3511
	// a=5: 2, 20771, 40487
	// a=6: 66161, 534851, 3152573
	// a=7: 5, 491531

	int64_t max = 10'000'000;
	int64_t base = 2;
	CLI::App app{"compute Wieferich primes"};
	app.add_option("--base", base);
	app.add_option("--max", max);
	CLI11_PARSE(app, argc, argv);
	if (max < 2)
		return 0;

	auto table = chalk::prime_table(max, 30);
	for (auto p : {2, 3, 5})
		if (p <= max && is_wieferich(p, base))
			fmt::print("{}\n", p);

	for (int64_t k = 0; k < (int64_t)table[1].size(); ++k)
		for (int64_t x : {1, 7, 11, 13, 17, 19, 23, 29})
			if (table[x][k])
			{
				auto p = k * 30 + x;
				if (p > max)
					return 0;
				if (is_wieferich(p, base))
					fmt::print("{}\n", p);
			}
}
