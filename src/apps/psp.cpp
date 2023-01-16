#include "CLI/CLI.hpp"
#include "chalk/numtheory.h"
#include "chalk/parser.h"
#include "fmt/format.h"
#include "util/lexer.h"
#include <algorithm>
#include <cassert>
using namespace chalk;

int main(int argc, char *argv[])
{
	int64_t base = 2;
	int64_t limit = 1e9;

	CLI::App app{"generate pseudo primes"};
	app.add_option("--base", base, "base to check");
	app.add_option("--limit", limit, "only check candidates <= limit");

	CLI11_PARSE(app, argc, argv);

	auto table = chalk::prime_table(limit, 2);
	for (int64_t k = 1; k < (int64_t)table[1].size(); ++k)
		if (!table[1][k])
		{
			auto n = 2 * k + 1;
			if (n > limit)
				return 0;
			if (is_prp(base, n))
				fmt::print("{}\n", n);
		}
}
