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
	uint64_t limit = INT64_MAX;
	bool strong = true;

	CLI::App app{"check probable primes"};
	app.add_option("--base", base, "base to check");
	app.add_option("--limit", limit, "only check candidates <= limit");
	app.add_flag("--strong", strong,
	             "only report strong pseudo primes (default = true)");

	CLI11_PARSE(app, argc, argv);

	for (uint64_t n = 3; n <= limit; n += 2)
	{

		if (!strong && !is_prp(base, n))
			continue;
		if (strong && !is_sprp(base, n))
			continue;
		if (is_prime(n))
			continue;

		fmt::print("{}\n", n);
	}
}
