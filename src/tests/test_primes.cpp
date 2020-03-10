#include "chalk/numtheory.h"
#include "fmt/format.h"
#include "util/random.h"
#include <cassert>

#include <random>
using namespace chalk;
using namespace util;

int main()
{
	auto ps = primes(10000);
	assert(ps.size() == 1229);

	size_t k = 0;
	for (int64_t n = 2; n <= ps.back(); ++n)
	{
		if (n == ps[k])
		{
			assert(isPrime(n));
			++k;
		}
		else
			assert(!isPrime(n));
	}

	{
		auto fs = factor(2L * 2 * 2 * 3 * 5 * 13);
		assert((fs == std::vector<int64_t>{2, 2, 2, 3, 5, 13}));
	}

	{
		auto fs = factor(1000000007L * 1000000009L);
		assert(fs[0] == 1000000007L && fs[1] == 1000000009L);
	}

	{
		auto dist = std::uniform_int_distribution<int64_t>(1, INT64_MAX);
		xoshiro256 rng = {};
		for (int i = 0; i < 100; ++i)
		{
			int64_t n = dist(rng);

			auto fs = factor(n);
			fmt::print("{}:", n);
			int64_t m = 1;

			for (auto f : fs)
			{
				assert(isPrime(f));
				m *= f;
				fmt::print(" {}", f);
			}
			assert(m == n);
			fmt::print("\n");
		}
	}

	fmt::print("all done\n");
}
