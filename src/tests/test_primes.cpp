#include "chalk/numtheory.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

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

	fmt::print("all done\n");
}
