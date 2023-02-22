#include "CLI/CLI.hpp"
#include "chalk/numtheory.h"
#include "fmt/format.h"
#include "util/threadpool.h"
#include <algorithm>
#include <cassert>
#include <deque>
using namespace chalk;

// md5sum up to 10M = 16fbfb07d90651519225ed7fd84c6713
//             100M = 631b5eb6b2385d6f540d5cfa8e0789a4
//               1B = 0f278d34fab1724313ef3d291fd2c6da

std::vector<uint64_t> find_psps(int64_t base, uint64_t offset, uint64_t count)
{
	std::vector<uint64_t> r;
	auto table = odd_prime_table(offset, count);
	for (uint64_t k = 0; k < count; ++k)
		if (!table[k] && is_prp(base, 2 * (k + offset) + 1))
			r.push_back(2 * (k + offset) + 1);
	return r;
}

int main(int argc, char *argv[])
{
	uint64_t base = 2;
	uint64_t limit = 1LL << 62; // practically infinite
	int nthreads = 1;

	CLI::App app{"generate pseudo primes"};
	app.add_option("--base", base, "base to check");
	app.add_option("--limit", limit, "only check candidates <= limit");
	app.add_option("--threads,-j", nthreads, "number of worker threads to use");
	CLI11_PARSE(app, argc, argv);

	util::ThreadPool threadpool(nthreads);

	std::deque<std::future<std::vector<uint64_t>>> q;

	uint64_t m = 1; // (exclusive) limit up to which we queued jobs already
	while (true)
	{
		while (m != (limit + 1) / 2 &&
		       q.size() < (size_t)threadpool.num_threads() * 4)
		{
			// chunk-size should be at least sqrt(current limit), otherwise the
			// the prime-sieve is not running at peak performance. Scaling it
			// dynamically seems nice to make it responsive in the beginning and
			// also low-overhead for long runs.
			uint64_t chunk_size = 100 * isqrt(m);
			chunk_size = std::max(chunk_size, uint64_t(1) << 16);
			chunk_size = std::min(chunk_size, (limit + 1) / 2 - m);
			q.push_back(threadpool.async(&find_psps, base, m, chunk_size));
			m += chunk_size;
		}
		if (q.empty())
			break;
		assert(q.front().valid());
		auto list = q.front().get();
		q.pop_front();
		for (uint64_t n : list)
			fmt::print("{}\n", n);
	}
}
