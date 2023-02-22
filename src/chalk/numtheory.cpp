#include "chalk/numtheory.h"

#include "util/bit_vector.h"
#include "util/lazy.h"
#include "util/random.h"
#include <algorithm>
#include <cassert>
#include <cmath>
using namespace util;

namespace chalk {

// note on storage requirements:
// a bit_vector indicating primeness of all odd numbers is more efficint than
// storing primes directly:
//     * The ~203M primes below 2^32 take ~775 MiB as a vector<uint32>, but only
//       256 MiB as a bit_vector.
//     * The ~4*10^17 primes below 2^64 would take ~3 ZiB as a vector<uint64>,
//       but "only" 1 ZiB as a bit_vector.
//     * Due to Pi(n)~n/log(n), this ~3:1 ratio holds true at all bit widths.
//     * With a larger wheel size, this advantage only grows.
//       For example with a (2,3,5,7) wheel, it is ~6.5:1

namespace {
auto odd_prime_table_19 = util::synchronized_lazy(
    []() { return odd_prime_table(1LL << 19); }); // 64 KiB
auto odd_prime_table_23 = util::synchronized_lazy(
    []() { return odd_prime_table(1LL << 23); }); // 1 MiB
auto odd_prime_table_27 = util::synchronized_lazy(
    []() { return odd_prime_table(1LL << 27); }); // 16 MiB
auto odd_prime_table_31 = util::synchronized_lazy(
    []() { return odd_prime_table(1LL << 31); }); // 256 MiB

util::bit_vector const &get_odd_prime_table(uint64_t n)
{
	if (n <= 1LL << 19)
		return *odd_prime_table_19;
	if (n <= 1LL << 23)
		return *odd_prime_table_23;
	if (n <= 1LL << 27)
		return *odd_prime_table_27;
	if (n <= 1LL << 31)
		return *odd_prime_table_31;
	assert(false);
}
} // namespace

util::bit_vector odd_prime_table(uint64_t o, uint64_t n)
{
	assert(o <= uint64_t(1) << 62);
	assert(n <= uint64_t(1) << 48);

	if (!o || !n)
		return odd_prime_table(n);

	uint64_t limit = isqrt((o + n) * 2) / 2 + 1;
	auto &base_table = get_odd_prime_table(limit);

	auto table = util::bit_vector(n, true);
	for (uint64_t k = 1; k < limit; ++k)
		if (base_table[k])
		{
			uint64_t p = 2 * k + 1;
			uint64_t start = p * p / 2;
			if (start < o)
				start += (o - start + p - 1) / p * p;
			for (uint64_t s = start; s < o + n; s += p)
				table[s - o] = false;
		}
	return table;
}

util::bit_vector odd_prime_table(uint64_t n)
{
	assert(n <= uint64_t(1) << 48);
	if (!n)
		return {};

	// completely standard odds-only Erathostenes sieve
	auto table = util::bit_vector(n, true);
	table[0] = false; // mark '1' as non-prime
	for (uint64_t k = 1; k * k + k < (n + 1) / 2; ++k)
		if (table[k])
			for (uint64_t s = 2 * k * (k + 1); s < n; s += 2 * k + 1)
				table[s] = false;
	return table;
}

std::vector<util::bit_vector> prime_table(int64_t n, int64_t w)
{
	// fmt::print("init...\n");
	assert(n >= 0);
	assert(w >= 2);                // w=1 case would be logical, but is broken
	assert(n <= int64_t(1) << 56); // just for sanity...

	// round n / n^0.5 / n^0.25 to multiples of w
	n = (n + w) / w * w;
	int64_t n2 = (isqrt(n) + w) / w * w;
	int64_t n4 = (isqrt(n2) + w) / w * w;
	assert(n4 <= n2 && n2 <= n);
	assert(n4 * n4 >= n2 && n2 * n2 >= n);

	auto r = std::vector<util::bit_vector>(w);
	auto inv = std::vector<int64_t>(w);
	for (int64_t i = 0; i < w; ++i)
		if (gcd(i, w) == 1)
		{
			r[i] = bit_vector(n / w, true);
			inv[i] = invmod(i, w);
		}
	r[1][0] = false; // mark '1' as non-prime

	// step 1: generate primes of to ~ n^0.25
	for (int64_t k = 0; k < n4 / w; ++k)
		for (int64_t i = 0; i < w; ++i)
		{
			if (!r[i].size())
				continue;
			if (!r[i][k])
				continue;
			auto p = k * w + i;

			for (int64_t j = 0; j < w; ++j)
			{
				if (!r[j].size())
					continue;
				int64_t f = submod(mulmod(j, inv[i], w), i, w);
				int64_t start = p * (p + f) / w;
				for (int64_t s = start; s < n / w; s += p)
					r[j][s] = false;
			}
		}

	// step 2: generate primes up to ~ n^0.5
	// (same operations, just better loop order)
	for (int64_t i = 0; i < w; ++i)
		for (int64_t j = 0; j < w; ++j)
		{
			if (!r[i].size())
				continue;
			if (!r[j].size())
				continue;
			int64_t f = submod(mulmod(j, inv[i], w), i, w);

			for (int64_t k = n4 / w; k < n2 / w; ++k)
			{
				if (!r[i][k])
					continue;
				auto p = k * w + i;

				int64_t start = p * (p + f) / w;
				for (int64_t s = start; s < n / w; s += p)
					r[j][s] = false;
			}
		}
	return r;
}

std::vector<int64_t> primes(int64_t n)
{
	constexpr int64_t w = 30;
	static constexpr int64_t xs[] = {1, 7, 11, 13, 17, 19, 23, 29};
	auto r = prime_table(n, 30);
	size_t count = 3;
	for (auto x : xs)
		count += r[x].count();

	std::vector<int64_t> primes;
	primes.reserve(count);
	primes.push_back(2);
	primes.push_back(3);
	primes.push_back(5);
	for (int64_t k = 0; k < (int64_t)r[1].size(); ++k)
		for (auto x : xs)
			if (r[x][k])
				primes.push_back(w * k + x);

	// now we might have computed slightly more primes than requested,
	// so we remove them again (simpler than doing it right in the beginning)
	while (!primes.empty() && primes.back() > n)
		primes.pop_back();

	return primes;
}

bool is_prp(int64_t a, int64_t n)
{
	assert(a >= 0);
	assert(n > 1 && n % 2 == 1);
	return powmod(a, n - 1, n) == 1;
}

bool is_sprp(int64_t a, int64_t n)
{
	assert(a >= 0);
	assert(n > 1 && n % 2 == 1);

	// the original definition assumes a < n-1, we use a straight-forward
	// generalization
	a %= n;
	if (a == 0 || a == 1 || a == n - 1)
		return true;

	int s = 0;
	int64_t d = n - 1;
	while (d % 2 == 0)
	{
		d /= 2;
		s += 1;
	}
	// now it is n-1 = d*2^s, t odd

	a = powmod(a, d, n);
	if (a == 1 || a == n - 1)
		return true;

	while (--s)
	{
		a = mulmod(a, a, n);
		if (a == n - 1)
			return true;
	}

	return false;
}

bool is_prime(int64_t n)
{
	// Idea: use SPRP with fixed bases for which the smalles composite pseudo-
	// prime is known. This gives a fast deterministic primality test for
	// small (<= 2^64) numbers. Seee http://priv.ckp.pl/wizykowski/sprp.pdf for
	// details and http://miller-rabin.appspot.com for the actual values used

	// largest product of primes fitting in 63/64 bit
	constexpr int64_t primorial =
	    2L * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47;

	if (n < 2)
		return false;

	if (n < 51)
		return n == 2 || n == 3 || n == 5 || n == 7 || n == 11 || n == 13 ||
		       n == 17 || n == 19 || n == 23 || n == 29 || n == 31 || n == 37 ||
		       n == 41 || n == 43 || n == 47;

	if (gcd(n, primorial) != 1)
		return false;

	if (n < 51 * 51)
		return true;

	if (n < 291831L)
		return is_sprp(126401071349994536L, n);

	if (n < 1050535501L)
		return is_sprp(336781006125L, n) && is_sprp(9639812373923155L, n);

	if (n < 273919523041L)
		return is_sprp(15L, n) && is_sprp(7363882082L, n) &&
		       is_sprp(992620450144556L, n);

	if (n < 47636622961201L)
		return is_sprp(2L, n) && is_sprp(2570940L, n) &&
		       is_sprp(211991001L, n) && is_sprp(3749873356L, n);

	if (n < 3770579582154547L)
		return is_sprp(2L, n) && is_sprp(880937L, n) && is_sprp(2570940L, n) &&
		       is_sprp(610386380L, n) && is_sprp(4130785767L, n);

	if (n < 585226005592931977L)
		return is_sprp(2L, n) && is_sprp(123635709730000L, n) &&
		       is_sprp(9233062284813009L, n) &&
		       is_sprp(43835965440333360L, n) &&
		       is_sprp(761179012939631437L, n) &&
		       is_sprp(1263739024124850375L, n);

	// if (n < UINT64_MAX) // should work for all n < 2^64
	return is_sprp(2L, n) && is_sprp(325L, n) && is_sprp(9375L, n) &&
	       is_sprp(28178L, n) && is_sprp(450775L, n) && is_sprp(9780504L, n) &&
	       is_sprp(1795265022L, n);

	// fall back to randomized test for very large numbers
	/*if (true)
	{
	    util::xoshiro256 rng(0);
	    for (int i = 0; i < 25; ++i)
	        if (!is_sprp(rng(), n))
	            return false;
	    return true;
	}*/
}

// Returns either a proper factor of n (not necessarily prime itself),
// or n if none was found. In the latter case try using a different value for c.
int64_t find_factor_pollard_rho(int64_t n, int64_t c)
{
	assert(n > 0);
	assert(0 < c && c < n);

	int64_t x = 0; // arbitrary start value
	int64_t runLength = 1;

	while (true)
	{
		int64_t y = x;

		for (int64_t i = 0; i < runLength; ++i)
		{
			x = mulmod(x, x, n);
			x = addmod(x, c, n);
			int64_t d = gcd(x - y, n);

			if (d != 1)
				return d;
		}

		runLength *= 2;
	}
}

int64_t find_factor(int64_t n)
{
	assert(n >= 2 && !is_prime(n));

	int64_t d = n;
	for (int64_t c = 1; d == n; ++c)
		d = find_factor_pollard_rho(n, c);
	return d;
}

std::vector<int64_t> factor(int64_t n)
{
	assert(n > 0);
	std::vector<int64_t> f;
	if (n == 1)
		return f;

	f.push_back(n);
	for (size_t i = 0; i < f.size(); ++i)
		while (!is_prime(f[i]))
		{
			int64_t d = find_factor(f[i]);
			f[i] /= d;
			f.push_back(d);
		}

	std::sort(f.begin(), f.end());
	return f;
}

int jacobi(int64_t a, int64_t n)
{
	assert(n > 0 && n % 2 == 1);

	a %= n;
	if (a < 0)
		a += n;
	if (n == 1)
		return 1;

	int r = 1;
	while (true)
	{
		if (a == 0)
			return 0;
		if (a == 1)
			return r;

		if (a % 2 == 0)
		{
			if (n % 8 == 3 || n % 8 == 5)
				r = -r;
			a /= 2;
		}
		else
		{
			if (a % 4 == 3 && n % 4 == 3)
				r = -r;
			std::swap(a, n);
			a %= n;
		}
	}
}

int64_t phi(int64_t n)
{
	assert(n >= 1);
	int64_t last = 0;
	int64_t r = 1;
	for (int64_t p : factor(n))
	{
		r *= p == last ? p : p - 1;
		last = p;
	}
	return r;
}

int64_t ordermod(int64_t a, int64_t m)
{
	assert(0 <= a && a < m);
	if (gcd(a, m) != 1)
		return 0;
	int64_t r = phi(m);
	assert(powmod(a, r, m) == 1);
	auto ps = factor(r);
	for (int64_t p : ps)
		if (powmod(a, r / p, m) == 1)
			r /= p;

	return r;
}

int64_t factorial(int64_t n)
{
	assert(0 <= n && n <= 20); // n=20 is the biggest possible in int64_t
	int64_t r = 1;
	for (int64_t i = 2; i <= n; ++i)
		r *= i;
	return r;
}

int64_t binomial(int64_t n, int64_t k)
{
	assert(n >= 0 && k >= 0);
	if (k > n)
		return 0;
	if (n - k < k)
		k = n - k;
	int64_t r = 1;
	for (int64_t i = 1; i <= k; ++i)
	{
		assert(r <= INT64_MAX / (n + 1 - i));
		r *= (n + 1 - i);
		r /= i;
	}
	return r;
}

} // namespace chalk
