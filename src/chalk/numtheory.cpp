#include "chalk/numtheory.h"

#include "util/bit_vector.h"
#include "util/random.h"
#include <algorithm>
#include <cassert>
#include <cmath>
using namespace util;

namespace chalk {

int64_t gcd(int64_t a, int64_t b)
{
	while (true)
	{
		if (a == 0)
			return b >= 0 ? b : -b;
		b %= a;

		if (b == 0)
			return a >= 0 ? a : -a;
		a %= b;
	}
}

int64_t lcm(int64_t a, int64_t b)
{
	if (a == 0 || b == 0)
		return 0;
	else
		return a / gcd(a, b) * b;
}

int128_t gcd(int128_t a, int128_t b)
{
	while (true)
	{
		if (a == 0)
			return b >= 0 ? b : -b;
		b %= a;

		if (b == 0)
			return a >= 0 ? a : -a;
		a %= b;
	}
}

int128_t lcm(int128_t a, int128_t b)
{
	if (a == 0 || b == 0)
		return 0;
	else
		return a / gcd(a, b) * b;
}

std::vector<int64_t> primes(int64_t n)
{
	// Idea: Sieve of Eratosthenes with multiples of 2 and 3 excluded

	if (n < 0)
		n = 0;

	// (exclusive) limit for relevant prime divisors
	int64_t limit = (int64_t)sqrt((double)n) + 1;
	assert(limit * limit > n);

	// excluding 2 and 3 as special cases, all primes have the form 6*k +- 1
	auto b5 = bit_vector(n / 6); // b5[k] represents 6*k+5
	auto b7 = bit_vector(n / 6); // b7[k] represents 6*k+7

	// mark all non-primes
	for (int64_t k = 0; k < n / 6; ++k)
	{
		if (!b5[k])
		{
			int64_t p = 6 * k + 5;

			if (p >= limit)
				break;

			for (long s = (p * (p + 2) - 5) / 6; s < (int64_t)b5.size(); s += p)
				b5[s] = true;

			for (long s = (p * p - 7) / 6; s < (int64_t)b7.size(); s += p)
				b7[s] = true;
		}

		if (!b7[k])
		{
			int64_t p = 6 * k + 7;

			if (p >= limit)
				break;

			for (long s = (p * (p + 4) - 5) / 6; s < (int64_t)b5.size(); s += p)
				b5[s] = true;

			for (long s = (p * p - 7) / 6; s < (int64_t)b7.size(); s += p)
				b7[s] = true;
		}
	}

	// collect primes into array
	std::vector<int64_t> primes;
	primes.reserve(b5.count(false) + b7.count(false) + 2);
	primes.push_back(2);
	primes.push_back(3);
	for (int64_t k = 0; k < n / 6; ++k)
	{
		if (!b5[k])
			primes.push_back(6 * k + 5);
		if (!b7[k])
			primes.push_back(6 * k + 7);
	}

	// now we might have computed slightly more primes than requested,
	// so we remove them again (simpler than doing it right in the beginning)
	while (!primes.empty() && primes.back() > n)
		primes.pop_back();

	return primes;
}

bool isSPRP(int64_t a, int64_t n)
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

bool isPrime(int64_t n)
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
		return isSPRP(126401071349994536L, n);

	if (n < 1050535501L)
		return isSPRP(336781006125L, n) && isSPRP(9639812373923155L, n);

	if (n < 273919523041L)
		return isSPRP(15L, n) && isSPRP(7363882082L, n) &&
		       isSPRP(992620450144556L, n);

	if (n < 47636622961201L)
		return isSPRP(2L, n) && isSPRP(2570940L, n) && isSPRP(211991001L, n) &&
		       isSPRP(3749873356L, n);

	if (n < 3770579582154547L)
		return isSPRP(2L, n) && isSPRP(880937L, n) && isSPRP(2570940L, n) &&
		       isSPRP(610386380L, n) && isSPRP(4130785767L, n);

	if (n < 585226005592931977L)
		return isSPRP(2L, n) && isSPRP(123635709730000L, n) &&
		       isSPRP(9233062284813009L, n) && isSPRP(43835965440333360L, n) &&
		       isSPRP(761179012939631437L, n) &&
		       isSPRP(1263739024124850375L, n);

	// if (n < UINT64_MAX) // should work for all n < 2^64
	return isSPRP(2L, n) && isSPRP(325L, n) && isSPRP(9375L, n) &&
	       isSPRP(28178L, n) && isSPRP(450775L, n) && isSPRP(9780504L, n) &&
	       isSPRP(1795265022L, n);

	// fall back to randomized test for very large numbers
	/*if (true)
	{
	    util::xoshiro256 rng(0);
	    for (int i = 0; i < 25; ++i)
	        if (!isSPRP(rng(), n))
	            return false;
	    return true;
	}*/
}

int64_t findFactorPollardRho(int64_t n, int64_t c)
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

int64_t findFactor(int64_t n)
{
	assert(n >= 2 && !isPrime(n));

	int64_t d = n;
	for (int64_t c = 1; d == n; ++c)
		d = findFactorPollardRho(n, c);
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
		while (!isPrime(f[i]))
		{
			int64_t d = findFactor(f[i]);
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
