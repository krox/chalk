#include "chalk/numtheory.h"

#include "util/bitset.h"
#include <algorithm>
#include <cassert>
#include <cmath>
using namespace util;

namespace chalk {

int64_t addmod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (b < m - a)
		return a + b;
	else
		return a + b - m;
}

int64_t submod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (a >= b)
		return a - b;
	else
		return a - b + m;
}

int64_t mulmod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (a > b)
		std::swap(a, b);

	int64_t r = 0;
	for (; a != 0; a >>= 1, b = addmod(b, b, m))
		if (a & 1)
			r = addmod(r, b, m);
	return r;
}

int64_t powmod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if (b < 0)
	{
		a = invmod(a, m);
		b = -b;
	}

	int64_t r = 1;
	for (; b != 0; b >>= 1, a = mulmod(a, a, m))
		if (b & 1)
			r = mulmod(r, a, m);
	return r;
}

int64_t invmod(int64_t a, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	int64_t a0 = m;
	int64_t a1 = a;
	int64_t b0 = 0;
	int64_t b1 = 1;

	while (a1 > 1)
	{
		int64_t q = a0 / a1;
		int64_t a2 = a0 - q * a1;
		int64_t b2 = b0 - q * b1;

		a0 = a1;
		a1 = a2;
		b0 = b1;
		b1 = b2;
	}

	if (b1 < 0)
		b1 += m;
	assert(0 <= b1 && b1 < m);
	return b1;
}

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

std::vector<int64_t> primes(int64_t n)
{
	// Idea: Sieve of Eratosthenes with multiples of 2 and 3 excluded

	if (n < 0)
		n = 0;

	// (exclusive) limit for relevant prime divisors
	int64_t limit = (int64_t)sqrt((double)n) + 1;
	assert(limit * limit > n);

	// excluding 2 and 3 as special cases, all primes have the form 6*k +- 1
	auto b5 = bitset(n / 6); // b5[k] represents 6*k+5
	auto b7 = bitset(n / 6); // b7[k] represents 6*k+7

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

	int64_t s = __builtin_ctz(n - 1);
	int64_t d = (n - 1) >> s;
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

	if (n < 2)
		return false;

	if (n % 2 == 0)
		return n == 2;

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

	if (true) // should work for all n < 2^64
		return isSPRP(2L, n) && isSPRP(325L, n) && isSPRP(9375L, n) &&
		       isSPRP(28178L, n) && isSPRP(450775L, n) && isSPRP(9780504L, n) &&
		       isSPRP(1795265022L, n);
}

int64_t findFactor(int64_t n, int64_t c)
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
			long d = gcd(x - y, n);

			if (d != 1)
				return d;
		}

		runLength *= 2;
	}
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
			int64_t d = f[i];
			for (int64_t c = 1; d == f[i]; ++c)
				d = findFactor(f[i], c);

			f[i] /= d;
			f.push_back(d);
		}

	std::sort(f.begin(), f.end());
	return f;
}

} // namespace chalk
