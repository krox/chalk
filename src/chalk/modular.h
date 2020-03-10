#ifndef CHALK_MODULAR_H
#define CHALK_MODULAR_H

/**
 * Basic modular arithmetic for 64-bit and 128-bit integers.
 *   - overflow safe
 *   - plenty of asserts
 */

#include <algorithm>
#include <cassert>
#include <cstdint>

namespace chalk {

inline int64_t negmod(int64_t a, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if (a == 0)
		return 0;
	else
		return m - a;
}

inline int64_t invmod(int64_t a, int64_t m)
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

inline int64_t addmod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (b < m - a)
		return a + b;
	else
		return a + b - m;
}

inline int64_t submod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (a >= b)
		return a - b;
	else
		return a - b + m;
}

inline int64_t mulmod(int64_t a, int64_t b, int64_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	__int128_t a_ = a;
	__int128_t b_ = b;
	__int128_t m_ = m;
	return (int64_t)((a_ * b_) % m_);
}

inline int64_t divmod(int64_t a, int64_t b, int64_t m)
{
	b = invmod(b, m);
	return mulmod(a, b, m);
}

inline int64_t powmod(int64_t a, int64_t b, int64_t m)
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

inline __int128_t negmod(__int128_t a, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if (a == 0)
		return 0;
	else
		return m - a;
}

inline __int128_t invmod(__int128_t a, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	__int128_t a0 = m;
	__int128_t a1 = a;
	__int128_t b0 = 0;
	__int128_t b1 = 1;

	while (a1 > 1)
	{
		__int128_t q = a0 / a1;
		__int128_t a2 = a0 - q * a1;
		__int128_t b2 = b0 - q * b1;

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

inline __int128_t addmod(__int128_t a, __int128_t b, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (b < m - a)
		return a + b;
	else
		return a + b - m;
}

inline __int128_t submod(__int128_t a, __int128_t b, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (a >= b)
		return a - b;
	else
		return a - b + m;
}

inline __int128_t mulmod(__int128_t a, __int128_t b, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if (a > b)
		std::swap(a, b);

	__int128_t r = 0;
	for (; a != 0; a >>= 1, b = addmod(b, b, m))
		if (a & 1)
			r = addmod(r, b, m);
	return r;
}

inline __int128_t divmod(__int128_t a, __int128_t b, __int128_t m)
{
	b = invmod(b, m);
	return mulmod(a, b, m);
}

inline __int128_t powmod(__int128_t a, __int128_t b, __int128_t m)
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if (b < 0)
	{
		a = invmod(a, m);
		b = -b;
	}

	__int128_t r = 1;
	for (; b != 0; b >>= 1, a = mulmod(a, a, m))
		if (b & 1)
			r = mulmod(r, a, m);
	return r;
}

} // namespace chalk

#endif
