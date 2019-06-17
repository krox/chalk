#ifndef CHALK_NUMTHEROY_H
#define CHALK_NUMTHEORY_H

#include <cstdint>
#include <vector>

namespace chalk {

/**
 * Basic modular arithmetic.
 * Faster than the naive expressions and safe from overflow.
 * Conditions: 0 <= a,b < m
 */
int64_t addmod(int64_t a, int64_t b, int64_t m); // (a + b) % m
int64_t submod(int64_t a, int64_t b, int64_t m); // (a - b) % m
int64_t mulmod(int64_t a, int64_t b, int64_t m); // (a * b) % m
int64_t powmod(int64_t a, int64_t b, int64_t m); // pow(a, b) % m
int64_t invmod(int64_t a, int64_t m);            // (1 / a) % m

/**
 * Greatest common divisor and least common multiple.
 * Sign of a and b is ignored, result is always >= 0.
 * Convention: gcd(0,x) = gcd(x,0) = abs(x)
 *             lcm(0,x) = lcm(x,0) = 0
 */
int64_t gcd(int64_t a, int64_t b);
int64_t lcm(int64_t a, int64_t b);

/** Compute all primes up to n (inclusive) */
std::vector<int64_t> primes(int64_t n);

/** Tests wether n (>= 3 ) is a strong probable prime to base a (>= 0). */
bool isSPRP(int64_t a, int64_t n);

/** Tests if n is prime number. */
bool isPrime(int64_t);

/**
 * Find a factor of n using the pollard rho method.
 * Returns either a proper factor of n (which is not necessarily prime itself),
 * or n if none was found. In the latter case try using a different value for c.
 */
int64_t findFactor(int64_t n, int64_t c);

/** Factor n (>= 1) into primes. Result is in ascending order. */
std::vector<int64_t> factor(int64_t n);

} // namespace chalk

#endif
