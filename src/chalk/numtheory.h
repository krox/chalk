#ifndef CHALK_NUMTHEROY_H
#define CHALK_NUMTHEORY_H

/**
 * Fast number theory for 64 bit integers.
 * Also some combinatorics for convenience.
 */

#include "chalk/modular.h"
#include <cstdint>
#include <vector>

namespace chalk {

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
bool isPrime(int64_t n);

/**
 * Find a factor of n using the pollard rho method.
 * Returns either a proper factor of n (which is not necessarily prime itself),
 * or n if none was found. In the latter case try using a different value for c.
 */
int64_t findFactorPollardRho(int64_t n, int64_t c);

/** Find a non-trivial factor of n */
int64_t findFactor(int64_t n);

/** Factor n (>= 1) into primes. Result is in ascending order. */
std::vector<int64_t> factor(int64_t n);

/** Jacobi symbol (a/n) defined for any odd integer n */
int jacobi(int64_t a, int64_t n);

/** factorial(n) = n! = n * (n-1) * ... * 2 * 1 */
int64_t factorial(int64_t n);

/** binomial(n,k) = (n choose k) = n! / (k! * (n-k)!) */
int64_t binomial(int64_t n, int64_t k);

} // namespace chalk

#endif
