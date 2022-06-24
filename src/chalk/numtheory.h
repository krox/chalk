#pragma once

// Fast number theory for 64 bit integers.
// Also some combinatorics for convenience.

#include "chalk/modular.h"
#include <cstdint>
#include <vector>

namespace chalk {

// cancel common factors and make a non-negative
void removeCommonFactor(int64_t &a, int64_t &b);
void removeCommonFactor(int128_t &a, int128_t &b);

// compute all primes up to n (inclusive)
std::vector<int64_t> primes(int64_t n);

// tests if n (>= 3 ) is a strong probable prime to base a (>= 0)
bool is_sprp(int64_t a, int64_t n);

// tests if n is prime number
bool is_prime(int64_t n);

// find a non-trivial factor of n
int64_t find_factor(int64_t n);

// Factor n (>= 1) into primes. Result is in ascending order.
std::vector<int64_t> factor(int64_t n);

// Jacobi symbol (a/n) defined for any odd integer n
int jacobi(int64_t a, int64_t n);

// Euler's totient function
int64_t phi(int64_t a);

// order of a mod m. returns 0 if gcd(a,b) != 1
int64_t ordermod(int64_t a, int64_t m);

// factorial(n) = n! = n * (n-1) * ... * 2 * 1
int64_t factorial(int64_t n);

// binomial(n,k) = (n choose k) = n! / (k! * (n-k)!)
int64_t binomial(int64_t n, int64_t k);

} // namespace chalk
