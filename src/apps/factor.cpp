#include "CLI/CLI.hpp"
#include "chalk/integer.h"
#include "chalk/numtheory.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

bool verbose = false;
std::vector<int> small_primes = {2,  3,  5,  7,  11, 13, 17, 19, 23,
                                 29, 31, 37, 41, 43, 47, 53, 59, 61,
                                 67, 71, 73, 79, 83, 89, 97};

Integer findFactorPollardRho(Integer n, int64_t c)
{
	// assert(0 < c && c < n);

	if (verbose)
		fmt::print("  starting pollard-rho with n = {}, c = {}\n", n, c);

	Integer x = Integer(0); // arbitrary start value
	int64_t runLength = 1;

	while (true)
	{
		Integer y = x;
		if (verbose)
			fmt::print("    trying run-length = {}\n", runLength);

		for (int64_t i = 0; i < runLength; ++i)
		{
			x = (x * x + c) % n;
			Integer d = gcd(x - y, n);

			if (!(d == 1))
			{
				if (verbose)
					fmt::print("    found factor d = {}\n", d);
				return d;
			}
		}

		runLength *= 2;
	}
}

Integer findFactor(Integer n)
{
	// assert(!(n < 2) && !is_prime(n));

	Integer d = n;
	for (int64_t c = 1; d == n; ++c)
		d = findFactorPollardRho(n, c);
	return d;
}

std::vector<Integer> factor(Integer n)
{
	if (verbose)
		fmt::print("starting factorization of n = {}\n", n);

	// assert(n > 0);
	std::vector<Integer> f;

	if (n.fits_int())
	{
		if (verbose)
			fmt::print("n is small, using non-GMP algorithm.\n");
		for (int64_t p : factor(n.to_int()))
			f.push_back(Integer(p));
		return f;
	}

	for (int p : small_primes)
		while (n % p == 0)
		{
			if (verbose)
				fmt::print("  small factor p = {}\n", p);
			f.emplace_back(p);
			n /= p;
		}

	if (n == 1)
		return f;

	f.push_back(n);
	for (size_t i = f.size() - 1; i < f.size(); ++i)
		while (!is_prime(f[i]))
		{
			Integer d = findFactor(f[i]);
			f[i] /= d;
			f.push_back(d);
		}

	std::sort(f.begin(), f.end());
	return f;
}

int main(int argc, char *argv[])
{
	// testing number: 2044000015922760011571288401882419756046486692
	// ( three 12-digit prime factors and some small ones )

	std::vector<Integer> numbers;

	CLI::App app{
	    "Factorizes positive integer(s) into primes. If none are given on "
	    "the command line, numbers are read from stdin. Usage and "
	    "output syntax are the same as 'factor' from GNU coreutils."};
	app.add_option("number", numbers, "number(s) to factorize")
	    ->check(CLI::PositiveNumber);
	app.add_flag("--verbose", verbose,
	             "print some debug infos on the used algorithm");
	CLI11_PARSE(app, argc, argv);

	auto process = [](Integer const &n) {
		fmt::print("{}:", n);
		for (auto &a : factor(n))
			fmt::print(" {}", a);
		fmt::print("\n");
	};

	if (numbers.empty())
	{
		std::string s;
		while (std::cin >> s)
		{
			auto n = Integer(s);
			if (!(n > 0))
			{
				fmt::print("ERROR: invalid number '{}'\n", s);
				return -1;
			}
			process(n);
		}
	}
	else
	{
		for (auto &n : numbers)
			process(n);
	}
}
