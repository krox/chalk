#include "catch2/catch_test_macros.hpp"

#include "chalk/fraction.h"
#include "chalk/integer.h"
#include "chalk/numtheory.h"
#include "util/random.h"
#include <fmt/format.h>
using namespace chalk;
using namespace util;

TEST_CASE("rational arithmetic", "[integer]")
{
	using R = Fraction<Integer>;
	auto x = R(Integer(3), Integer(5));
	CHECK(fmt::format("{}", x) == "3/5");
	CHECK(fmt::format("{}", x * x) == "9/25");
	CHECK(fmt::format("{}", x * x * x * x) == "81/625");
	CHECK(fmt::format("{}", x * x + 1) == "34/25");
	CHECK(fmt::format("{}", x * x + x) == "24/25");
}

TEST_CASE("small rational arithmetic", "[integer]")
{
	using R = Fraction<int64_t>;
	CHECK(fmt::format("{}", R(6, 2)) == "3");
	CHECK(fmt::format("{}", inverse(R(6, 2))) == "1/3");
	CHECK(fmt::format("{}", R(12, -9)) == "-4/3");
	CHECK(fmt::format("{}", inverse(R(12, -9))) == "-3/4");
	CHECK(R(5, 6) + R(1, 6) == R(-7, -7));
	CHECK(-R(3, 8) - R(1, 16) == R(7, -16));
}

TEST_CASE("primes", "[numtheory][integer]")
{
	auto ps = primes(10000);
	REQUIRE(ps.size() == 1229);

	size_t k = 0;
	for (int64_t n = 2; n <= ps.back(); ++n)
	{
		if (n == ps[k])
		{
			CHECK(is_prime(n));
			++k;
		}
		else
			CHECK(!is_prime(n));
	}
}

TEST_CASE("factorization", "[numtheory][integer]")
{
	{
		auto fs = factor(2L * 2 * 2 * 3 * 5 * 13);
		CHECK(fs == std::vector<int64_t>{2, 2, 2, 3, 5, 13});
	}

	{
		auto fs = factor(1000000007L * 1000000009L);
		CHECK(fs == std::vector<int64_t>{1000000007L, 1000000009L});
	}

	{
		auto dist = std::uniform_int_distribution<int64_t>(1, INT64_MAX);
		xoshiro256 rng = {};
		for (int i = 0; i < 100; ++i)
		{
			int64_t n = dist(rng);

			auto fs = factor(n);
			int64_t m = 1;

			for (auto f : fs)
			{
				assert(is_prime(f));
				m *= f;
			}
			CHECK(m == n);
		}
	}
}

TEST_CASE("square numbers", "[numtheory][integer]")
{
	auto check = [](int64_t n) {
		auto s = isqrt(n);
		CHECK(s * s <= n);
		CHECK(s * s > n - 1 - 2 * s);
	};

	CHECK(!is_square(-1));
	CHECK(is_square(0));
	CHECK(is_square(1));
	CHECK(!is_square(2));
	CHECK(!is_square(3));
	CHECK(is_square(4));
	CHECK(!is_square(5));
	for (int i = 0; i < 1000; ++i)
	{
		check(i);
		check(INT64_MAX - i);
		check(i * i);
		if (i)
			check(i * i - 1);
		check((3037000499 - i) * (3037000499 - i));
		check((3037000499 - i) * (3037000499 - i) - 1);
	}

	xoshiro256 rng;
	for (int i = 0; i < 1000; ++i)
	{
		check(rng.uniform<int64_t>(0, INT64_MAX));
		auto x = rng.uniform<int64_t>(0, 3037000499);
		check(x * x);
		check(x * x - 1);
	}
}

TEST_CASE("ipow, iroot, is_power", "[numtheory][integer][power]")
{
	CHECK(ipow(0, 1) == 0);
	CHECK(ipow(1, 0) == 1);
	CHECK(ipow(2, 0) == 1);
	CHECK(ipow(3, 5) == 3 * 3 * 3 * 3 * 3);
	CHECK(is_power(0));
	CHECK(is_power(1));
	CHECK(!is_power(2));
	CHECK(!is_power(3));
	CHECK(is_power(4));
	CHECK(!is_power(5));
	CHECK(!is_power(6));
	CHECK(!is_power(7));
	CHECK(is_power(8));

	for (int64_t base = 2; base < 100; ++base)
		for (int64_t x = base * base, e = 2; x <= INT64_MAX / base;
		     x *= base, e += 1)
		{
			CHECK(ipow(base, e) == x);
			CHECK(iroot(x, e) == base);
			CHECK(iroot(x - 1, e) == base - 1);
			CHECK(is_power(x));
		}
}

TEST_CASE("modular arithmetic", "[numtheory][integer]")
{
	auto a = IntMod(3, 499);
	CHECK(a * inverse(a) == 1);
	CHECK(pow(a, 499 - 1) == 1);
	CHECK(pow(a, (499 - 1) / 2) == -1); // only bc 3 mod 499 is cyclic
}

TEST_CASE("misc numtheory", "[numtheory][integer]")
{
	CHECK(phi(2 * 2 * 3 * 5 * 19) == 2 * 2 * 4 * 18);

	auto orderNaive = [](int64_t a, int64_t m) -> int64_t {
		assert(0 <= a && a < m);
		int64_t b = a;
		int64_t r = 1;
		while (b != 1)
		{
			b = mulmod(b, a, m);
			++r;
			if (r > m)
				return 0;
		}
		return r;
	};

	for (int64_t m = 2; m < 300; ++m)
		for (int64_t a = 0; a < m; ++a)
			CHECK(ordermod(a, m) == orderNaive(a, m));
}
