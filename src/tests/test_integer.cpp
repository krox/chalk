#include "chalk/fraction.h"
#include "chalk/integer.h"
#include "chalk/numtheory.h"
#include "util/random.h"
#include <fmt/format.h>
#include <gtest/gtest.h>
using namespace chalk;
using namespace util;

TEST(Integer, RationalArithmetic)
{
	using R = Fraction<Integer>;
	auto x = R(Integer(3), Integer(5));
	EXPECT_EQ(fmt::format("{}", x), "3/5");
	EXPECT_EQ(fmt::format("{}", x * x), "9/25");
	EXPECT_EQ(fmt::format("{}", x * x * x * x), "81/625");
	EXPECT_EQ(fmt::format("{}", x * x + 1), "34/25");
	EXPECT_EQ(fmt::format("{}", x * x + x), "24/25");
}

TEST(Integer, SmallRationalArithmetic)
{
	using R = Fraction<int64_t>;
	EXPECT_EQ(fmt::format("{}", R(6, 2)), "3");
	EXPECT_EQ(fmt::format("{}", inverse(R(6, 2))), "1/3");
	EXPECT_EQ(fmt::format("{}", R(12, -9)), "-4/3");
	EXPECT_EQ(fmt::format("{}", inverse(R(12, -9))), "-3/4");
	EXPECT_EQ(R(5, 6) + R(1, 6), R(-7, -7));
	EXPECT_EQ(-R(3, 8) - R(1, 16), R(7, -16));
}

TEST(NumberTheory, Primes)
{
	auto ps = primes(10000);
	EXPECT_EQ(ps.size(), 1229);

	size_t k = 0;
	for (int64_t n = 2; n <= ps.back(); ++n)
	{
		if (n == ps[k])
		{
			ASSERT_TRUE(isPrime(n));
			++k;
		}
		else
			ASSERT_FALSE(isPrime(n));
	}
}

TEST(NumberTheory, Factorization)
{
	{
		auto fs = factor(2L * 2 * 2 * 3 * 5 * 13);
		EXPECT_EQ(fs, (std::vector<int64_t>{2, 2, 2, 3, 5, 13}));
	}

	{
		auto fs = factor(1000000007L * 1000000009L);
		EXPECT_EQ(fs, (std::vector<int64_t>{1000000007L, 1000000009L}));
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
				assert(isPrime(f));
				m *= f;
			}
			EXPECT_EQ(m, n);
		}
	}
}

TEST(NumberTheory, ModularArithmetic)
{
	auto a = IntMod(3, 499);
	EXPECT_EQ(a * inverse(a), 1);
	EXPECT_EQ(pow(a, 499 - 1), 1);
	EXPECT_EQ(pow(a, (499 - 1) / 2), -1); // only bc 3 mod 499 is cyclic
}
