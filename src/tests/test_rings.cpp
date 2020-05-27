#include "chalk/fraction.h"
#include "chalk/free_group.h"
#include "chalk/sparse_polynomial.h"
#include <fmt/format.h>
#include <gtest/gtest.h>
using namespace chalk;

TEST(Polynomial, Univariate)
{
	using R = Fraction<int64_t>;
	auto x = Polynomial<R>::generator();
	assert(fmt::format("{}", (x + 1) * (x - 1)) == "x^2 - 1");
	assert(fmt::format("{}", (x + 1) * (x + 1)) == "x^2 + 2*x + 1");
}

TEST(Polynomial, Multivariate)
{
	auto R = PolynomialRing<Fraction<int64_t>, 2>({"x", "y"});
	auto x = R.generator(0);
	auto y = R.generator(1);
	auto z = (x + y) * (x - y) + R(1);
	assert(fmt::format("{}", z) == "x^2 - y^2 + 1");
	assert(fmt::format("{}", -z) == "-x^2 + y^2 - 1");
	assert(fmt::format("{}", z * 2) == "2*x^2 - 2*y^2 + 2");
	assert(fmt::format("{}", z / 2) == "1/2*x^2 - 1/2*y^2 + 1/2");

	auto R2 = PolynomialRing<Fraction<int64_t>, 4>();
	assert(R2("(x+1)^2") == R2("x^2+x*2+5-4"));
}

TEST(FreeGroup, BasicGroupOperations)
{
	auto a = chalk::FreeProduct::generator(0);
	auto b = chalk::FreeProduct::generator(1);
	auto c = chalk::FreeProduct::generator(2);
	auto d = chalk::FreeProduct::generator(3);
	EXPECT_EQ(fmt::format("{}", a * b * c), "abc");
	EXPECT_EQ(fmt::format("{}", inverse(a * b * c)), "CBA");
	EXPECT_EQ(fmt::format("{}", a * b * c / (d * b * c)), "aD");
}
