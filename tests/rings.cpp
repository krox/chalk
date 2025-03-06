#include "catch2/catch_test_macros.hpp"

#include "chalk/fraction.h"
#include "chalk/free_group.h"
#include "chalk/sparse_polynomial.h"
#include <fmt/format.h>
using namespace chalk;

TEST_CASE("univariate polynomials", "[polynomial]")
{
	using R = Fraction<int64_t>;
	auto x = Polynomial<R>::generator();
	auto y = Polynomial<Polynomial<R, 'y'>>(Polynomial<R, 'y'>::generator());
	CHECK(fmt::format("{}", (x + 1) * (x - 1)) == "x^2 - 1");
	CHECK(fmt::format("{}", (y + 1) * (y + 1)) == "y^2 + 2*y + 1");
	CHECK(fmt::format("{}", (y + x) * (x - y)) == "x^2 - y^2");
}

TEST_CASE("multivariate polys", "[polynomial]")
{
	auto R = PolynomialRing<Fraction<int64_t>, 2>({"x", "y"});
	auto x = R.generator(0);
	auto y = R.generator(1);
	auto z = (x + y) * (x - y) + R(1);
	CHECK(fmt::format("{}", z) == "x^2 - y^2 + 1");
	CHECK(fmt::format("{}", -z) == "-x^2 + y^2 - 1");
	CHECK(fmt::format("{}", z * 2) == "2*x^2 - 2*y^2 + 2");
	CHECK(fmt::format("{}", z / 2) == "1/2*x^2 - 1/2*y^2 + 1/2");

	auto R2 = PolynomialRing<Fraction<int64_t>, 4>();
	CHECK(R2("(x+1)^2") == R2("x^2+x*2+5-4"));
}

TEST_CASE("free group", "")
{
	auto a = chalk::FreeProduct::generator("a");
	auto b = chalk::FreeProduct::generator("b");
	auto c = chalk::FreeProduct::generator("c");
	auto d = chalk::FreeProduct::generator("d");
	CHECK(fmt::format("{}", a * b * c) == "a * b * c");
	CHECK(fmt::format("{}", inverse(a * b * c)) == "c^-1 * b^-1 * a^-1");
	CHECK(fmt::format("{}", a * b * c / (d * b * c)) == "a * d^-1");
}

TEST_CASE("free abelian group", "")
{
	auto a = chalk::FreeAbelianProduct::generator("a");
	auto b = chalk::FreeAbelianProduct::generator("b");

	CHECK(fmt::format("{}", a * b * a) == "a^2 * b");
}
