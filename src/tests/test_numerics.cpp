#include "chalk/floating.h"
#include "chalk/fraction.h"
#include "chalk/integer.h"
#include "chalk/numerics.h"
#include "chalk/polynomial.h"
#include "gtest/gtest.h"
#include <fmt/format.h>
using namespace chalk;

TEST(Numerics, PolyRootFinding)
{
	using Poly = Polynomial<double>;
	auto x = Poly::generator();

	auto poly = (x * x + 1) * (x - 1) * (x + 1);
	auto zeros = roots(poly);
	ASSERT_EQ(zeros.size(), 2);
	EXPECT_DOUBLE_EQ(zeros[0], -1.);
	EXPECT_DOUBLE_EQ(zeros[1], +1.);
}

TEST(Numerics, LegendrePolys)
{
	using Rat = Fraction<Integer>;
	using Poly = Polynomial<Rat>;
	auto x = Poly::generator();

	EXPECT_EQ(fmt::format("{}", jacobiPolynomial(Rat(0), Rat(0), 0, 0, x)),
	          "1");
	EXPECT_EQ(fmt::format("{}", jacobiPolynomial(Rat(0), Rat(0), 1, 0, x)),
	          "x");
	EXPECT_EQ(fmt::format("{}", jacobiPolynomial(Rat(0), Rat(0), 2, 0, x)),
	          "3/2*x^2 - 1/2");
	EXPECT_EQ(fmt::format("{}", jacobiPolynomial(Rat(0), Rat(0), 3, 0, x)),
	          "5/2*x^3 - 3/2*x");
	EXPECT_EQ(fmt::format("{}", jacobiPolynomial(Rat(0), Rat(0), 4, 0, x)),
	          "35/8*x^4 - 15/4*x^2 + 3/8");
	EXPECT_EQ(fmt::format("{}", jacobiPolynomial(Rat(0), Rat(0), 5, 0, x)),
	          "63/8*x^5 - 35/4*x^3 + 15/8*x");
	EXPECT_EQ(fmt::format("{}", jacobiPolynomial(Rat(0), Rat(0), 6, 0, x)),
	          "231/16*x^6 - 315/16*x^4 + 105/16*x^2 - 5/16");
}
