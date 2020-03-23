#include "chalk/fraction.h"
#include "chalk/polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

int main()
{
	{
		using R = Fraction<int64_t>;
		assert(fmt::format("{}", R(6, 2)) == "3");
		assert(fmt::format("{}", inverse(R(6, 2))) == "1/3");
		assert(fmt::format("{}", R(12, -9)) == "-4/3");
		assert(fmt::format("{}", inverse(R(12, -9))) == "-3/4");
		assert(R(5, 6) + R(1, 6) == R(-7, -7));
		assert(-R(3, 8) - R(1, 16) == R(7, -16));
	}

	{
		auto R = PolynomialRing<Fraction<int64_t>, 2>({"x", "y"});
		auto x = R.generator(0);
		auto y = R.generator(1);
		auto z = (x + y) * (x - y) + R(1);
		assert(fmt::format("{}", z) == "x^2 - y^2 + 1");
		assert(fmt::format("{}", -z) == "-x^2 + y^2 - 1");
		assert(fmt::format("{}", z * 2) == "2*x^2 - 2*y^2 + 2");
		assert(fmt::format("{}", z / 2) == "1/2*x^2 - 1/2*y^2 + 1/2");
	}

	fmt::print("done\n");
}
