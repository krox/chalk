#include "chalk/fraction.h"
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

	fmt::print("done\n");
}
