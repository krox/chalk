#include "chalk/fraction.h"
#include "chalk/integer.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

int main()
{

	using R = Fraction<Integer>;
	auto x = R(Integer(3), Integer(5));
	assert(fmt::format("{}", x) == "3/5");
	assert(fmt::format("{}", x * x) == "9/25");
	assert(fmt::format("{}", x * x * x * x) == "81/625");
	assert(fmt::format("{}", x * x + 1) == "34/25");
	assert(fmt::format("{}", x * x + x) == "24/25");

	fmt::print("done\n");
}
