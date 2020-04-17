#include "chalk/numerics.h"
#include "chalk/polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

int main()
{
	using Poly = Polynomial<double>;
	auto x = Poly::generator();

	// auto poly = (x - 5) * (x + 7);
	auto poly = (x * x + 1) * (x - 1) * (x + 1);
	fmt::print("{}\n", poly);
	for (auto x : roots(poly))
		fmt::print("{}\n", x);
	roots(poly);
}
