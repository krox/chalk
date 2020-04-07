#include "chalk/covariant.h"
#include "chalk/fraction.h"
#include "chalk/polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

int main()
{

	using R = Fraction<int64_t>;
	using Cov = Covariant<R>;
	auto f = Cov("eta", 1) / 2 + Cov("S", 1, 2, 2);
	fmt::print("{}\n", f);
	fmt::print("{}\n", f * f);
	fmt::print("{}\n", f * (Cov("eta", 1) / 2 - Cov("S", 1, 2, 2)) + Cov(7));
	fmt::print("done\n");
}
