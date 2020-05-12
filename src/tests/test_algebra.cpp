#include "chalk/bch.h"
#include "chalk/fraction.h"
#include "chalk/free_group.h"
#include "chalk/group_ring.h"
#include "chalk/polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

using Rat = Fraction<int64_t>;
using Algebra = FreeAlgebra<Rat>;
using R = Polynomial<Algebra>;

R exp(R const &a)
{
	assert(a.max_order() < 100);
	auto inc = R(1);
	R r = 0;
	for (int i = 1; !(inc == 0); ++i)
	{
		r += inc;
		inc *= a;
		inc /= i;
	}
	return r;
}

R comm(R const &a, R const &b) { return a * b - b * a; }

int main()
{
	// basic tests of free group
	{
		auto a = FreeProduct::generator(0);
		auto b = FreeProduct::generator(1);
		auto c = FreeProduct::generator(2);
		auto d = FreeProduct::generator(3);
		assert(fmt::format("{}", a * b * c) == "abc");
		assert(fmt::format("{}", inverse(a * b * c)) == "CBA");
		assert(fmt::format("{}", a * b * c / (d * b * c)) == "aD");
	}

	// test BCH formula by explicitly expanding everything
	for (int order = 1; order <= 8; ++order)
	{
		auto a = truncate(R::generator() * Algebra::generator(0), order);
		auto b = truncate(R::generator() * Algebra::generator(1), order);
		auto error = exp(a) * exp(b) - exp(bch(a, b, order, comm));
		assert(fmt::format("{}", error) ==
		       fmt::format("0 + O(x^{})", order + 1));
	}

	fmt::print("all done\n");
}
