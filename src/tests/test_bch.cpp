#include "chalk/bch.h"
#include "chalk/fraction.h"
#include "chalk/free_group.h"
#include "chalk/group_ring.h"
#include "chalk/polynomial.h"
#include "gtest/gtest.h"
#include <fmt/format.h>
using namespace chalk;

using Rat = Fraction<int64_t>;
using Algebra = FreeAlgebra<Rat>;
using R = Polynomial<Algebra>;

namespace {

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
} // namespace

TEST(FreeGroup, BCH)
{
	// test BCH formula by explicitly expanding everything
	for (int order = 1; order <= 8; ++order)
	{
		auto a = truncate(R::generator() * Algebra::generator(0), order);
		auto b = truncate(R::generator() * Algebra::generator(1), order);
		auto error = exp(a) * exp(b) - exp(bch(a, b, order, comm));
		EXPECT_EQ(fmt::format("{}", error),
		          fmt::format("0 + O(x^{})", order + 1));
	}
}
