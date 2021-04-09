#include "chalk/bch.h"
#include "chalk/fraction.h"
#include "chalk/free_group.h"
#include "chalk/group_ring.h"
#include "chalk/series.h"
#include "gtest/gtest.h"
#include <fmt/format.h>
using namespace chalk;

using Rat = Fraction<int64_t>;
using Algebra = FreeAlgebra<Rat>;

namespace {

template <int N> Series<Algebra, N> exp(Series<Algebra, N> const &a)
{
	auto inc = Series<Algebra, N>(1);
	auto r = Series<Algebra, N>(0);
	for (int i = 1; !(inc == 0); ++i)
	{
		r += inc;
		inc *= a;
		inc /= i;
	}
	return r;
}

template <int N>
Series<Algebra, N> comm(Series<Algebra, N> const &a,
                        Series<Algebra, N> const &b)
{
	return a * b - b * a;
}

template <int order> void test()
{
	// test BCH formula by explicitly expanding everything

	using R = Series<Algebra, order>;
	Series<Algebra, order> a = R::generator() * Scalar(Algebra::generator(0));
	Series<Algebra, order> b = R::generator() * Scalar(Algebra::generator(1));
	auto error = exp(a) * exp(b) - exp(bch(a, b, order, comm<order>));
	EXPECT_EQ(fmt::format("{}", error), fmt::format("0 + O(x^{})", order + 1));
	// fmt::print("{}\n", exp(a) * exp(b));
}

} // namespace

TEST(FreeGroup, BCH)
{
	test<1>();
	test<2>();
	test<3>();
	test<4>();
	test<5>();
	test<6>();
	test<7>();
	test<8>();
}
