#include "chalk/ddouble.h"
#include "chalk/floating.h"
#include "util/random.h"
#include "util/stopwatch.h"
#include "gtest/gtest.h"
#include <fmt/format.h>
using namespace chalk;

namespace {

// prevent optimization while benchmarking
static void escape(void *p) { asm volatile("" : : "g"(p) : "memory"); }
static void clobber() { asm volatile("" : : : "memory"); }

FloatingOctuple to_octuple(ddouble a)
{
	return FloatingOctuple(a.high()) + a.low();
}
ddouble to_ddouble(FloatingOctuple const &a)
{
	double high = (double)a;
	double low = (double)(a - high);
	return ddouble::sum(high, low);
} // namespace

template <typename F> void test_unary(F &&f, double min = 0, double max = 1)
{
	auto rng = util::xoshiro256(0);
	double worst = 0.0;
	ddouble worst_x = ddouble(0.0);
	for (int i = 0; i < 10000; ++i)
	{
		ddouble x = ddouble::random(rng) * (max - min) + min;
		ddouble result = f(x);
		FloatingOctuple exact = f(to_octuple(x));
		double error =
		    std::abs((double)(to_octuple(result) - exact) / (double)exact);
		if (error > worst)
		{
			worst = error;
			worst_x = x;
		}
	}
	// fmt::print("worst x = {}\n", worst_x.high());
	EXPECT_LE(worst, 1e-28);
	// fmt::print("worst relative error = {}\n", worst);
}

template <typename F>
void benchmark_unary(std::string const &label, F &&f, double min = 0,
                     double max = 1)
{
	auto rng = util::xoshiro256(0);
	util::Stopwatch sw;
	sw.start();
	for (int i = 0; i < 100000; ++i)
	{
		ddouble x = ddouble::random(rng) * (max - min) + min;
		for (int j = 0; j < 100; ++j)
		{
			escape(&x);
			ddouble result = f(x);
			escape(&result);
		}
	}
	fmt::print("{}: {:.3f} seconds\n", label, sw.secs());
}

template <typename F> void test_binary(F &&f)
{
	auto rng = util::xoshiro256(0);
	double worst = 0.0;
	for (int i = 0; i < 10000; ++i)
	{
		ddouble x = ddouble::random(rng);
		ddouble y = ddouble::random(rng);
		ddouble result = f(x * 2 - 1, y * 2 - 1);
		FloatingOctuple exact = f(to_octuple(x) * 2 - 1, to_octuple(y) * 2 - 1);
		double error =
		    std::abs((double)(to_octuple(result) - exact) / (double)exact);
		worst = std::max(worst, error);
	}
	EXPECT_LE(worst, 1e-28);
	// fmt::print("worst relative error = {}\n", worst);
}

void print_ddouble(std::string const &name, FloatingOctuple a)
{
	double high = (double)a;
	double low = (double)(a - high);
	fmt::print("static constexpr ddouble {}(){{return ddouble({:a}, {:a});}}\n",
	           name, high, low);

	a = 1.0 / a;
	high = (double)a;
	low = (double)(a - high);
	fmt::print(
	    "static constexpr ddouble inv_{}(){{return ddouble({:a}, {:a});}}\n",
	    name, high, low);
}

TEST(Numerics, ddouble)
{
	test_binary([](auto a, auto b) { return a + b; });
	test_binary([](auto a, auto b) { return a - b; });
	test_binary([](auto a, auto b) { return a * b; });
	test_binary([](auto a, auto b) { return a / b; });

	test_unary([](auto a) { return sqrt(a); });
	test_unary([](auto a) { return cbrt(a); });
	test_unary([](auto a) { return rec_sqrt(a); }, 0.1, 1000000000.0);

	test_unary([](auto a) { return exp(a); }, -5.0, 5.0);
	test_unary([](auto a) { return log(a); }, 1e-10, 1e10);

	test_unary([](auto a) { return sin(a); }, -3, 3);

	/*
	print_ddouble("e", exp(FloatingOctuple(1)));
	print_ddouble("egamma", FloatingOctuple::euler());
	print_ddouble("pi", FloatingOctuple::pi());
	print_ddouble("sqrt2", sqrt(FloatingOctuple(2)));
	print_ddouble("sqrt3", sqrt(FloatingOctuple(3)));
	print_ddouble("ln2", log(FloatingOctuple(2)));
	print_ddouble("ln10", log(FloatingOctuple(10)));
	print_ddouble("log10e", 1.0 / log(FloatingOctuple(10)));
	print_ddouble("log2e", 1.0 / log(FloatingOctuple(2)));
	print_ddouble("phi", (sqrt(FloatingOctuple(5)) + 1) / 2);
	*/

	// ASSERT_EQ(zeros.size(), 2);
	// EXPECT_DOUBLE_EQ(zeros[0], -1.);
	// EXPECT_DOUBLE_EQ(zeros[1], +1.);
}

} // namespace
