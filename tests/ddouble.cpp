#include "catch2/catch_test_macros.hpp"

#include "chalk/floating.h"
#include "util/ddouble.h"
#include "util/random.h"
#include "util/stopwatch.h"
#include <Eigen/Dense>
#include <fmt/format.h>
using namespace chalk;
using util::ddouble;

namespace {

// prevent optimization while benchmarking
static void escape(void *p) { asm volatile("" : : "g"(p) : "memory"); }
static void clobber() { asm volatile("" : : : "memory"); }

FloatingOctuple to_octuple(ddouble a)
{
	return FloatingOctuple(a.high()) + a.low();
}
FloatingOctuple to_octuple(double a) { return FloatingOctuple(a); }
ddouble to_ddouble(FloatingOctuple const &a)
{
	double high = (double)a;
	double low = (double)(a - high);
	return ddouble::sum(high, low);
} // namespace

template <typename F>
void test_unary([[maybe_unused]] std::string const &label, F f, double min,
                double max, double eps, bool relative)
{
	auto rng = util::xoshiro256(0);
	double worst = 0.0;
	ddouble worst_x = ddouble(0.0);
	for (int i = 0; i < 10000; ++i)
	{
		ddouble x = ddouble::random(rng) * (max - min) + min;
		ddouble result = f(x);
		FloatingOctuple exact = f(to_octuple(x));
		double error = (double)(to_octuple(result) - exact);

		error = relative ? std::abs(error / (double)exact) : std::abs(error);

		if (error > worst)
		{
			worst = error;
			worst_x = x;
		}
	}
	// fmt::print("worst x = {}\n", worst_x.high());
	CHECK(worst < eps);
	/*fmt::print("{:8} : max error ({}) = {:5.3e} (at x = {})\n", label,
	           relative ? "relative" : "absolute", worst, (double)worst_x);*/
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
		ddouble x = ddouble::random(rng) * 2 - 1;
		ddouble y = ddouble::random(rng) * 2 - 1;
		{
			ddouble result = f(x, y);
			FloatingOctuple exact = f(to_octuple(x), to_octuple(y));
			double error =
			    std::abs((double)(to_octuple(result) - exact) / (double)exact);
			worst = std::max(worst, error);
		}
		{
			ddouble result = f(x, y.high());
			FloatingOctuple exact = f(to_octuple(x), to_octuple(y.high()));
			double error =
			    std::abs((double)(to_octuple(result) - exact) / (double)exact);
			worst = std::max(worst, error);
		}
		{
			ddouble result = f(x.high(), y);
			FloatingOctuple exact = f(to_octuple(x.high()), to_octuple(y));
			double error =
			    std::abs((double)(to_octuple(result) - exact) / (double)exact);
			worst = std::max(worst, error);
		}
	}
	CHECK(worst < 1e-28);
	// fmt::print("worst relative error = {}\n", worst);
}

void print_ddouble(std::string const &name, FloatingOctuple a)
{
	double high = (double)a;
	double low = (double)(a - high);
	fmt::print("CHALK_DDOUBLE_CONST({}, {:a}, {:a}) // {:.4f}...\n", name, high,
	           low, high);

	a = 1.0 / a;
	high = (double)a;
	low = (double)(a - high);
	fmt::print("CHALK_DDOUBLE_CONST(inv_{}, {:a}, {:a}) // {:.4f}...\n", name,
	           high, low, high);
}

TEST_CASE("ddouble arithmetic", "[ddouble][numerics]")
{
	test_binary([](auto a, auto b) { return a + b; });
	test_binary([](auto a, auto b) { return a - b; });
	test_binary([](auto a, auto b) { return a * b; });
	test_binary([](auto a, auto b) { return a / b; });

	test_unary(
	    "sqrt", [](auto a) { return sqrt(a); }, 1e-10, 1e10, 1e-29, true);
	test_unary(
	    "cbrt", [](auto a) { return cbrt(a); }, 1e-10, 1e10, 1e-29, true);
	test_unary(
	    "rec_sqrt", [](auto a) { return rec_sqrt(a); }, 1e-10, 1e10, 1e-29,
	    true);
	test_unary(
	    "x^17", [](auto a) { return pow(a, 17); }, 1e-5, 1e5, 1e-29, true);
	test_unary(
	    "x^-3", [](auto a) { return pow(a, -3); }, 1e-5, 1e5, 1e-29, true);

	test_unary(
	    "exp", [](auto a) { return exp(a); }, -5.0, 5.0, 1e-29, true);
	test_unary(
	    "log", [](auto a) { return log(a); }, 1e-10, 1e10, 1e-29, true);

	test_unary(
	    "sin", [](auto a) { return sin(a); }, -10, 10, 1e-29, false);
	test_unary(
	    "cos", [](auto a) { return cos(a); }, -10, 10, 1e-29, false);

	/*benchmark_unary(
	    "sin", [](auto a) { return sin(a); }, -10, 10);
	benchmark_unary(
	    "cos", [](auto a) { return cos(a); }, -10, 10);*/

	/*print_ddouble("e", exp(FloatingOctuple(1)));
	print_ddouble("egamma", FloatingOctuple::euler());
	print_ddouble("pi", FloatingOctuple::pi());
	print_ddouble("sqrt2", sqrt(FloatingOctuple(2)));
	print_ddouble("sqrt3", sqrt(FloatingOctuple(3)));
	print_ddouble("ln2", log(FloatingOctuple(2)));
	print_ddouble("ln10", log(FloatingOctuple(10)));
	print_ddouble("log10e", 1.0 / log(FloatingOctuple(10)));
	print_ddouble("log2e", 1.0 / log(FloatingOctuple(2)));
	print_ddouble("phi", (sqrt(FloatingOctuple(5)) + 1) / 2);

	print_ddouble("fac3", FloatingOctuple(2) * 3);
	print_ddouble("fac4", FloatingOctuple(2) * 3 * 4);
	print_ddouble("fac5", FloatingOctuple(2) * 3 * 4 * 5);
	print_ddouble("fac6", FloatingOctuple(2) * 3 * 4 * 5 * 6);
	print_ddouble("fac7", FloatingOctuple(2) * 3 * 4 * 5 * 6 * 7);
	print_ddouble("fac8", FloatingOctuple(2) * 3 * 4 * 5 * 6 * 7 * 8);
	print_ddouble("fac9", FloatingOctuple(2) * 3 * 4 * 5 * 6 * 7 * 8 * 9);
	print_ddouble("fac10", FloatingOctuple(2) * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10);
	print_ddouble("fac11",
	              FloatingOctuple(2) * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11);
	print_ddouble("fac12", FloatingOctuple(2) * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 *
	                           11 * 12);
	print_ddouble("fac13", FloatingOctuple(2) * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 *
	                           11 * 12 * 13);
	print_ddouble("fac14", FloatingOctuple(2) * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 *
	                           11 * 12 * 13 * 14);*/

	// ASSERT_EQ(zeros.size(), 2);
	// EXPECT_DOUBLE_EQ(zeros[0], -1.);
	// EXPECT_DOUBLE_EQ(zeros[1], +1.);
}

TEST_CASE("combining ddouble with eigen library", "[ddouble][eigen]")
{
	// matrix type to use
	using Matrix = Eigen::Matrix<ddouble, Eigen::Dynamic, Eigen::Dynamic>;

	// get some random matrices
	auto rng = util::xoshiro256(0);
	auto dist = std::uniform_real_distribution<double>(-1.0, 1.0);
	size_t n = 20;
	Matrix A(n, n);
	Matrix B(n, n);
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
		{
			A(i, j) = dist(rng);
			B(i, j) = dist(rng);
		}

	// test basic linear algebra
	Matrix M = A * B;
	Matrix Minv = M.inverse();
	CHECK((double)(M * Minv - Matrix::Identity(n, n)).norm() < 1.e-28);
	CHECK((double)(Minv - B.inverse() * A.inverse()).norm() < 1.e-28);
}

TEST_CASE("misc stuff", "[ddouble]")
{
	// actually, rounding in the last digit is wrong... (TODO)
	CHECK(fmt::format("{:.25e}", ddouble::pi() * 1.e12) ==
	      "3.1415926535897932384626433e+12");
}

} // namespace
