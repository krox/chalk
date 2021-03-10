#include "chalk/floating.h"
#include "util/random.h"
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <fmt/format.h>

using namespace chalk;

namespace {

TEST(Numerics, mpfr_eigen)
{
	// matrix type to use
	using Matrix =
	    Eigen::Matrix<FloatingOctuple, Eigen::Dynamic, Eigen::Dynamic>;

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
	EXPECT_TRUE((M * Minv - Matrix::Identity(n, n)).norm() < 1.e-60);
	EXPECT_TRUE((Minv - B.inverse() * A.inverse()).norm() < 1.e-60);
}

} // namespace
