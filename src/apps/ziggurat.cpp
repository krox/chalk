#include "chalk/floating.h"
#include "fmt/format.h"
#include "util/gnuplot.h"
#include "util/random.h"
#include "util/stats.h"
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
using std::min, std::max;

using Real = chalk::FloatingOctuple;
using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

Real pdf(Real x) { return exp(-x * x * 0.5); } // normalization irrelevant
Real pdf_diff(Real x) { return -x * exp(-x * x * 0.5); }
constexpr int n = 16; // number of layers
auto low = Real(0.0);
auto high = Real(9.0); // only ~2^-64 of pdf is outside of 9-sigma

int main()
{
	/**
	NOTES:
	This uses a simple Newton algorithm (with some stabilization) to make all
	areas equal size. Therefore eventual convergence is quick, but overall it
	is still quite slow because:
	    * the starting values (equal spacing) are kinda bad
	    * we ask Eigen for a dense matrix decomposition, even though the
	      matrix has a narrow band-structure
	    * we do everything in MPFR, even before we converged at all
	All of these points could be improved, but do we actually care?
	*/

	// first guess: equal spacing (x(0) and x(n) will stay fixed of course)
	Vector x = Vector::Zero(n + 1);
	for (int i = 0; i <= n; ++i)
		x(i) = low + (high - low) * i / n;

	Vector a;
	for (int iter = 0;; ++iter)
	{
		// function values and its derivatives
		Vector y = Vector::Zero(n + 1);
		Vector dy = Vector::Zero(n + 1);
		for (int i = 0; i <= n; ++i)
		{
			y(i) = pdf(x(i));
			dy(i) = pdf_diff(x(i));
		}

		// areas and derivatives thereof
		a = Vector::Zero(n);
		Matrix da = Matrix::Zero(n, n + 1);
		for (int i = 0; i < n; ++i)
		{
			a(i) = x(i + 1) * (y(i) - y(i + 1));
			da(i, i) = x(i + 1) * dy(i);
			da(i, i + 1) = y(i) - y(i + 1) - x(i + 1) * dy(i + 1);
		}

		// after convergence is pretty-good, do one more iteration
		auto error = double((a.maxCoeff() - a.minCoeff()) / a.mean());
		fmt::print("iter = {}, error = {}, first/last are = {}, {}\n", iter,
		           error, double(a(0)), double(a(n - 1)));

		// std::cout << "---- points:" << std::endl << x << std::endl;
		// std::cout << "---- areas:" << std::endl << a << std::endl;

		// differences of adjacent areas and their derivatives
		Matrix v = a.segment(0, n - 1) - a.segment(1, n - 1);
		Matrix M = da.block(0, 1, n - 1, n - 1) - da.block(1, 1, n - 1, n - 1);

		// take a newton step
		Vector x_new = x;
		x_new.segment(1, n - 1) -= M.colPivHouseholderQr().solve(v);

		// for stability, clamp new values proportional to old
		for (int i = 1; i < n; ++i)
			x(i) = min(max(x_new(i), (x(i - 1) + 2 * x(i)) / 3),
			           (x(i + 1) + 2 * x(i)) / 3);

		if (error < 1.0e-12)
			break;
	}

	fmt::print("constexpr double table_x[{}] = {{\n", n + 1);
	for (int i = 0; i <= n; ++i)
		fmt::print("\t{},\n", double(x(i)));
	fmt::print("}};\n");
	fmt::print("constexpr double table_y[{}] = {{\n", n + 1);
	for (int i = 0; i <= n; ++i)
		fmt::print("\t{},\n", double(pdf(x(i))));
	fmt::print("}};\n");
}
