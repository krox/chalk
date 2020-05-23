#include "chalk/floating.h"
#include "chalk/numerics.h"
#include "chalk/polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

int main()
{
	{
		using Poly = Polynomial<double>;
		auto x = Poly::generator();

		auto poly = (x * x + 1) * (x - 1) * (x + 1);
		auto zeros = roots(poly);
		assert(zeros.size() == 2);
		assert(std::abs(zeros[0] - (-1)) < 1.0e-12);
		assert(std::abs(zeros[1] - (+1)) < 1.0e-12);
	}

	{
		// alternating harmonic series = ln(2)

		auto exact = Floating::log2();
		auto x = Floating(0);
		std::vector<Floating> xs;
		for (int i = 1; i < 10; ++i)
		{
			if (i % 2)
				x += 1 / Floating(i);
			else
				x -= 1 / Floating(i);
			xs.push_back(x);
			fmt::print("{} (error = {})\n", x, (double)(x - exact));
		}

		while (xs.size() > 1)
		{
			fmt::print("\n");
			for (size_t i = 0; i < xs.size() - 1; ++i)
			{
				xs[i] = 0.5 * (xs[i] + xs[i + 1]);
				fmt::print("{} (error = {})\n", xs[i], (double)(xs[i] - exact));
			}
			xs.pop_back();
		}
	}
}
