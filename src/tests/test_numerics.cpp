#include "chalk/floating.h"
#include "chalk/fraction.h"
#include "chalk/integer.h"
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
		using Rat = Fraction<Integer>;
		using Poly = Polynomial<Rat>;
		auto x = Poly::generator();

		assert(fmt::format("{}", jacobi(Rat(0), Rat(0), 0, 0, x)) == "1");
		assert(fmt::format("{}", jacobi(Rat(0), Rat(0), 0, 0, x)) == "x");
		assert(fmt::format("{}", jacobi(Rat(0), Rat(0), 0, 0, x)) ==
		       "3/2*x^2 - 1/2");
		assert(fmt::format("{}", jacobi(Rat(0), Rat(0), 0, 0, x)) ==
		       "5/2*x^3 - 3/2*x");
		assert(fmt::format("{}", jacobi(Rat(0), Rat(0), 0, 0, x)) ==
		       "35/8*x^4 - 15/4*x^2 + 3/8");
		assert(fmt::format("{}", jacobi(Rat(0), Rat(0), 0, 0, x)) ==
		       "63/8*x^5 - 35/4*x^3 + 15/8*x");
		assert(fmt::format("{}", jacobi(Rat(0), Rat(0), 0, 0, x)) ==
		       "231/16*x^6 - 315/16*x^4 + 105/16*x^2 - 5/16");

		fmt::print("{}\n", jacobi<Poly, Rat>(Rat(0), Rat(0), 0, 0, x));
		fmt::print("{}\n", jacobi<Poly, Rat>(Rat(0), Rat(0), 1, 0, x));
		fmt::print("{}\n", jacobi<Poly, Rat>(Rat(0), Rat(0), 2, 0, x));
		fmt::print("{}\n", jacobi<Poly, Rat>(Rat(0), Rat(0), 3, 0, x));
		fmt::print("{}\n", jacobi<Poly, Rat>(Rat(0), Rat(0), 4, 0, x));
		fmt::print("{}\n", jacobi<Poly, Rat>(Rat(0), Rat(0), 5, 0, x));
		fmt::print("{}\n", jacobi<Poly, Rat>(Rat(0), Rat(0), 6, 0, x));
	}

	/*{
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
	}*/
}
