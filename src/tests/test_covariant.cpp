#include "chalk/covariant.h"
#include "chalk/fraction.h"
#include "chalk/ideal.h"
#include "chalk/sparse_polynomial.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

int main()
{
	int max_order = 2; // orders of epsilon we want everything in

	/** polynomial ring of coefficients */
	auto ring = PolynomialRing<Fraction<int64_t>, 7>();
	auto seps = ring.generator("seps", max_order * 2);
	auto cA = ring.generator("cA");
	auto eps = seps * seps;
	auto k1 = ring.generator("k1");
	auto k2 = ring.generator("k2");
	auto k3 = ring.generator("k3");
	auto k4 = ring.generator("k4");
	auto k5 = ring.generator("k5");

	using R = SparsePolynomial<Fraction<int64_t>, 7>;

	/** indexed expressions */
	using Cov = Covariant<R>;

	/** extract the coefficient of eps^k */
	auto getEpsOrder = [&](Cov const &a, int k) {
		return map_coefficients(
		    a, [&](R const &poly) { return get_coefficient(poly, 0, 2 * k); });
	};

	/** extract coefficient of cA^k */
	auto getcAOrder = [&](Cov const &a, int k) {
		return map_coefficients(
		    a, [&](R const &poly) { return get_coefficient(poly, 1, k); });
	};

	/** diff(A*exp(-S), i)*exp(S) */
	auto myDiff = [&](Cov const &a, int k) {
		return diff(a, k) - a * Cov("S", k);
	};

	/** the forces */
	auto eta = Cov("eta", 1);
	auto S0 = Cov("S", 1);

	fmt::print("building force term...\n");
	auto f1 = S0 * eps * k1 + eta * seps * k2;
	auto S1 = taylor(S0, -f1, 2 * max_order);
	fmt::print("f1 = {}\n", f1);

	auto f2 = (S0 * k3 + S1 * k4) * eps + eta * seps + S1 * k5 * cA * eps * eps;
	// auto S2 = taylor(S0, -f2, 2 * max_order);
	fmt::print("f2 = {}\n", f2);

	auto f = f2;

	fmt::print("building conditions...\n");
	std::vector<R> cond_list;
	for (int order = 1; order <= max_order; ++order)
	{
		fmt::print("\nat order eps^{}:\n", order);
		auto cond = Cov(0); // create T at order epsilon^order
		for (int k = 1; k <= order * 2; ++k)
		{
			auto ev = Cov(1); // tmp = <f^k> (with indices 1...k)
			for (int i = 1; i <= k; ++i)
				ev *= rename_index(f, 1, i);

			ev = wick_contract(ev, "eta", R(2));
			ev = getEpsOrder(ev, order);

			fmt::print("<f^{}> = {}\n", k, ev);

			auto tmp = ev;
			for (int i = 1; i <= k; ++i)
			{
				tmp = myDiff(tmp, i);
				tmp = tmp / i;
			}
			cond += tmp;
		}
		cond = simplify_lie(cond, "S", cA);

		for (int k = 0; k <= 1; ++k)
		{
			auto tmp = getcAOrder(cond, k);
			fmt::print("\nall conditions at order cA^{}:\n", k);
			// dump_summary(tmp);
			dump(tmp);
			for (auto term : tmp.terms())
				cond_list.push_back(term.coefficient);
		}
	}

	auto ideal = Ideal(cond_list);
	fmt::print("\n(general form)\n");
	ideal.groebner();
	dump(ideal);

	fmt::print("\n(bf / trapezoid)\n");
	cond_list.push_back(k2 - 1);
	auto ideal2 = Ideal(cond_list);
	cond_list.pop_back();
	ideal2.groebner();
	dump(ideal2);

	fmt::print("\n(torrero / minimal step)\n");
	cond_list.push_back(k3);
	auto ideal3 = Ideal(cond_list);
	ideal3.groebner();
	dump(ideal3);
}
