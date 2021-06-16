#ifndef CHALK_BCH_H
#define CHALK_BCH_H

#include <cassert>
#include <vector>

namespace chalk {

namespace bch_detail {

struct BchStage
{
	int i1, i2;
	int za, zb;
};

// tables from http://www.ehu.eus/ccwmuura/bch.html
inline BchStage stages[] = {
    {0, 0, 0, 0},

    // order 1: terms 1-2
    {1, 0, 1, 1},
    {2, 0, 1, 1},

    // order 2: terms 3
    {2, 1, -1, 2},

    // order 3: terms 4-5
    {3, 1, 1, 12},
    {3, 2, -1, 12},

    // order 4: terms 6-8
    {4, 1, 0, 1},
    {4, 2, 1, 24},
    {5, 2, 0, 1},

    // order 5: terms 9-14
    {6, 1, -1, 720},
    {6, 2, -1, 180},
    {7, 2, 1, 180},
    {8, 2, 1, 720},
    {4, 3, -1, 120},
    {5, 3, -1, 360},

    // order 6: 15-23
    {9, 1, 0, 1},
    {9, 2, -1, 1440},
    {10, 2, -1, 360},
    {11, 2, -1, 1440},
    {12, 2, 0, 1},
    {6, 3, 0, 1},
    {7, 3, -1, 240},
    {8, 3, -1, 720},
    {5, 4, 1, 240},

    // order 7: terms 24-41 WRONG
    {15, 1, 1, 30240},
    {15, 2, 1, 5040},
    {16, 2, 1, 3780},
    {17, 2, -1, 3780},
    {18, 2, -1, 5040},
    {19, 2, -1, 30240},
    {9, 3, 1, 2016},
    {10, 3, 23, 15120},
    {11, 3, 1, 5040},
    {12, 3, -1, 10080},
    {13, 3, 1, 1260},
    {14, 3, 1, 5040},
    {6, 4, 1, 5040},
    {7, 4, -1, 10080},
    {8, 4, 1, 1680},
    {6, 5, 13, 15120},
    {7, 5, -1, 1120},
    {8, 5, -1, 5040},

    // order 8: terms 42-71
    {24, 1, 0, 1},
    {24, 2, 1, 60480},
    {25, 2, 1, 10080},
    {26, 2, 23, 120960},
    {27, 2, 1, 10080},
    {28, 2, 1, 60480},
    {29, 2, 0, 1},
    {15, 3, 0, 1},
    {16, 3, 1, 4032},
    {17, 3, 23, 30240},
    {18, 3, 1, 2240},
    {19, 3, 1, 15120},
    {20, 3, 0, 1},
    {21, 3, 1, 2520},
    {22, 3, 1, 10080},
    {9, 4, 0, 1},
    {10, 4, 1, 10080},
    {11, 4, -1, 20160},
    {12, 4, -1, 20160},
    {13, 4, 0, 1},
    {14, 4, -1, 2520},
    {9, 5, 1, 4032},
    {10, 5, 1, 840},
    {11, 5, 1, 1440},
    {12, 5, 1, 12096},
    {13, 5, 1, 1260},
    {14, 5, 1, 10080},
    {7, 6, -1, 10080},
    {8, 6, -13, 30240},
    {8, 7, -1, 3360},
};

} // namespace bch_detail

/**
 * Computes log(exp(a)*exp(b)) using a Baker–Campbell–Hausdorff formula
 *   - comm can be set to  '[](auto& a, auto& b){ return a*b-b*a; }' usually
 *   - order <= 8 is currently implemented. Up to order 20 is available at
 *     www.ehu.eus/ccwmuura/bch.html ( but the required tables are enormous )
 */
template <typename R, typename F>
R bch(R const &a, R const &b, int order, F &&comm)
{
	using namespace bch_detail;
	std::vector<R> hall;
	hall.push_back({});
	hall.push_back(a);
	hall.push_back(b);
	R r = a + b;
	assert(1 <= order && order <= 8);
	// clang-format off
	int max_term = order == 1 ? 2
	             : order == 2 ? 3
	             : order == 3 ? 5
				 : order == 4 ? 8
				 : order == 5 ? 14
				 : order == 6 ? 23
				 : order == 7 ? 41
				 : 71;
	// clang-format on
	for (int i = 3; i <= max_term; ++i)
	{
		auto &x = hall[stages[i].i1];
		auto &y = hall[stages[i].i2];
		hall.push_back(comm(x, y));
		r += hall.back() * stages[i].za / stages[i].zb;
	}
	return r;
}
} // namespace chalk

#endif
