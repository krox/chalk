#ifndef CHALK_NUMERICS_H
#define CHALK_NUMERICS_H

#include "chalk/polynomial.h"
#include <vector>

namespace chalk {
std::vector<double> roots(Polynomial<double> const &poly);
}

#endif
