#include "chalk/numtheory.h"
#include "util/lexer.h"
#include <cassert>
#include <fmt/format.h>
using namespace chalk;

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		fmt::print("usage: factor [options] <number(s)>\n");
		return -1;
	}

	for (auto p : primes(util::parse_int(argv[1])))
		fmt::print("{}\n", p);
}
