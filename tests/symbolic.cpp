#include "catch2/catch_test_macros.hpp"

#include "chalk/symbolic.h"
#include "util/hash_map.h"
#include <fmt/format.h>
using namespace chalk;
using namespace util;

TEST_CASE("expression parsing - constants", "[expression]")
{
	CHECK(Expression::parse("0").is_constant());
	CHECK(Expression::parse("42").is_constant());
	CHECK(Expression::parse("-17").is_constant());

	CHECK(fmt::format("{}", Expression::parse("0")) == "0");
	CHECK(fmt::format("{}", Expression::parse("42")) == "42");
	CHECK(fmt::format("{}", Expression::parse("-17")) == "-17");
}

TEST_CASE("expression parsing - variables", "[expression]")
{
	auto x = Expression::parse("x");
	auto y = Expression::parse("y");
	auto abc = Expression::parse("abc");

	CHECK(x.is_variable());
	CHECK(y.is_variable());
	CHECK(abc.is_variable());

	CHECK(x.as_variable().name == "x");
	CHECK(y.as_variable().name == "y");
	CHECK(abc.as_variable().name == "abc");

	CHECK(fmt::format("{}", x) == "x");
	CHECK(fmt::format("{}", y) == "y");
	CHECK(fmt::format("{}", abc) == "abc");
}

TEST_CASE("expression parsing - basic arithmetic", "[expression]")
{
	// Addition
	CHECK(fmt::format("{}", Expression::parse("x + y")) == "x + y");
	CHECK(fmt::format("{}", Expression::parse("1 + 2")) == "3");
	CHECK(fmt::format("{}", Expression::parse("x + 0")) == "x");
	CHECK(fmt::format("{}", Expression::parse("0 + x")) == "x");

	// Subtraction
	CHECK(fmt::format("{}", Expression::parse("x - y")) == "x - y");
	CHECK(fmt::format("{}", Expression::parse("5 - 3")) == "2");
	CHECK(fmt::format("{}", Expression::parse("x - 0")) == "x");
	CHECK(fmt::format("{}", Expression::parse("0 - x")) == "-x");

	// Multiplication
	CHECK(fmt::format("{}", Expression::parse("x * y")) == "x * y");
	CHECK(fmt::format("{}", Expression::parse("2 * 3")) == "6");
	CHECK(fmt::format("{}", Expression::parse("x * 1")) == "x");
	CHECK(fmt::format("{}", Expression::parse("1 * x")) == "x");
	CHECK(fmt::format("{}", Expression::parse("x * 0")) == "0");
	CHECK(fmt::format("{}", Expression::parse("0 * x")) == "0");

	// Division
	CHECK(fmt::format("{}", Expression::parse("x / y")) == "x * y^-1");
	CHECK(fmt::format("{}", Expression::parse("6 / 2")) == "3");
	CHECK(fmt::format("{}", Expression::parse("x / 1")) == "x");
}

TEST_CASE("expression parsing - powers", "[expression]")
{
	CHECK(fmt::format("{}", Expression::parse("x^2")) == "x^2");
	CHECK(fmt::format("{}", Expression::parse("x**3")) == "x^3");
	CHECK(fmt::format("{}", Expression::parse("2^3")) == "8");
	CHECK(fmt::format("{}", Expression::parse("x^1")) == "x");
	CHECK(fmt::format("{}", Expression::parse("x^0")) == "1");
	CHECK(fmt::format("{}", Expression::parse("x^-1")) == "x^-1");
	CHECK(fmt::format("{}", Expression::parse("-x^2")) == "-x^2");
}

TEST_CASE("expression parsing - parentheses", "[expression]")
{
	CHECK(fmt::format("{}", Expression::parse("(x)")) == "x");
	CHECK(fmt::format("{}", Expression::parse("(x + y)")) == "x + y");
	CHECK(fmt::format("{}", Expression::parse("(x + y) * z")) == "(x + y) * z");
	CHECK(fmt::format("{}", Expression::parse("x * (y + z)")) == "x * (y + z)");
	CHECK(fmt::format("{}", Expression::parse("(x + y) * (a + b)")) ==
	      "(x + y) * (a + b)");
}

TEST_CASE("expression parsing - function calls", "[expression]")
{
	auto sin_x = Expression::parse("sin(x)");
	auto cos_xy = Expression::parse("cos(x, y)");
	auto f_abc = Expression::parse("f(a, b, c)");

	CHECK(sin_x.is_function());
	CHECK(cos_xy.is_function());
	CHECK(f_abc.is_function());

	CHECK(sin_x.as_function().name == "sin");
	CHECK(sin_x.as_function().arguments.size() == 1);
	CHECK(cos_xy.as_function().name == "cos");
	CHECK(cos_xy.as_function().arguments.size() == 2);
	CHECK(f_abc.as_function().name == "f");
	CHECK(f_abc.as_function().arguments.size() == 3);

	CHECK(fmt::format("{}", sin_x) == "sin(x)");
	CHECK(fmt::format("{}", cos_xy) == "cos(x, y)");
	CHECK(fmt::format("{}", f_abc) == "f(a, b, c)");
}

TEST_CASE("expression simplification - constant folding", "[expression]")
{
	// Addition
	CHECK(is_zero(Expression::parse("0 + 0")));
	CHECK(Expression::parse("1 + 2") == Expression(3));
	CHECK(Expression::parse("5 + (-3)") == Expression(2));

	// Multiplication
	CHECK(is_zero(Expression::parse("0 * 42")));
	CHECK(is_zero(Expression::parse("42 * 0")));
	CHECK(Expression::parse("2 * 3") == Expression(6));
	CHECK(Expression::parse("(-2) * 3") == Expression(-6));

	// Power
	CHECK(Expression::parse("2^3") == Expression(8));
	CHECK(Expression::parse("5^0") == Expression(1));
	CHECK(Expression::parse("0^5") == Expression(0));
}

TEST_CASE("expression simplification - identity elements", "[expression]")
{
	auto x = Expression::variable("x");

	// Additive identity
	CHECK((x + Expression(0)) == x);
	CHECK((Expression(0) + x) == x);

	// Multiplicative identity
	CHECK((x * Expression(1)) == x);
	CHECK((Expression(1) * x) == x);

	// Multiplicative zero
	CHECK(is_zero(x * Expression(0)));
	CHECK(is_zero(Expression(0) * x));

	// Power identities
	CHECK(pow(x, 0) == Expression(1));
	CHECK(pow(x, 1) == x);
}

TEST_CASE("expression simplification - sum normalization", "[expression]")
{
	auto x = Expression::variable("x");
	auto y = Expression::variable("y");

	// Collect like terms
	CHECK((x + x) == Expression(2) * x);
	CHECK((x + x + x) == Expression(3) * x);
	CHECK((x + y + x) == Expression(2) * x + y);

	// Remove zero terms
	CHECK((x + Expression(0)) == x);
	CHECK((Expression(0) + x + Expression(0)) == x);

	// Single term sum becomes just the term
	auto single_term_sum = Expression::sum(std::array{x});
	CHECK(single_term_sum == x);

	// Empty sum becomes zero
	auto empty_sum = Expression::sum(std::span<const Expression>{});
	CHECK(is_zero(empty_sum));
}

TEST_CASE("expression simplification - product normalization", "[expression]")
{
	auto x = Expression::variable("x");
	auto y = Expression::variable("y");

	// Collect constants
	auto expr1 = Expression(2) * Expression(3) * x;
	CHECK(expr1 == Expression(6) * x);

	// Handle zero
	CHECK(is_zero(Expression(0) * x * y));
	CHECK(is_zero(x * Expression(0) * y));

	// Remove ones
	CHECK((Expression(1) * x) == x);
	CHECK((x * Expression(1)) == x);

	// Single term product (this may not create a ProductExpression due to
	// normalization)
	auto single_term_product = Expression::product(std::array{x});
	CHECK(single_term_product == x);
}

TEST_CASE("expression simplification - flattening", "[expression]")
{
	auto x = Expression::variable("x");
	auto y = Expression::variable("y");
	auto z = Expression::variable("z");

	// Sum flattening: (x + y) + z should flatten to x + y + z
	auto sum1 = x + y;
	auto sum2 = sum1 + z;
	CHECK(sum2.is_sum());
	CHECK(sum2.as_sum().terms.size() == 3);

	// Product flattening: (x * y) * z should flatten to x * y * z
	auto prod1 = x * y;
	auto prod2 = prod1 * z;
	CHECK(prod2.is_product());
	CHECK(prod2.as_product().terms.size() == 3);
}

TEST_CASE("expression operations - algebraic properties", "[expression]")
{
	auto x = Expression::variable("x");
	auto y = Expression::variable("y");
	auto z = Expression::variable("z");

	// Unary operations
	CHECK(+x == x);
	CHECK(-(-x) == x);

	// Distributivity (not automatically applied, but operations should work)
	auto left = x * (y + z);
	auto right = x * y + x * z;
	// Note: These won't be equal without explicit distribution rules
	CHECK(left !=
	      right); // Expected since automatic distribution isn't implemented

	// But manual distribution should work
	auto manual_dist = Expression::sum(std::array{x * y, x * z});
	CHECK(manual_dist == (x * y + x * z));
}

TEST_CASE("expression substitute", "[expression]")
{
	auto x = Expression::variable("x");

	// Simple substitution
	util::hash_map<std::string, Expression> table1;
	table1["x"] = Expression(42);

	auto result1 = substitute(x, table1);
	CHECK(result1 == Expression(42));
}

TEST_CASE("expression comparison and equality", "[expression]")
{
	auto x = Expression::variable("x");
	auto y = Expression::variable("y");

	// Syntactic equality
	CHECK(x == x);
	CHECK(y == y);
	CHECK(x != y);

	// Constant equality
	CHECK(Expression(42) == Expression(42));
	CHECK(Expression(42) != Expression(43));

	// Complex expression equality
	auto expr1 = x + y;
	auto expr2 = x + y;
	auto expr3 = y + x; // Note: may not be equal due to term ordering

	CHECK(expr1 == expr2);
	// expr1 and expr3 comparison depends on internal term ordering
}

TEST_CASE("expression formatting", "[expression]")
{
	auto x = Expression::variable("x");
	auto y = Expression::variable("y");

	// Basic formatting
	CHECK(fmt::format("{}", Expression(42)) == "42");
	CHECK(fmt::format("{}", x) == "x");
	CHECK(fmt::format("{}", x + y) == "x + y");
	CHECK(fmt::format("{}", x * y) == "x * y");
	CHECK(fmt::format("{}", pow(x, 2)) == "x^2");
}

TEST_CASE("expression type predicates", "[expression]")
{
	auto x = Expression::variable("x");
	auto sum_expr = x + Expression(1);
	auto prod_expr = x * Expression(2);
	auto pow_expr = pow(x, 2);
	auto func_expr = Expression::function("sin", std::array{x});

	// is_zero and is_one
	CHECK(is_zero(Expression(0)));
	CHECK(!is_zero(Expression(1)));
	CHECK(!is_zero(x));

	CHECK(is_one(Expression(1)));
	CHECK(!is_one(Expression(0)));
	CHECK(!is_one(x));

	// Type checks
	CHECK(Expression(42).is_constant());
	CHECK(x.is_variable());
	CHECK(sum_expr.is_sum());
	CHECK(prod_expr.is_product());
	CHECK(pow_expr.is_power());
	CHECK(func_expr.is_function());

	// Negative checks
	CHECK(!x.is_constant());
	CHECK(!Expression(42).is_variable());
	CHECK(!x.is_sum());
	CHECK(!x.is_product());
	CHECK(!x.is_power());
	CHECK(!x.is_function());
}