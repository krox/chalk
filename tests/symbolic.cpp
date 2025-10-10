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

TEST_CASE("symbolic context - pattern matching basics",
          "[expression][context][pattern]")
{
	SymbolicContext ctx;

	// Simple function patterns
	ctx("let sin(0) = 0");
	CHECK(ctx("sin(0)") == Expression(0));

	// Should not affect other expressions
	auto sin_x = ctx("sin(x)");
	CHECK(sin_x.is_function());
	CHECK(sin_x.as_function().name == "sin");

	// Constant patterns
	ctx("let foo(1) = 42");
	CHECK(ctx("foo(1)") == Expression(42));
	CHECK(ctx("foo(2)").is_function()); // Should not match
}

TEST_CASE("symbolic context - placeholder patterns",
          "[expression][context][pattern]")
{
	SymbolicContext ctx;

	// Single placeholder
	ctx("let square(_x) = _x^2");
	CHECK(ctx("square(5)") == Expression(25));
	CHECK(ctx("square(x)") == pow(Expression::variable("x"), 2));

	// Multiple placeholders
	ctx("let add(_a, _b) = _a + _b");
	auto result = ctx("add(3, 4)");
	CHECK(result == Expression(7));

	// Same placeholder used multiple times
	ctx("let double(_x) = _x + _x");
	CHECK(ctx("double(5)") == Expression(10));
}

TEST_CASE("symbolic context - nested patterns",
          "[expression][context][pattern]")
{
	SymbolicContext ctx;

	// Nested function patterns
	ctx("let foo(bar(_x)) = sin(_x)");
	auto result = ctx("foo(bar(0))");
	CHECK(result.is_function());
	CHECK(result.as_function().name == "sin");
	CHECK(result.as_function().arguments[0] == Expression(0));

	// Pattern with arithmetic
	ctx("let y = 42");
	auto simplified = ctx("y + 0");
	CHECK(simplified == Expression(42));
}

TEST_CASE("symbolic context - pattern precedence",
          "[expression][context][pattern]")
{
	SymbolicContext ctx;

	// Rules are applied in the order they were defined
	ctx("let f(_x) = _x"); // Generic pattern
	ctx("let f(0) = 100"); // Specific pattern

	CHECK(ctx("f(0)") ==
	      Expression(0)); // First rule matches (f(_x) = _x with _x=0)
	CHECK(ctx("f(5)") ==
	      Expression(5)); // First rule matches (f(_x) = _x with _x=5)

	// If we want specific patterns to take precedence, define them first
	SymbolicContext ctx2;
	ctx2("let g(0) = 100"); // Specific pattern first
	ctx2("let g(_x) = _x"); // Generic pattern second

	CHECK(ctx2("g(0)") == Expression(100)); // Specific pattern matches
	CHECK(ctx2("g(5)") == Expression(5));   // Generic pattern matches
}

TEST_CASE("symbolic context - complex pattern expressions",
          "[expression][context][pattern]")
{
	SymbolicContext ctx;

	// Pattern with operations
	ctx("let derivative(sin(_x)) = cos(_x)");
	auto result = ctx("derivative(sin(x))");
	CHECK(result.is_function());
	CHECK(result.as_function().name == "cos");

	// Pattern with powers
	ctx("let integrate(_x^2) = _x^3 / 3");
	ctx("let x = x"); // Define symbolic x
	auto integral = ctx("integrate(x^2)");
	// Should be x^3 / 3, which is x^3 * (1/3)^(-1) = x^3 * 3^(-1)
}

TEST_CASE("symbolic context - pattern matching with variables",
          "[expression][context][pattern]")
{
	SymbolicContext ctx;

	// Set up some variables
	ctx("let a = 5");
	ctx("let b = 3");

	// Pattern that should work with variables
	ctx("let max(_x, _x) = _x"); // max(x,x) = x
	CHECK(ctx("max(a, a)") == Expression(5));

	// This should not match because different values
	auto max_ab = ctx("max(a, b)");
	CHECK(max_ab.is_function()); // Should remain as function call
}

TEST_CASE("symbolic context - rule application order",
          "[expression][context][pattern]")
{
	SymbolicContext ctx;

	// Rules should be applied in the order they were defined
	ctx("let f(_x) = first");
	ctx("let f(_x) = second");

	// Should use the first rule that matches
	auto result = ctx("f(anything)");
	CHECK(result.is_variable());
	CHECK(result.as_variable().name == "first");
}

TEST_CASE("symbolic context - mixed variables and patterns",
          "[expression][context][pattern]")
{
	SymbolicContext ctx;

	// Define some variables
	ctx("let pi = 3");
	ctx("let e = 2");

	// Define some patterns
	ctx("let ln(e) = 1"); // This pattern expects literal 'e', not the value 2
	ctx("let sin(pi/2) = 1"); // This pattern expects literal 'pi/2', not 3/2

	// For patterns to work with literal variable names, we need to define them
	// before the variables, or create an expression that evaluates to the
	// literal

	// Test that variable substitution happens before pattern matching
	auto result1 =
	    ctx("ln(e)"); // e becomes 2, so this is ln(2), no pattern match
	CHECK(result1.is_function()); // Should remain as ln(2)

	auto result2 =
	    ctx("sin(pi/2)"); // pi becomes 3, so this is sin(3/2), no pattern match
	CHECK(result2.is_function()); // Should remain as sin(3/2)
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

TEST_CASE("symbolic context - basic functionality", "[expression][context]")
{
	SymbolicContext ctx;

	// Basic expression parsing without variables
	auto result1 = ctx("2 + 3");
	CHECK(result1 == Expression(5));

	// Variable definition with let
	auto result2 = ctx("let x = 42");
	CHECK(result2 == Expression(42));
	CHECK(ctx.has_variable("x"));

	// Using previously defined variable
	auto result3 = ctx("x + 1");
	CHECK(result3 == Expression(43));
}

TEST_CASE("symbolic context - let expressions", "[expression][context]")
{
	SymbolicContext ctx;

	// Define multiple variables
	ctx("let a = 5");
	ctx("let b = 3");
	ctx("let c = a + b");

	CHECK(ctx.has_variable("a"));
	CHECK(ctx.has_variable("b"));
	CHECK(ctx.has_variable("c"));

	// Check values
	CHECK(ctx("a") == Expression(5));
	CHECK(ctx("b") == Expression(3));
	CHECK(ctx("c") == Expression(8));

	// Define variable in terms of other variables
	ctx("let sum = a + b + c");
	CHECK(ctx("sum") == Expression(16));

	// Redefine existing variable
	ctx("let a = 10");
	CHECK(ctx("a") == Expression(10));
	// Note: previously defined variables that depended on 'a' won't update
	// automatically
	CHECK(ctx("c") == Expression(8)); // still 5 + 3, not 10 + 3
}

TEST_CASE("symbolic context - symbolic variables", "[expression][context]")
{
	SymbolicContext ctx;

	// Define symbolic variables
	ctx("let x = x"); // This creates a symbolic variable x
	ctx("let y = y");

	// Use in expressions
	auto poly = ctx("x^2 + 2*x + 1");
	CHECK(poly.is_sum());

	// Substitute symbolic with numeric
	ctx("let x = 3");
	auto result = ctx("x + 1");
	CHECK(result == Expression(4));
}

TEST_CASE("symbolic context - error handling", "[expression][context]")
{
	SymbolicContext ctx;

	// Invalid let syntax - missing identifier
	CHECK_THROWS_AS(ctx("let = 42"), util::ParseError);

	// Invalid let syntax - missing assignment
	CHECK_THROWS_AS(ctx("let x 42"), util::ParseError);

	// Trailing input
	CHECK_THROWS_AS(ctx("42 + 3 extra"), util::ParseError);

	// Invalid expression
	CHECK_THROWS_AS(ctx("2 + + 3"), util::ParseError);
}

TEST_CASE("symbolic context - symbol table access", "[expression][context]")
{
	SymbolicContext ctx;

	// Initially empty
	CHECK(ctx.symbol_table().empty());
	CHECK(!ctx.has_variable("x"));

	// Add variables
	ctx("let x = 42");
	ctx("let y = x + 1");

	// Check symbol table
	auto const &table = ctx.symbol_table();
	CHECK(table.size() == 2);
	CHECK(table.contains("x"));
	CHECK(table.contains("y"));
	CHECK(table.at("x") == Expression(42));
	CHECK(table.at("y") == Expression(43));

	// Manual variable setting
	ctx.set_variable("z", Expression(100));
	CHECK(ctx.has_variable("z"));
	CHECK(ctx("z") == Expression(100));

	// Clear symbol table
	ctx.clear();
	CHECK(ctx.symbol_table().empty());
	CHECK(!ctx.has_variable("x"));
	CHECK(!ctx.has_variable("y"));
	CHECK(!ctx.has_variable("z"));
}

TEST_CASE("symbolic context - function expressions", "[expression][context]")
{
	SymbolicContext ctx;

	// Define variables and use in functions
	ctx("let x = 5");
	auto func_result = ctx("sin(x + 1)");

	CHECK(func_result.is_function());
	CHECK(func_result.as_function().name == "sin");
	CHECK(func_result.as_function().arguments.size() == 1);
	CHECK(func_result.as_function().arguments[0] == Expression(6));

	// Functions with multiple arguments
	ctx("let a = 1");
	ctx("let b = 2");
	auto multi_func = ctx("f(a, b, a + b)");

	CHECK(multi_func.is_function());
	CHECK(multi_func.as_function().name == "f");
	CHECK(multi_func.as_function().arguments.size() == 3);
	CHECK(multi_func.as_function().arguments[0] == Expression(1));
	CHECK(multi_func.as_function().arguments[1] == Expression(2));
	CHECK(multi_func.as_function().arguments[2] == Expression(3));
}

TEST_CASE("symbolic context - calculator workflow", "[expression][context]")
{
	SymbolicContext ctx;

	// Simulate a calculator session
	CHECK(ctx("let pi = 3") == Expression(3));
	CHECK(ctx("let r = 5") == Expression(5));

	// Compute area using previously defined variables
	auto area = ctx("let area = pi * r^2");
	CHECK(area == Expression(75));

	// Use the stored variable
	CHECK(ctx("area") == Expression(75));

	// Symbolic computation
	ctx("let x = x"); // Define symbolic variable
	auto poly = ctx("let poly = x^2 + 2*x + 1");
	CHECK(poly.is_sum());

	// Substitute and recompute
	ctx("let x = 10");
	auto poly_val = ctx("poly"); // This should still be the symbolic polynomial

	// Define new expression using current x value
	auto new_poly = ctx("x^2 + 2*x + 1");
	CHECK(new_poly == Expression(121)); // 100 + 20 + 1
}