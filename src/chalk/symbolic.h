#pragma once

#include "chalk/fraction.h"
#include "chalk/parser.h"
#include "fmt/format.h"
#include "util/functional.h"
#include "util/hash_map.h"
#include "util/json.h" // for 'overloaded'. should be moved eventually
#include "util/memory.h"
#include "util/string.h"
#include <compare>
#include <memory>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace chalk {

class Expression;

struct ConstantExpression
{
	Rational value;

	std::strong_ordering operator<=>(ConstantExpression const &) const;
	bool operator==(ConstantExpression const &) const;
};

struct VariableExpression
{
	std::string name;

	std::strong_ordering operator<=>(VariableExpression const &) const;
	bool operator==(VariableExpression const &) const;
};

struct SumExpression
{
	std::vector<std::pair<Rational, Expression>> terms;

	std::strong_ordering operator<=>(SumExpression const &) const;
	bool operator==(SumExpression const &) const;
};

struct ProductExpression
{
	std::vector<Expression> terms;

	std::strong_ordering operator<=>(ProductExpression const &) const;
	bool operator==(ProductExpression const &) const;
};

struct PowerExpression
{
	util::value_ptr<Expression> base; // not null
	int exp;

	std::strong_ordering operator<=>(PowerExpression const &) const;
	bool operator==(PowerExpression const &) const;

	PowerExpression(Expression const &base, int exp)
	    : base(util::make_value<Expression>(base)), exp(exp)
	{}
};

struct FunctionExpression
{
	std::string name;
	std::vector<Expression> arguments;

	std::strong_ordering operator<=>(FunctionExpression const &) const;
	bool operator==(FunctionExpression const &) const;
};

using ExpressionVariant =
    std::variant<ConstantExpression, VariableExpression, SumExpression,
                 ProductExpression, PowerExpression, FunctionExpression>;

class Expression
{
	ExpressionVariant data_;

	Expression(ExpressionVariant &&data) : data_(std::move(data)) {}

  public:
	// default = ConstantExpression(0)
	Expression() = default;

	Expression(Expression const &) = default;
	Expression &operator=(Expression const &) = default;

	// moved-from object is 'ConstantExpression(0)'
	Expression(Expression &&) noexcept;
	Expression &operator=(Expression &&) noexcept;

	// convenience conversions
	explicit Expression(int value);
	explicit Expression(Rational const &value);

	// factory functions
	// note: due to normalization rules, the result might not actually be of
	// that type. E.g.:
	//   sum({0, X}) -> X
	static Expression variable(std::string_view name);
	static Expression integer(std::string_view s);
	static Expression sum(std::span<const Expression> terms);
	static Expression product(std::span<const Expression> terms);
	static Expression power(Expression const &base, int exp);
	static Expression function(std::string_view fun,
	                           std::span<const Expression> args);

	template <class T> bool is() const
	{
		return std::holds_alternative<T>(data_);
	}
	bool is_constant() const { return is<ConstantExpression>(); }
	bool is_variable() const { return is<VariableExpression>(); }
	bool is_sum() const { return is<SumExpression>(); }
	bool is_product() const { return is<ProductExpression>(); }
	bool is_power() const { return is<PowerExpression>(); }
	bool is_function() const { return is<FunctionExpression>(); }

	template <class T> T const &as() const { return std::get<T>(data_); }
	auto const &as_constant() const { return as<ConstantExpression>(); }
	auto const &as_variable() const { return as<VariableExpression>(); }
	auto const &as_sum() const { return as<SumExpression>(); }
	auto const &as_product() const { return as<ProductExpression>(); }
	auto const &as_power() const { return as<PowerExpression>(); }
	auto const &as_function() const { return as<FunctionExpression>(); }

	// throws if not the full string is consumed
	static Expression parse(std::string_view expr);

	// consumes as much as possible, possibly leaving trailing tokens
	static Expression parse(util::Parser &parser);

	// this defines an arbitrary syntactic order on expressions, nothing
	// mathematically meaningful
	std::strong_ordering operator<=>(Expression const &other) const;

	// syntactic equality
	// (i.e. sufficient but not necessary for mathematical equivalence)
	bool operator==(Expression const &other) const;

	template <class R = void> R visit(auto &&visitor) const
	{
		return std::visit(visitor, data_);
	}
};

Expression substitute(Expression const &expr,
                      util::hash_map<std::string, Expression> const &table);

// Pattern matching and substitution for expressions
// Returns a substitution table if the pattern matches, std::nullopt otherwise
std::optional<util::hash_map<std::string, Expression>>
match_pattern(Expression const &pattern, Expression const &expr);

// Apply pattern substitution rules to an expression
Expression
apply_rules(Expression const &expr,
            std::vector<std::pair<Expression, Expression>> const &rules);

// A context for symbolic computation that maintains a symbol table and provides
// parsing with automatic variable substitution and "let" expression handling.
class SymbolicContext
{
  private:
	util::hash_map<std::string, Expression> symbol_table_;
	std::vector<std::pair<Expression, Expression>> rules_;

  public:
	SymbolicContext() = default;

	// Parse an expression string and apply variable substitution and pattern
	// rules. If the expression starts with "let", it defines a new variable or
	// rule. Format: "let var = expression" or "let pattern = expression" or
	// just "expression"
	Expression operator()(std::string_view source);

	// Get the current symbol table (read-only access)
	util::hash_map<std::string, Expression> const &symbol_table() const
	{
		return symbol_table_;
	}

	// Get the current rule table (read-only access)
	std::vector<std::pair<Expression, Expression>> const &rules() const
	{
		return rules_;
	}

	// Set a variable in the symbol table
	void set_variable(std::string const &name, Expression const &expr)
	{
		symbol_table_[name] = expr;
	}

	// Add a pattern rule
	void add_rule(Expression const &pattern, Expression const &replacement)
	{
		rules_.emplace_back(pattern, replacement);
	}

	// Clear the symbol table and rules
	void clear()
	{
		symbol_table_.clear();
		rules_.clear();
	}

	// Check if a variable exists in the symbol table
	bool has_variable(std::string const &name) const
	{
		return symbol_table_.contains(name);
	}

	// Execute statements from a file
	void execute_file(std::string const &filename);

  private:
	// Execute statements from a string, separated by newlines or semicolons
	void execute_statements(std::string const &content);
};

// basic math operations
Expression operator+(Expression const &a);
Expression operator-(Expression const &a);
Expression operator+(Expression const &a, Expression const &b);
Expression operator-(Expression const &a, Expression const &b);
Expression operator*(Expression const &a, Expression const &b);
Expression operator/(Expression const &a, Expression const &b);
Expression pow(Expression const &base, int exp);

bool is_zero(Expression const &);
bool is_one(Expression const &);

// helpers for pretty printing
bool is_negative(Expression const &);
bool need_parens_product(Expression const &);
bool need_parens_power(Expression const &);
void print(Expression const &a,
           util::function_view<void(std::string_view)> sink,
           format_type format = {}, int depth = INT_MAX);
} // namespace chalk

template <> struct fmt::formatter<chalk::Expression>
{
  private:
	chalk::format_type type_ = chalk::format_type::text;
	int depth_ = INT_MAX;

  public:
	constexpr auto parse(auto &ctx) -> decltype(ctx.begin())
	{
		auto it = ctx.begin();
		auto end = ctx.end();

		if (it != end && '0' <= *it && *it <= '9')
		{
			depth_ = 0;
			while (it != end && '0' <= *it && *it <= '9')
			{
				depth_ = depth_ * 10 + (*it - '0');
				++it;
			}
		}

		if (it != end)
		{
			if (*it == 't')
			{
				type_ = chalk::format_type::text;
				++it;
			}
			else if (*it == 'd')
			{
				type_ = chalk::format_type::debug;
				++it;
			}
			/*else if (*it == 'l')
			{
			    type_ = format_type::latex;
			    ++it;
			}*/
		}
		// check that the end is reached
		if (it != end && *it != '}')
			throw fmt::format_error("invalid format");
		return it;
	}

	auto format(chalk::Expression const &a,
	            auto &ctx) const -> decltype(ctx.out())
	{
		auto it = ctx.out();
		auto sink = [&it](std::string_view s) {
			for (char c : s)
				*it++ = c;
		};
		print(a, sink, type_, depth_);
		return it;
	}
};
