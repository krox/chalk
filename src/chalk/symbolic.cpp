#include "chalk/symbolic.h"

namespace chalk {

std::strong_ordering
ConstantExpression::operator<=>(ConstantExpression const &) const = default;
bool ConstantExpression::operator==(ConstantExpression const &) const = default;
std::strong_ordering
VariableExpression::operator<=>(VariableExpression const &) const = default;
bool VariableExpression::operator==(VariableExpression const &) const = default;
std::strong_ordering
SumExpression::operator<=>(SumExpression const &) const = default;
bool SumExpression::operator==(SumExpression const &) const = default;
std::strong_ordering
ProductExpression::operator<=>(ProductExpression const &) const = default;
bool ProductExpression::operator==(ProductExpression const &) const = default;
std::strong_ordering
PowerExpression::operator<=>(PowerExpression const &other) const
{
	if (auto cmp = base->operator<=>(*other.base); cmp != 0)
		return cmp;
	return exp <=> other.exp;
}
bool PowerExpression::operator==(PowerExpression const &other) const
{
	return exp == other.exp && *base == *other.base;
}
std::strong_ordering
FunctionExpression::operator<=>(FunctionExpression const &) const = default;
bool FunctionExpression::operator==(FunctionExpression const &) const = default;

std::strong_ordering
Expression::operator<=>(Expression const &) const = default;
bool Expression::operator==(Expression const &) const = default;

Expression::Expression(Expression &&other) noexcept
    : data_(std::exchange(other.data_, {}))
{}

Expression &Expression::operator=(Expression &&other) noexcept
{
	if (this != &other)
		data_ = std::exchange(other.data_, ConstantExpression{Rational(0)});
	return *this;
}

Expression::Expression(int value)
    : Expression(ConstantExpression{Rational(value)})
{}

Expression::Expression(Rational const &value)
    : Expression(ConstantExpression{value})
{}

Expression Expression::variable(std::string_view name)
{
	assert(!name.empty());
	return Expression(VariableExpression(std::string(name)));
}

Expression Expression::integer(std::string_view s)
{
	return Expression(Rational(util::parse_int(s)));
}

Expression Expression::sum(std::span<const Expression> summands)
{
	using Term = std::pair<Rational, Expression>;
	std::vector<Term> terms;

	for (auto const &s : summands)
		s.visit(util::overloaded{
		    [&](SumExpression const &ex) {
			    terms.insert(terms.end(), ex.terms.begin(), ex.terms.end());
		    },
		    [&](ProductExpression const &ex) {
			    assert(ex.terms.size() >= 2);
			    if (ex.terms[0].is_constant())
			    {
				    auto coeff = ex.terms[0].as_constant().value;
				    auto tail =
				        Expression::product(std::span(ex.terms).subspan(1));
				    terms.emplace_back(coeff, std::move(tail));
			    }
			    else
				    terms.emplace_back(Rational(1), s);
		    },
		    [&](ConstantExpression const &ex) {
			    terms.emplace_back(ex.value, Expression(1));
		    },
		    [&](auto const &ex) {
			    terms.emplace_back(Rational(1), Expression(ex));
		    }

		});

	// sort terms, combine them, and remove zeros
	std::ranges::sort(terms, {}, &Term::second);
	size_t j = 0;
	for (size_t i = 0; i < terms.size(); ++i)
	{
		if (j != 0 && terms[j - 1].second == terms[i].second)
			terms[j - 1].first += terms[i].first;

		// else new term
		else if (i != j)
			terms[j++] = std::move(terms[i]);
		else
			j++;

		// if coeff = 0 -> remove term
		assert(j);
		if (is_zero(terms[j - 1].first))
			--j;
	}
	terms.resize(j);

	// convention: 'SumExpression' is only used for >= 2 terms.
	if (terms.size() == 0)
		return Expression(0);
	else if (terms.size() == 1)
		return Expression(terms[0].first) * terms[0].second;
	else
		return Expression(SumExpression(std::move(terms)));
}

Expression Expression::product(std::span<const Expression> terms)
{
	// flatten product of products
	std::vector<Expression> r;
	for (auto const &t : terms)
		t.visit(util::overloaded{[&](ProductExpression const &p) {
			                         r.insert(r.end(), p.terms.begin(),
			                                  p.terms.end());
		                         },
		                         [&](auto const &) { r.push_back(t); }});

	auto coeff = Rational(1);
	std::erase_if(r, [&](Expression const &x) {
		if (x.is_constant())
		{
			coeff *= x.as_constant().value;
			return true;
		}
		return false;
	});

	if (is_zero(coeff))
		return Expression(0);
	if (!is_one(coeff))
		r.insert(r.begin(), Expression(coeff));
	if (r.size() == 1)
		return std::move(r[0]);
	return Expression(ProductExpression(std::move(r)));
}

Expression Expression::power(Expression const &base, int exp)
{
	if (exp == 0)
		return Expression(1);
	else if (exp == 1)
		return base;
	else if (base.is_constant())
	{
		// special case: (constant)^n
		return Expression(pow(base.as_constant().value, exp));
	}
	else
		return Expression(PowerExpression(base, exp));
}

Expression Expression::function(std::string_view fun,
                                std::span<const Expression> args)
{
	return Expression(FunctionExpression(
	    std::string(fun), std::vector<Expression>(args.begin(), args.end())));
}

Expression operator+(Expression const &a) { return a; }

Expression operator-(Expression const &a)
{
	return Expression::product(std::array{Expression(-1), a});
}

Expression operator+(Expression const &a, Expression const &b)
{
	return Expression::sum(std::array{a, b});
}

Expression operator-(Expression const &a, Expression const &b)
{
	return Expression::sum(
	    std::array{a, Expression::product(std::array{Expression(-1), b})});
}

Expression operator*(Expression const &a, Expression const &b)
{
	return Expression::product(std::array{a, b});
}

Expression operator/(Expression const &a, Expression const &b)
{
	return Expression::product(std::array{a, Expression::power(b, -1)});
}

Expression pow(Expression const &base, int exp)
{
	return Expression::power(base, exp);
}

bool is_zero(Expression const &a)
{
	return a.visit<bool>(util::overloaded{
	    [](ConstantExpression const &x) { return is_zero(x.value); },
	    [](auto const &) { return false; }});
}
bool is_one(Expression const &a)
{
	return a.visit<bool>(util::overloaded{
	    [](ConstantExpression const &x) { return is_one(x.value); },
	    [](auto const &) { return false; }});
}
bool is_negative(Expression const &a)
{
	return a.visit<bool>(util::overloaded{
	    [](ConstantExpression const &x) { return is_negative(x.value); },
	    [](ProductExpression const &p) {
		    assert(p.terms.size() >= 2);
		    if (p.terms[0].is_constant() &&
		        is_negative(p.terms[0].as_constant().value))
			    return true;
		    return false;
	    },
	    [](auto const &) { return false; }});
}

bool need_parens_product(Expression const &a)
{
	return a.visit<bool>(
	    util::overloaded{[](SumExpression const &) { return true; },
	                     [](auto const &) { return false; }});
}
bool need_parens_power(Expression const &a)
{
	return a.visit<bool>(util::overloaded{
	    [](ConstantExpression const &b) { return need_parens_power(b.value); },
	    [](VariableExpression const &) { return false; },
	    [](PowerExpression const &) { return true; },
	    [](ProductExpression const &) { return true; },
	    [](SumExpression const &) { return true; },
	    [](FunctionExpression const &) { return false; }});
}

namespace {
void print_impl(ConstantExpression const &x,
                util::function_view<void(std::string_view)> sink,
                format_type format, int)
{
	print(x.value, sink, format);
}

void print_impl(VariableExpression const &s,
                util::function_view<void(std::string_view)> sink, format_type,
                int)
{
	sink(s.name);
}

void print_impl(SumExpression const &sum,
                util::function_view<void(std::string_view)> sink,
                format_type format, int depth)
{
	if (format == format_type::debug)
	{
		sink("Sum[");
		bool written = false;

		for (auto const &[coeff, term] : sum.terms)
		{
			if (written)
				sink(", ");
			written = true;
			print(coeff, sink, format);
			sink(", ");
			print(term, sink, format, depth);
		}
		sink("]");
	}
	else
	{
		bool written = false;

		for (auto const &[coeff, term] : sum.terms)
		{
			auto c = coeff;
			if (is_negative(c))
			{
				c = -c;
				if (written)
					sink(" - ");
				else
					sink("-");
			}
			else if (written)
				sink(" + ");
			written = true;

			if (is_one(term))
			{
				print(c, sink, format);
			}
			else
			{
				if (!is_one(c))
				{
					print(c, sink, format);
					sink(" * ");
				}
				if (need_parens_product(term))
					sink("(");
				print(term, sink, format, depth);
				if (need_parens_product(term))
					sink(")");
			}
		}
	}
}

void print_impl(ProductExpression const &prod,
                util::function_view<void(std::string_view)> sink,
                format_type format, int depth)
{
	if (format == format_type::debug)
		sink("Product[");
	bool written = false;

	std::string_view sep = (format == format_type::debug) ? ", " : " * ";

	for (auto const &term : prod.terms)
	{
		if (written)
			sink(sep);
		written = true;

		if (format != format_type::debug && need_parens_product(term))
			sink("(");
		print(term, sink, format, depth);
		if (format != format_type::debug && need_parens_product(term))
			sink(")");
	}
	if (format == format_type::debug)
		sink("]");
}
void print_impl(PowerExpression const &pow,
                util::function_view<void(std::string_view)> sink,
                format_type format, int depth)
{
	if (format == format_type::debug)
	{
		sink("Power[");
		print(*pow.base, sink, format, depth);
		sink(", ");
		sink(std::to_string(pow.exp));
		sink("]");
	}
	else
	{
		if (need_parens_power(*pow.base))
			sink("(");
		print(*pow.base, sink, format, depth);
		if (need_parens_power(*pow.base))
			sink(")");
		sink("^");
		sink(std::to_string(pow.exp));
	}
}

void print_impl(FunctionExpression const &f,
                util::function_view<void(std::string_view)> sink,
                format_type format, int depth)
{
	if (format == format_type::debug)
	{
		sink("Function[");
		sink(f.name);
		sink(", ");
	}
	else
	{
		sink(f.name);
		sink("(");
	}

	if (depth <= 0)
	{
		sink("...");
	}
	else
	{
		bool written = false;
		for (auto const &arg : f.arguments)
		{
			if (written)
				sink(", ");
			written = true;
			print(arg, sink, format, depth - 1);
		}
	}

	if (format == format_type::debug)
		sink("]");
	else
		sink(")");
}

} // namespace

void print(Expression const &a,
           util::function_view<void(std::string_view)> sink, format_type format,
           int depth)
{
	a.visit([&](auto const &x) { print_impl(x, sink, format, depth); });
}

Expression Expression::parse(std::string_view expr)
{
	auto parser = util::Parser(expr);
	auto result = parse(parser);
	parser.expect_end();
	return result;
}

namespace {

Expression parse_sum(util::Parser &parser);
Expression parse_product(util::Parser &parser);
Expression parse_power(util::Parser &parser);
Expression parse_primary(util::Parser &parser);

Expression parse_primary(util::Parser &parser)
{
	if (parser.match('('))
	{
		auto r = parse_sum(parser);
		parser.expect(')');
		return r;
	}
	else if (auto tok = parser.integer(); !tok.empty())
	{
		return Expression::integer(tok);
	}
	/*else if (auto tok = parser.try_match(Tok::floating()); tok)
	{
	    return T(FloatingOctuple(tok->value));
	}*/
	else if (auto ident = parser.ident(); !ident.empty())
	{
		if (parser.match('(')) // function call
		{
			std::vector<Expression> args;
			while (true)
			{
				args.push_back(parse_sum(parser));
				if (parser.match(')'))
					break;
				else
					parser.expect(',');
			}

			return Expression::function(ident, args);
		}
		else // variable / constant
		{
			return Expression::variable(ident);
		}
	}
	else
		throw util::ParseError(fmt::format("unexpected token in input"));
}

Expression parse_power(util::Parser &parser)
{
	bool neg = (bool)parser.match('-');
	auto r = parse_primary(parser);
	if (parser.match('^') || parser.match("**"))
	{
		if (parser.match('-'))
			r = pow(r, -util::parse_int(parser.expect_integer()));
		else
			r = pow(r, util::parse_int(parser.expect_integer()));
	}
	return neg ? -r : r;
}

Expression parse_product(util::Parser &parser)
{
	auto r = parse_power(parser);
	while (true)
		if (parser.match('*'))
			r = r * parse_power(parser);
		else if (parser.match('/'))
			r = r / parse_power(parser);
		else
			return r;
}

Expression parse_sum(util::Parser &parser)
{
	auto r = parse_product(parser);
	while (true)
		if (parser.match('+'))
			r = r + parse_product(parser);
		else if (parser.match('-'))
			r = r - parse_product(parser);
		else
			return r;
}
} // namespace

Expression Expression::parse(util::Parser &parser) { return parse_sum(parser); }

Expression substitute(Expression const &expr,
                      util::hash_map<std::string, Expression> const &table)
{
	return expr.visit<Expression>(util::overloaded{
	    [&](ConstantExpression const &) { return expr; },
	    [&](VariableExpression const &v) {
		    if (auto it = table.find(v.name); it != table.end())
			    return it->second;
		    else
			    return expr;
	    },
	    [&](SumExpression const &s) {
		    std::vector<Expression> terms;
		    for (auto const &[coeff, term] : s.terms)
			    terms.push_back(Expression(coeff) * substitute(term, table));
		    return Expression::sum(terms);
	    },
	    [&](ProductExpression const &p) {
		    std::vector<Expression> terms;
		    for (auto const &term : p.terms)
			    terms.push_back(substitute(term, table));
		    return Expression::product(terms);
	    },
	    [&](PowerExpression const &p) {
		    return pow(substitute(*p.base, table), p.exp);
	    },
	    [&](FunctionExpression const &f) {
		    std::vector<Expression> args;
		    for (auto const &arg : f.arguments)
			    args.push_back(substitute(arg, table));
		    return Expression::function(f.name, args);
	    }});
}

} // namespace chalk
