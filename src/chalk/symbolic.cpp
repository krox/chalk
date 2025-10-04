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

Expression::Expression(int value) : Expression(Rational(value)) {}

Expression::Expression(Rational const &value)
    : variant(ConstantExpression{value})
{}

Expression Expression::variable(std::string_view name)
{
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
	{
		if (auto p = std::get_if<SumExpression>(&s))
			terms.insert(terms.end(), p->terms.begin(), p->terms.end());
		else if (auto c = std::get_if<ConstantExpression>(&s))
			terms.emplace_back(c->value, Expression(1));
		else if (auto p = std::get_if<ProductExpression>(&s))
		{
			assert(p->terms.size() >= 2);
			if (auto c = std::get_if<ConstantExpression>(&p->terms[0]))
				terms.emplace_back(
				    c->value,
				    Expression::product(std::span(p->terms).subspan(1)));
			else
				terms.emplace_back(Rational(1), Expression(s));
		}
		else
			terms.emplace_back(Rational(1), Expression(s));
	}

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
		if (auto p = std::get_if<ProductExpression>(&t))
			r.insert(r.end(), p->terms.begin(), p->terms.end());
		else
			r.push_back(t);

	auto coeff = Rational(1);
	std::erase_if(r, [&](Expression const &x) {
		if (auto c = std::get_if<ConstantExpression>(&x))
		{
			coeff *= c->value;
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
	else if (auto c = std::get_if<ConstantExpression>(&base))
	{
		// special case: (constant)^n
		return Expression(pow(c->value, exp));
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
	auto v = util::overloaded{
	    [](ConstantExpression const &x) { return is_zero(x.value); },
	    [](auto const &) { return false; }};
	return std::visit(v, a);
}
bool is_one(Expression const &a)
{
	auto v = util::overloaded{
	    [](ConstantExpression const &x) { return is_one(x.value); },
	    [](auto const &) { return false; }};
	return std::visit(v, a);
}
bool is_negative(Expression const &a)
{
	auto v = util::overloaded{
	    [](ConstantExpression const &x) { return is_negative(x.value); },
	    [](ProductExpression const &p) {
		    assert(p.terms.size() >= 2);
		    if (auto c = std::get_if<ConstantExpression>(&p.terms[0]))
			    return is_negative(c->value);
		    return false;
	    },
	    [](auto const &) { return false; }};
	return std::visit(v, a);
}

bool need_parens_product(Expression const &a)
{
	auto v = util::overloaded{[](SumExpression const &) { return true; },
	                          [](auto const &) { return false; }};
	return std::visit(v, a);
}
bool need_parens_power(Expression const &a)
{
	auto v = util::overloaded{
	    [](ConstantExpression const &b) { return need_parens_power(b.value); },
	    [](VariableExpression const &) { return false; },
	    [](PowerExpression const &) { return true; },
	    [](ProductExpression const &) { return true; },
	    [](SumExpression const &) { return true; },
	    [](FunctionExpression const &) { return false; }};
	return std::visit(v, a);
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
	std::visit([&](auto const &x) { print_impl(x, sink, format, depth); }, a);
}

} // namespace chalk
