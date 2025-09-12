#include "chalk/symbolic.h"

namespace chalk {

namespace {

std::shared_ptr<const ExpressionVariant> make_expr(auto &&v)
{
	// custom allocater would be nice here?
	return std::make_shared<const ExpressionVariant>(
	    std::forward<decltype(v)>(v));
}

bool lex_equ(Expression const &a, Expression const &b);
bool lex_cmp(Expression const &a, Expression const &b);

bool lex_equ_impl(Rational const &a, Rational const &b)
{
	return a.num() == b.num() && a.denom() == b.denom();
}

bool lex_equ_impl(std::string const &a, std::string const &b) { return a == b; }

bool lex_equ_impl(SumExpression const &a, SumExpression const &b)
{
	if (a.terms.size() != b.terms.size())
		return false;
	for (size_t i = 0; i < a.terms.size(); ++i)
		if (a.terms[i].first != b.terms[i].first ||
		    !lex_equ(a.terms[i].second, b.terms[i].second))
			return false;
	return true;
}

bool lex_equ_impl(ProductExpression const &a, ProductExpression const &b)
{
	if (a.terms.size() != b.terms.size())
		return false;
	for (size_t i = 0; i < a.terms.size(); ++i)
		if (!lex_equ(a.terms[i], b.terms[i]))
			return false;
	return true;
}

bool lex_equ_impl(PowerExpression const &a, PowerExpression const &b)
{
	return a.exp == b.exp && lex_equ(a.base, b.base);
}

bool lex_equ_impl(FunctionExpression const &a, FunctionExpression const &b)
{
	if (a.name != b.name || a.arguments.size() != b.arguments.size())
		return false;
	for (size_t i = 0; i < a.arguments.size(); ++i)
		if (!lex_equ(a.arguments[i], b.arguments[i]))
			return false;
	return true;
}

bool lex_equ(Expression const &a, Expression const &b)
{
	if (a.data() == b.data())
		return true;
	if (a.data()->index() != b.data()->index())
		return false;
	return a.visit<bool>([&]<class T>(T const &x) {
		auto const &y = std::get<T>(*b.data());
		return lex_equ_impl(x, y);
	});
}

bool lex_cmp_impl(Rational const &a, Rational const &b)
{
	if (a.num() != b.num())
		return a.num() < b.num();
	return a.denom() < b.denom();
}

bool lex_cmp_impl(std::string const &a, std::string const &b) { return a < b; }

bool lex_cmp_impl(SumExpression const &a, SumExpression const &b)
{
	if (a.terms.size() != b.terms.size())
		return a.terms.size() < b.terms.size();
	for (size_t i = 0; i < a.terms.size(); ++i)
	{
		if (lex_cmp(a.terms[i].second, b.terms[i].second))
			return true;
		else if (lex_cmp(b.terms[i].second, a.terms[i].second))
			return false;

		if (a.terms[i].first != b.terms[i].first)
			return a.terms[i].first < b.terms[i].first;
	}
	return false;
}

bool lex_cmp_impl(ProductExpression const &a, ProductExpression const &b)
{
	if (a.terms.size() != b.terms.size())
		return a.terms.size() < b.terms.size();
	for (size_t i = 0; i < a.terms.size(); ++i)
		if (lex_cmp(a.terms[i], b.terms[i]))
			return true;
		else if (lex_cmp(b.terms[i], a.terms[i]))
			return false;
	return false;
}

bool lex_cmp_impl(PowerExpression const &a, PowerExpression const &b)
{
	if (lex_cmp(a.base, b.base))
		return true;
	else if (lex_cmp(b.base, a.base))
		return false;
	return a.exp < b.exp;
}

bool lex_cmp_impl(FunctionExpression const &a, FunctionExpression const &b)
{
	if (a.name != b.name)
		return a.name < b.name;
	if (a.arguments.size() != b.arguments.size())
		return a.arguments.size() < b.arguments.size();
	for (size_t i = 0; i < a.arguments.size(); ++i)
		if (lex_cmp(a.arguments[i], b.arguments[i]))
			return true;
		else if (lex_cmp(b.arguments[i], a.arguments[i]))
			return false;
	return false;
}

bool lex_cmp(Expression const &a, Expression const &b)
{
	if (a.data() == b.data())
		return false;
	if (a.data()->index() != b.data()->index())
		return a.data()->index() < b.data()->index();
	return a.visit<bool>([&]<class T>(T const &x) {
		auto const &y = std::get<T>(*b.data());
		return lex_cmp_impl(x, y);
	});
}

std::pair<Rational, Expression> split_factor(Expression const &a)
{
	using P = std::pair<Rational, Expression>;
	return a.visit<P>(util::overloaded{
	    [](Rational const &x) -> P { return {x, Expression(1)}; },
	    [&](ProductExpression const &x) -> P {
		    assert(x.terms.size() >= 2);
		    if (auto c = std::get_if<Rational>(x.terms[0].data()))
		    {
			    if (x.terms.size() == 2)
				    return P{*c, x.terms[1]};
			    else
				    return P{
				        *c, Expression::product(std::span(x.terms).subspan(1))};
		    }
		    else
			    return P{Rational(1), a};
	    },
	    [&](auto const &) -> P { return P{Rational(1), a}; }});
}

} // namespace

Expression::Expression() : ptr_(make_expr(Rational(0))) {}

Expression::Expression(std::shared_ptr<const ExpressionVariant> p)
    : ptr_(std::move(p))
{
	assert(ptr_);
}

Expression::Expression(int value) : ptr_(make_expr(Rational(value))) {}

Expression::Expression(Rational const &value) : ptr_(make_expr(value)) {}

Expression Expression::variable(std::string_view name)
{
	return Expression(make_expr(std::string(name)));
}

Expression Expression::integer(std::string_view s)
{
	return Expression(Rational(util::parse_int(s)));
}

Expression Expression::sum(std::span<const Expression> terms)
{
	// flatten sum of sums and extract pre-factors
	std::vector<std::pair<Rational, Expression>> r;
	for (auto const &t : terms)
		if (auto s = std::get_if<SumExpression>(&*t.ptr_))
			r.insert(r.end(), s->terms.begin(), s->terms.end());
		else
			r.push_back(split_factor(t));

	// sort terms
	std::sort(r.begin(), r.end(), [](auto const &a, auto const &b) {
		return lex_cmp(a.second, b.second);
	});

	// combine terms and remove zeros
	size_t j = 0;
	for (size_t i = 0; i < r.size(); ++i)
	{
		if (j != 0 && lex_equ(r[j - 1].second, r[i].second))
			r[j - 1].first += r[i].first;

		// else new term
		else if (i != j)
			r[j++] = std::move(r[i]);
		else
			j++;

		// if coeff = 0 -> remove term
		assert(j);
		if (r[j - 1].first == 0)
			--j;
	}
	r.resize(j);

	// convention: 'SumExpression' is only used for >= 2 terms.
	if (r.size() == 0)
		return Expression(0);
	else if (r.size() == 1)
		return Expression(r[0].first) * r[0].second;
	else
		return Expression(make_expr(SumExpression{std::move(r)}));
}

Expression Expression::product(std::span<const Expression> terms)
{
	// flatten product of products
	std::vector<Expression> r;
	for (auto const &t : terms)
		if (auto p = std::get_if<ProductExpression>(&*t.ptr_))
			r.insert(r.end(), p->terms.begin(), p->terms.end());
		else
			r.push_back(t);

	auto coeff = Rational(1);
	std::erase_if(r, [&](auto const &x) {
		if (auto c = std::get_if<Rational>(x.data()))
		{
			coeff *= *c;
			return true;
		}
		return false;
	});

	if (r.size() == 0)
		return Expression(coeff);
	else if (terms.size() == 1 && is_one(coeff))
		return terms[0];
	if (!is_one(coeff))
		r.insert(r.begin(), Expression(coeff));
	return Expression(make_expr(ProductExpression{std::move(r)}));
}

Expression Expression::power(Expression const &base, int exp)
{
	if (exp == 0)
		return Expression(1);
	else if (exp == 1)
		return base;
	else
		return Expression(make_expr(PowerExpression{base, exp}));
}

Expression Expression::function(std::string_view fun,
                                std::span<const Expression> args)
{
	return Expression(make_expr(FunctionExpression{
	    std::string(fun), std::vector<Expression>(args.begin(), args.end())}));
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
	return a.visit<bool>(
	    util::overloaded{[](Rational const &x) { return is_zero(x); },
	                     [](auto const &) { return false; }});
}
bool is_one(Expression const &a)
{
	return a.visit<bool>(
	    util::overloaded{[](Rational const &x) { return x == 1; },
	                     [](auto const &) { return false; }});
}
bool is_negative(Expression const &a)
{
	return a.visit<bool>(util::overloaded{
	    [](Rational const &x) { return is_negative(x); },
	    [](std::string const &) { return false; },
	    [](PowerExpression const &) { return false; },
	    [](ProductExpression const &p) {
		    assert(p.terms.size() >= 2);
		    if (auto c = std::get_if<Rational>(p.terms[0].data()))
			    return is_negative(*c);
		    return false;
	    },
	    [](SumExpression const &) { return false; },
	    [](FunctionExpression const &) { return false; },
	});
}

bool need_parens_product(Expression const &a)
{
	return a.visit<bool>(util::overloaded{
	    [](Rational const &) { return false; },
	    [](std::string const &) { return false; },
	    [](PowerExpression const &) { return false; },
	    [](ProductExpression const &) { return false; },
	    [](SumExpression const &) { return true; },
	    [](FunctionExpression const &) { return false; },
	});
}
bool need_parens_power(Expression const &a)
{
	return a.visit<bool>(util::overloaded{
	    [](Rational const &b) { return need_parens_power(b); },
	    [](std::string const &) { return false; },
	    [](PowerExpression const &) { return true; },
	    [](ProductExpression const &) { return true; },
	    [](SumExpression const &) { return true; },
	    [](FunctionExpression const &) { return false; },
	});
}

namespace {
void print_impl(Rational const &x,
                util::function_view<void(std::string_view)> sink,
                format_type format, int)
{
	print(x, sink, format);
}

void print_impl(std::string const &s,
                util::function_view<void(std::string_view)> sink, format_type,
                int)
{
	sink(s);
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
		print(pow.base, sink, format, depth);
		sink(", ");
		sink(std::to_string(pow.exp));
		sink("]");
	}
	else
	{
		if (need_parens_power(pow.base))
			sink("(");
		print(pow.base, sink, format, depth);
		if (need_parens_power(pow.base))
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

} // namespace chalk
