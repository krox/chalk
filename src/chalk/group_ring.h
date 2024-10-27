#pragma once

#include "chalk/free_group.h"
#include "chalk/rings.h"
#include "fmt/format.h"
#include <algorithm>
#include <vector>

namespace chalk {

template <class R, class G> struct MultinomialTerm
{
	R coefficient;
	G index;

	friend bool operator==(MultinomialTerm const &,
	                       MultinomialTerm const &) = default;

	MultinomialTerm operator*(MultinomialTerm const &b) const
	{
		return {coefficient * b.coefficient, index * b.index};
	}

	void operator*=(MultinomialTerm const &b)
	{
		coefficient *= b.coefficient;
		index *= b.index;
	}
};

// If more customization is required, might add a trait-type template parameter
template <class R, class G> class Multinomial
{
  public:
	using Term = MultinomialTerm<R, G>;

  private:
	std::vector<Term> terms_;

	void cleanup()
	{
		std::sort(terms_.begin(), terms_.end(),
		          [](auto &a, auto &b) { return a.index < b.index; });

		size_t j = 0;
		for (size_t i = 0; i < terms_.size(); ++i)
		{
			// if index == previous index -> join
			if (j != 0 && terms_[j - 1].index == terms_[i].index)
				terms_[j - 1].coefficient += terms_[i].coefficient;

			// else new term
			else if (i != j)
				terms_[j++] = std::move(terms_[i]);
			else
				j++;

			// if coeff = 0 -> remove term
			assert(j);
			if (terms_[j - 1].coefficient == 0)
				--j;
		}
		terms_.resize(j);
	}

  public:
	// default constructor initializes to zero (= empty sum)
	Multinomial() = default;

	explicit Multinomial(std::vector<Term> terms) : terms_{std::move(terms)}
	{
		cleanup();
	}

	explicit Multinomial(R const &a) : terms_{Term{a, G(1)}} { cleanup(); }

	explicit Multinomial(G const &b) : terms_{Term{R(1), b}} { cleanup(); }

	explicit Multinomial(R const &a, G const &b)
	    : terms_{MultinomialTerm<R, G>{a, b}}
	{
		cleanup();
	}

	explicit Multinomial(int a) : terms_{Term{R(a), G(1)}} { cleanup(); }

	static Multinomial generator(int k)
	{
		return Multinomial({Term{R(1), G::generator(k)}});
	}

	static Multinomial generator(std::string_view s)
	{
		return Multinomial({Term{R(1), G::generator(s)}});
	}

	// direct access to the terms
	std::vector<Term> const &terms() const { return terms_; }

	// inplace operations that actually gain efficiency by being inplace and
	// modifying '.terms_' directly

	void operator+=(Multinomial const &b)
	{
		// NOTE: simply calling '.reserve(... + ...)' here would not be quite
		//       correct, efficiency-wise
		assert(this != &b);
		for (auto &term : b.terms())
			terms_.push_back(term);
		cleanup();
	}
	void operator-=(Multinomial const &b)
	{
		assert(this != &b);
		for (auto &term : b.terms())
			terms_.push_back({-term.coefficient, term.index});
		cleanup();
	}

	void operator*=(R const b)
	{
		for (auto &term : terms_)
			term.coefficient *= b;
		cleanup();
	}
	void operator/=(R const b)
	{
		for (auto &term : terms_)
			term.coefficient /= b;
		cleanup();
	}
	void operator*=(int b)
	{
		for (auto &term : terms_)
			term.coefficient *= b;
		cleanup();
	}
	void operator/=(int b)
	{
		for (auto &term : terms_)
			term.coefficient /= b;
		cleanup();
	}
};

template <class R, class G>
Multinomial<R, G> operator-(Multinomial<R, G> const &a)
{
	auto r = a.terms();
	for (auto &term : r)
		term.coefficient = -term.coefficient;
	return Multinomial<R, G>(std::move(r));
}

template <class R, class G>
Multinomial<R, G> operator+(Multinomial<R, G> const &a,
                            Multinomial<R, G> const &b)
{
	std::vector<MultinomialTerm<R, G>> r;
	r.reserve(a.terms().size() + b.terms().size());
	r.insert(r.end(), a.terms().begin(), a.terms().end());
	r.insert(r.end(), b.terms().begin(), b.terms().end());
	return Multinomial<R, G>(std::move(r));
}

template <class R, class G>
Multinomial<R, G> operator-(Multinomial<R, G> const &a,
                            Multinomial<R, G> const &b)
{
	std::vector<MultinomialTerm<R, G>> r;
	r.reserve(a.terms().size() + b.terms().size());
	r.insert(r.end(), a.terms().begin(), a.terms().end());
	for (auto &term : b.terms())
		r.push_back({-term.coefficient, term.index});
	return Multinomial<R, G>(std::move(r));
}

template <class R, class G>
Multinomial<R, G> operator*(Multinomial<R, G> const &a,
                            Multinomial<R, G> const &b)
{
	std::vector<MultinomialTerm<R, G>> r;
	r.reserve(a.terms().size() * b.terms().size());
	for (auto &t : a.terms())
		for (auto &s : b.terms())
			r.push_back(t * s);
	return Multinomial<R, G>(std::move(r));
}

// apply a function to all coefficients
template <class R, class G, class F>
Multinomial<R, G> map_coefficients(F f, Multinomial<R, G> const &a)
{
	// TODO: this function should be able to change the coefficient type

	auto terms = a.terms();
	for (auto &t : terms)
		t.coefficient = f(t.coefficient);
	return Multinomial<R, G>(std::move(terms));
}

template <class R, class G>
Multinomial<R, G> operator*(Multinomial<R, G> const &a, R const &b)
{
	return map_coefficients([&b](auto &x) { return x * b; }, a);
}
template <class R, class G>
Multinomial<R, G> operator/(Multinomial<R, G> const &a, R const &b)
{
	return map_coefficients([&b](auto &x) { return x / b; }, a);
}
template <class R, class G>
Multinomial<R, G> operator*(R const &a, Multinomial<R, G> const &b)
{
	return map_coefficients([&a](auto &x) { return a * x; }, b);
}
template <class R, class G>
Multinomial<R, G> operator*(Multinomial<R, G> const &a, int b)
{
	return map_coefficients([b](auto &x) { return x * b; }, a);
}
template <class R, class G>
Multinomial<R, G> operator/(Multinomial<R, G> const &a, int b)
{
	return map_coefficients([b](auto &x) { return x / b; }, a);
}
template <class R, class G>
void operator*=(Multinomial<R, G> &a, Multinomial<R, G> const &b)
{
	a = a * b;
}

template <class R, class G>
bool operator==(Multinomial<R, G> const &a, Multinomial<R, G> const &b)
{
	return a.terms() == b.terms();
}

template <class R, class G> bool operator==(Multinomial<R, G> const &a, int b)
{
	if (a.terms().empty())
		return b == 0;
	if (a.terms().size() > 1)
		return false;
	return a.terms()[0].index == 1 && a.terms()[0].coefficient == b;
}

/** print in multi-line (hopefully more human-readable for huge sums) */
template <class R, class G> void dump(Multinomial<R, G> const &a)
{
	for (auto &term : a.terms())
		fmt::print("{} :: {}\n", term.index, term.coefficient);
}

/** ditto, but shorter */
template <class R, class G> void dump_summary(Multinomial<R, G> const &a)
{
	for (auto &term : a.terms())
		fmt::print("{} :: {} + ... ( {} terms total )\n", term.index,
		           term.coefficient.lm(), term.coefficient.terms().size());
}

template <class R, class G> struct RingTraits<Multinomial<R, G>>
{
	static bool is_zero(Multinomial<R, G> const &a)
	{
		return a.terms().empty();
	}
	static bool is_one(Multinomial<R, G> const &a) { return a == 1; }
	static bool is_negative(Multinomial<R, G> const &) { return false; }

	/** these two are not precisely correct, but at least not wrong */
	static bool need_parens_product(Multinomial<R, G> const &) { return true; }
	static bool need_parens_power(Multinomial<R, G> const &) { return true; }
};

template <typename R> using FreeAlgebra = Multinomial<R, FreeProduct>;

} // namespace chalk

template <class R, class G> struct fmt::formatter<chalk::Multinomial<R, G>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Multinomial<R, G> &a,
	            FormatContext &ctx) const -> decltype(ctx.out())
	{
		// empty sum -> '0'
		if (a.terms().empty())
			return fmt::format_to(ctx.out(), "0");

		// otherwise -> terms with '+' inbetween
		auto it = ctx.out();
		for (size_t i = 0; i < a.terms().size(); ++i)
		{
			auto const &term = a.terms()[i];

			R coeff = term.coefficient;
			if (i != 0)
			{
				if (chalk::is_negative(coeff))
				{
					coeff = -coeff;
					it = fmt::format_to(it, " - ");
				}
				else
					it = fmt::format_to(it, " + ");
			}

			if (chalk::is_one(coeff))
				it = fmt::format_to(it, "{}", term.index);
			else if (term.index == 1)
				it = fmt::format_to(it, "{}", coeff);
			else if (chalk::need_parens_product(coeff))
				it = fmt::format_to(it, "({})*{}", coeff, term.index);
			else
				it = fmt::format_to(it, "{}*{}", coeff, term.index);
		}
		return it;
	}
};
