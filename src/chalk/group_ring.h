#pragma once

#include "chalk/free_group.h"
#include "chalk/rings.h"
#include "fmt/format.h"
#include <vector>

namespace chalk {

template <typename R, typename G> struct MultinomialTerm
{
	R coefficient;
	G index;

	bool operator==(MultinomialTerm const &other) const
	{
		return coefficient == other.coefficient && index == other.index;
	}

	bool operator<(MultinomialTerm const &other) const
	{
		return index < other.index;
	}

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

template <typename R, typename G> class Multinomial
{
	std::vector<MultinomialTerm<R, G>> terms_;

	void cleanup()
	{
		std::sort(terms_.begin(), terms_.end());

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
	Multinomial() = default;
	explicit Multinomial(std::vector<MultinomialTerm<R, G>> terms)
	    : terms_{std::move(terms)}
	{
		cleanup();
	}
	explicit Multinomial(R const &a) : terms_{MultinomialTerm<R, G>{a, G(1)}}
	{
		cleanup();
	}
	explicit Multinomial(G const &b) : terms_{MultinomialTerm<R, G>{R(1), b}}
	{
		cleanup();
	}
	explicit Multinomial(R const &a, G const &b)
	    : terms_{MultinomialTerm<R, G>{a, b}}
	{
		cleanup();
	}
	explicit Multinomial(int a) : terms_{MultinomialTerm<R, G>{R(a), G(1)}}
	{
		cleanup();
	}

	static Multinomial generator(int k)
	{
		return Multinomial({MultinomialTerm<R, G>{R(1), G::generator(k)}});
	}

	std::vector<MultinomialTerm<R, G>> const &terms() const { return terms_; }

	void operator+=(Multinomial const &b)
	{
		assert(this != &b);
		terms_.reserve(terms_.size() + b.terms_.size());
		for (auto &term : b.terms())
			terms_.push_back(term);
		cleanup();
	}
	void operator-=(Multinomial const &b)
	{
		assert(this != &b);
		terms_.reserve(terms_.size() + b.terms_.size());
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

template <typename R, typename G> bool isNegative(Multinomial<R, G> const &)
{
	return false;
}

template <typename R, typename G>
bool operator==(Multinomial<R, G> const &a, Multinomial<R, G> const &b)
{
	return a.terms() == b.terms();
}

template <typename R, typename G>
bool operator==(Multinomial<R, G> const &a, int b)
{
	if (a.terms().empty())
		return b == 0;
	if (a.terms().size() > 1)
		return false;
	return a.terms()[0].index == 1 && a.terms()[0].coefficient == b;
}

template <typename R, typename G>
Multinomial<R, G> operator-(Multinomial<R, G> const &a)
{
	auto r = a.terms();
	for (auto &term : r)
		term.coefficient = -term.coefficient;
	return Multinomial<R, G>(std::move(r));
}

template <typename R, typename G>
Multinomial<R, G> operator+(Multinomial<R, G> const &a,
                            Multinomial<R, G> const &b)
{
	std::vector<MultinomialTerm<R, G>> r;
	r.reserve(a.terms().size() + b.terms().size());
	for (auto &term : a.terms())
		r.push_back(term);
	for (auto &term : b.terms())
		r.push_back(term);
	return Multinomial<R, G>(std::move(r));
}

template <typename R, typename G>
Multinomial<R, G> operator-(Multinomial<R, G> const &a,
                            Multinomial<R, G> const &b)
{
	std::vector<MultinomialTerm<R, G>> r;
	r.reserve(a.terms().size() + b.terms().size());
	for (auto &term : a.terms())
		r.push_back(term);
	for (auto &term : b.terms())
		r.push_back({-term.coefficient, term.index});
	return Multinomial<R, G>(std::move(r));
}

template <typename R, typename G>
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

template <typename R, typename G>
Multinomial<R, G> operator*(Multinomial<R, G> const &a, R const &b)
{
	auto r = a.terms();
	for (auto &term : r)
		term.coefficient *= b;
	return Multinomial<R, G>(std::move(r));
}
template <typename R, typename G>
Multinomial<R, G> operator/(Multinomial<R, G> const &a, R const &b)
{
	auto r = a.terms();
	for (auto &term : r)
		term.coefficient /= b;
	return Multinomial<R, G>(std::move(r));
}
template <typename R, typename G>
Multinomial<R, G> operator*(R const &a, Multinomial<R, G> const &b)
{
	auto r = b.terms();
	for (auto &term : r)
		term.coefficient = a * term.coefficient;
	return Multinomial<R, G>(std::move(r));
}
template <typename R, typename G>
Multinomial<R, G> operator*(Multinomial<R, G> const &a, int b)
{
	auto r = a.terms();
	for (auto &term : r)
		term.coefficient *= b;
	return Multinomial<R, G>(std::move(r));
}
template <typename R, typename G>
Multinomial<R, G> operator/(Multinomial<R, G> const &a, int b)
{
	auto r = a.terms();
	for (auto &term : r)
		term.coefficient /= b;
	return Multinomial<R, G>(std::move(r));
}
template <typename R, typename G>
void operator*=(Multinomial<R, G> &a, Multinomial<R, G> const &b)
{
	a = a * b;
}

/** apply a function to all coefficients */
template <typename R, typename G, typename F>
Multinomial<R, G> mapCoefficients(F f, Multinomial<R, G> const &a)
{
	// TODO: this function should be able to change the coefficient type

	auto terms = a.terms();
	for (auto &t : terms)
		t.coefficient = f(t.coefficient);
	return Multinomial<R, G>(std::move(terms));
}

/** print in multi-line (hopefully more human-readable for huge sums) */
template <typename R, typename G> void dump(Multinomial<R, G> const &a)
{
	for (auto &term : a.terms())
		fmt::print("{} :: {}\n", term.index, term.coefficient);
}

/** ditto, but shorter */
template <typename R, typename G> void dump_summary(Multinomial<R, G> const &a)
{
	for (auto &term : a.terms())
		fmt::print("{} :: {} + ... ( {} terms total )\n", term.index,
		           term.coefficient.lm(), term.coefficient.terms().size());
}

template <typename R, typename G> struct RingTraits<Multinomial<R, G>>
{
	static bool isZero(Multinomial<R, G> const &a)
	{
		return a.terms().empty();
	}
	static bool isOne(Multinomial<R, G> const &a) { return a == 1; }
	static bool isNegative(Multinomial<R, G> const &) { return false; }

	/** these two are not precisely correct, but at least not wrong */
	static bool needParensProduct(Multinomial<R, G> const &) { return true; }
	static bool needParensPower(Multinomial<R, G> const &) { return true; }
};

template <typename R> using FreeAlgebra = Multinomial<R, FreeProduct>;

} // namespace chalk

template <typename R, typename G>
struct fmt::formatter<chalk::Multinomial<R, G>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Multinomial<R, G> &a, FormatContext &ctx)
	    -> decltype(ctx.out())
	{
		// empty sum -> '0'
		if (a.terms().empty())
			return format_to(ctx.out(), "0");

		// otherwise -> terms with '+' inbetween
		auto it = ctx.out();
		for (size_t i = 0; i < a.terms().size(); ++i)
		{
			if (i != 0)
				it = format_to(it, " + ");
			auto const &term = a.terms()[i];

			if (chalk::isOne(term.coefficient))
				it = format_to(it, "{}", term.index);
			else if (term.index == 1)
				it = format_to(it, "{}", term.coefficient);
			else if (chalk::needParensProduct(term.coefficient))
				it = format_to(it, "({})*{}", a.terms()[i].coefficient,
				               term.index);
			else
				it = format_to(it, "{}*{}", term.coefficient, term.index);
		}
		return it;
	}
};
