#pragma once

#include "chalk/rings.h"
#include "fmt/format.h"
#include "util/string_id.h"
#include "util/vector.h"
#include <cassert>
#include <stdexcept>

namespace chalk {

// should be moved to a more central location to facilitate interoperability
// between different symbolic types within chalk.
inline util::StringPool symbol_pool;

// Elements of the (non-abelian) free group.
// Generators are the symbols in chalk::symbol_pool,
class FreeProduct
{
  public:
	using Term = std::pair<util::string_id, int16_t>;

  protected:
	util::small_vector<Term, 3> terms_;

	// collect adjacent terms
	void cleanup()
	{
		size_t j = 0;
		for (size_t i = 0; i < terms_.size(); ++i)
		{
			// no particular reason, but empty symbols are probably bogus
			assert(terms_[i].first);

			if (j != 0 && terms_[i].first == terms_[j - 1].first)
				terms_[j - 1].second += terms_[i].second;
			else
				terms_[j++] = terms_[i];

			if (j != 0 && terms_[j - 1].second == 0)
				j--;
		}
		terms_.resize(j);
	}

  public:
	// default = empty word = "one"
	constexpr FreeProduct() = default;
	constexpr FreeProduct(one_t) noexcept {}

	// main constructor
	explicit FreeProduct(std::span<const Term> terms)
	    : terms_(terms.begin(), terms.end())
	{
		cleanup();
	}

	// single generator, e.g: "x^1"
	static FreeProduct generator(util::string_id sym)
	{
		assert(sym);
		FreeProduct r;
		r.terms_.emplace_back(sym, 1);
		return r;
	}

	// ditto
	static FreeProduct generator(std::string_view sym)
	{
		return generator(symbol_pool(sym));
	}

	// read-only access to the terms
	std::span<const Term> terms() const { return terms_; }
};

class FreeAbelianProduct : public FreeProduct
{
	void cleanup()
	{
		std::ranges::sort(terms_, {}, &Term::first);
		FreeProduct::cleanup();
	}

  public:
	// default = empty word = "one"
	constexpr FreeAbelianProduct() = default;
	constexpr FreeAbelianProduct(one_t) noexcept {}

	// main constructor
	explicit FreeAbelianProduct(std::span<const Term> terms)
	    : FreeProduct(terms)
	{
		// TODO: optimize. currently does cleanup twice...
		cleanup();
	}

	// single generator, e.g: "x^1"
	static FreeAbelianProduct generator(util::string_id sym)
	{
		assert(sym);
		FreeAbelianProduct r;
		r.terms_.emplace_back(sym, 1);
		return r;
	}

	// ditto
	static FreeAbelianProduct generator(std::string_view sym)
	{
		return generator(symbol_pool(sym));
	}
};

inline bool operator==(FreeProduct const &a, FreeProduct const &b)
{
	if (a.terms().size() != b.terms().size())
		return false;
	for (size_t i = 0; i < a.terms().size(); ++i)
		if (a.terms()[i] != b.terms()[i])
			return false;
	return true;
}
inline bool operator==(FreeProduct const &a, one_t)
{
	return a.terms().empty();
}

// Lexicographic ordering of words.
// Mathematically not very meaningful, but practically useful.
inline std::strong_ordering operator<=>(FreeProduct const &a,
                                        FreeProduct const &b)
{
	return std::lexicographical_compare_three_way(
	    a.terms().begin(), a.terms().end(), b.terms().begin(), b.terms().end());
}

// TODO: optimize multiplication/division/inverse of Free[Abelian]Product. The
// current "full expression + full cleanup()" is quite wasteful.

inline FreeProduct operator*(FreeProduct const &a, FreeProduct const &b)
{
	util::small_vector<FreeProduct::Term, 63> terms;
	terms.insert(terms.end(), a.terms().begin(), a.terms().end());
	terms.insert(terms.end(), b.terms().begin(), b.terms().end());
	return FreeProduct(terms);
}

inline FreeProduct inverse(FreeProduct const &a)
{
	util::small_vector<FreeProduct::Term, 63> terms;
	for (auto it = a.terms().rbegin(); it != a.terms().rend(); ++it)
		terms.emplace_back(it->first, -it->second);
	return FreeProduct(terms);
}

inline FreeProduct operator/(FreeProduct const &a, FreeProduct const &b)
{
	return a * inverse(b);
}

inline FreeAbelianProduct operator*(FreeAbelianProduct const &a,
                                    FreeAbelianProduct const &b)
{
	util::small_vector<FreeProduct::Term, 63> terms;
	terms.insert(terms.end(), a.terms().begin(), a.terms().end());
	terms.insert(terms.end(), b.terms().begin(), b.terms().end());
	return FreeAbelianProduct(terms);
}

inline FreeAbelianProduct inverse(FreeAbelianProduct const &a)
{
	util::small_vector<FreeProduct::Term, 63> terms;
	for (auto it = a.terms().rbegin(); it != a.terms().rend(); ++it)
		terms.emplace_back(it->first, -it->second);
	return FreeAbelianProduct(terms);
}

inline FreeAbelianProduct operator/(FreeAbelianProduct const &a,
                                    FreeAbelianProduct const &b)
{
	return a * inverse(b);
}

} // namespace chalk

template <> struct fmt::formatter<chalk::FreeProduct>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	auto format(const chalk::FreeProduct &prod, auto &ctx) const
	    -> decltype(ctx.out()) const
	{
		// empty product -> '1'
		if (prod == chalk::one)
			return fmt::format_to(ctx.out(), "1");

		// otherwise print the word (uppercase = inverse)
		auto it = ctx.out();
		bool first = true;
		for (auto [x, k] : prod.terms())
		{
			if (!first)
				it = format_to(it, " * ");
			first = false;
			it = format_to(it, "{}", chalk::symbol_pool(x));
			if (k != 1)
				it = format_to(it, "^{}", k);
		}
		return it;
	}
};

template <>
struct fmt::formatter<chalk::FreeAbelianProduct>
    : public fmt::formatter<chalk::FreeProduct>
{};
