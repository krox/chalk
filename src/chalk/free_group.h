#ifndef CHALK_FREE_GROUP_H
#define CHALK_FREE_GROUP_H

#include "fmt/format.h"
#include <stdexcept>
#include <vector>

namespace chalk {

/** elements of the free group <a,b,...> */
class FreeProduct
{
	std::vector<int> word_; // even = generators, odd = inverse

	/* (...,a,A,...) -> (...) */
	void cleanup()
	{
		size_t j = 0;
		for (size_t i = 0; i < word_.size(); ++i)
		{
			if (j != 0 && word_[i] == (word_[j - 1] ^ 1))
				--j;
			else
				word_[j++] = word_[i];
		}
		word_.resize(j);
	}

  public:
	/* constructor */
	FreeProduct() = default;
	explicit FreeProduct(std::vector<int> word) : word_{std::move(word)}
	{
		cleanup();
	}
	explicit FreeProduct(int a)
	{
		if (a != 1)
			throw std::runtime_error("invalid use of FreeProduct(int)");
	}

	/** array-like access */
	std::vector<int> const &word() const { return word_; }
	auto begin() const { return word_.begin(); }
	auto end() const { return word_.end(); }
	int operator[](size_t i) const { return word_.at(i); }
	size_t size() const { return word_.size(); }
	bool empty() const { return word_.empty(); }

	static FreeProduct generator(int i)
	{
		return FreeProduct(std::vector<int>{i * 2});
	}

	bool operator==(FreeProduct const &other) const
	{
		return word_ == other.word_;
	}

	/** only for checking 'x == 1' */
	bool operator==(int b) const
	{
		if (b != 1)
			throw std::runtime_error("invalid use of 'FreeProduct == int'");
		return empty();
	}

	/** arbitrary order */
	bool operator<(FreeProduct const &other) const
	{
		return word_ < other.word_;
	}

	void operator*=(FreeProduct const &b)
	{
		word_.reserve(size() + b.size());
		for (int a : b)
			if (!word_.empty() && word_.back() == (a ^ 1))
				word_.pop_back();
			else
				word_.push_back(a);
	}

	void operator/=(FreeProduct const &b)
	{
		word_.reserve(size() + b.size());
		for (size_t i = 0; i < b.size(); ++i)
			if (!word_.empty() && word_.back() == b[b.size() - 1 - i])
				word_.pop_back();
			else
				word_.push_back(b[b.size() - 1 - i] ^ 1);
	}

	FreeProduct operator*(FreeProduct const &b) const
	{
		FreeProduct r;
		r.word_.reserve(size() + b.size());
		r.word_ = word_;
		r *= b;
		return r;
	}

	FreeProduct operator/(FreeProduct const &b) const
	{
		FreeProduct r;
		r.word_.reserve(size() + b.size());
		r.word_ = word_;
		r /= b;
		return r;
	}
};

inline FreeProduct inverse(FreeProduct const &a)
{
	std::vector<int> r;
	r.resize(a.size());
	for (size_t i = 0; i < r.size(); ++i)
		r[i] = a[a.size() - 1 - i] ^ 1;
	return FreeProduct(std::move(r));
}

} // namespace chalk

template <> struct fmt::formatter<chalk::FreeProduct>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::FreeProduct &prod, FormatContext &ctx)
	    -> decltype(ctx.out())
	{
		// empty product -> '1'
		if (prod.empty())
			return format_to(ctx.out(), "1");

		// otherwise print the word (uppercase = inverse)
		auto it = ctx.out();
		for (int a : prod)
		{
			if (a < 0 || a >= 52)
				throw std::runtime_error("format(FreeProduct) only supported "
				                         "for up to 26 generators");
			it++ = "aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ"[a];
		}
		return it;
	}
};

#endif
