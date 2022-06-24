#if 0

#include "chalk/symbolic.h"

namespace chalk {

namespace {
// hide concrete Expression types in anonymous namespace

class VariableExpression : public BaseExpression
{
	std::string name_;

  public:
	explicit VariableExpression(std::string_view name) : name_(name) {}

	std::string const &name() const { return name_; }
	void to_string(fmt::memory_buffer &buf, int parens) const override
	{
		(void)parens;
		format_to(buf, "{}", name_);
	}
};

class ConstantExpression : public BaseExpression
{
	Rational value_;

  public:
	explicit ConstantExpression(int value) : value_(Rational(value)) {}
	explicit ConstantExpression(Integer value)
	    : value_(Rational(std::move(value)))
	{}
	explicit ConstantExpression(Rational value) : value_(std::move(value)) {}

	Rational const &value() const { return value_; }
	void to_string(fmt::memory_buffer &buf, int parens) const override
	{
		if (parens >= 2 && !(value_.denom() == 1))
			format_to(buf, "({})", value_);
		else
			format_to(buf, "{}", value_);
	}
};

class SumExpression : public BaseExpression
{
  private:
	std::vector<ExpressionPtr> terms_;

  public:
	explicit SumExpression(std::vector<ExpressionPtr> terms)
	    : terms_(std::move(terms))
	{
		// the expression should already be fully simplified at this point
		assert(terms_.size() >= 2);
	}

	std::vector<ExpressionPtr> const &terms() const { return terms_; }
	void to_string(fmt::memory_buffer &buf, int parens) const override
	{
		if (parens >= 1)
			format_to(buf, "(");
		for (size_t i = 0; i < terms_.size(); ++i)
		{
			if (i != 0)
				format_to(buf, " + ");
			terms_[i]->to_string(buf, 1);
		}
		if (parens >= 1)
			format_to(buf, ")");
	}
};

class ProductExpression : public BaseExpression
{
  private:
	std::vector<ExpressionPtr> terms_;

  public:
	explicit ProductExpression(std::vector<ExpressionPtr> terms)
	    : terms_(std::move(terms))
	{
		assert(terms_.size() >= 2);
	}

	std::vector<ExpressionPtr> const &terms() const { return terms_; }
	void to_string(fmt::memory_buffer &buf, int parens) const override
	{
		if (parens >= 2)
			format_to(buf, "(");
		for (size_t i = 0; i < terms_.size(); ++i)
		{
			if (i != 0)
				format_to(buf, "*");
			terms_[i]->to_string(buf, 2);
		}
		if (parens >= 2)
			format_to(buf, ")");
	}
};

class PowerExpression : public BaseExpression
{
  private:
	ExpressionPtr base_, exponent_;

  public:
	explicit PowerExpression(ExpressionPtr base, ExpressionPtr exponent)
	    : base_(std::move(base)), exponent_(std::move(exponent))
	{
		assert(base_);
		assert(exponent_);
		// assert(exponent != 0 && exponent != 1);
		// assert(base != 0);
	}

	ExpressionPtr const &base() const { return base_; }
	ExpressionPtr const &exponent() const { return exponent_; }

	void to_string(fmt::memory_buffer &buf, int parens) const override
	{
		if (parens >= 3)
			format_to(buf, "(");
		base_->to_string(buf, 3);
		format_to(buf, "^");
		exponent_->to_string(buf, 2);
		if (parens >= 3)
			format_to(buf, ")");
	}
};

class FunctionExpression : public BaseExpression
{
  private:
	std::string name_;
	ExpressionPtr argument_;

  public:
	explicit FunctionExpression(std::string_view name, ExpressionPtr argument)
	    : name_(name), argument_(std::move(argument))
	{
		assert(name_.size() >= 1);
		assert(argument_ != nullptr);
	}
	std::string const &name() const { return name_; }
	ExpressionPtr const &argument() const { return argument_; }

	void to_string(fmt::memory_buffer &buf, int parens) const override
	{
		(void)parens;
		format_to(buf, "{}(", name_);
		argument_->to_string(buf, 0);
		format_to(buf, ")");
	}
};

// some shortcuts to make dynamic casts more readable
/*std::shared_ptr<const VariableExpression> as_variable(ExpressionPtr const &a)
{
    return std::dynamic_pointer_cast<const VariableExpression>(a);
}
std::shared_ptr<const ConstantExpression> as_constant(ExpressionPtr const &a)
{
    return std::dynamic_pointer_cast<const ConstantExpression>(a);
}*/
std::shared_ptr<const SumExpression> as_sum(ExpressionPtr const &a)
{
	return std::dynamic_pointer_cast<const SumExpression>(a);
}
std::shared_ptr<const ProductExpression> as_product(ExpressionPtr const &a)
{
	return std::dynamic_pointer_cast<const ProductExpression>(a);
}
/*std::shared_ptr<const PowerExpression> as_power(ExpressionPtr const &a)
{
    return std::dynamic_pointer_cast<const PowerExpression>(a);
}
std::shared_ptr<const FunctionExpression> as_function(ExpressionPtr const &a)
{
    return std::dynamic_pointer_cast<const FunctionExpression>(a);
}*/

} // namespace

ExpressionPtr make_variable(std::string_view name)
{
	return std::make_shared<VariableExpression>(name);
}

ExpressionPtr make_constant(int value)
{
	return std::make_shared<ConstantExpression>(value);
}

ExpressionPtr make_constant(Rational const &value)
{
	return std::make_shared<ConstantExpression>(value);
}

ExpressionPtr make_sum(std::vector<ExpressionPtr> const &terms_raw)
{
	// flatten sum of sum
	std::vector<ExpressionPtr> terms;
	for (auto &term : terms_raw)
	{
		if (auto sum = as_sum(term))
			for (auto &term2 : sum->terms())
				terms.push_back(term2);
		else
			terms.push_back(term);
	}

	// sort terms

	// collect terms and remove zeros

	// extract constants factors from terms that are products themselves

	return {std::make_shared<SumExpression>(std::move(terms))};
}

ExpressionPtr make_product(std::vector<ExpressionPtr> const &terms_raw)
{
	// flatten product of product
	std::vector<ExpressionPtr> terms;
	for (auto &term : terms_raw)
	{
		if (auto prod = as_product(term))
		{
			for (auto &term2 : prod->terms())
				terms.push_back(term2);
		}
		else
			terms.push_back(term);
	}

	return std::make_shared<ProductExpression>(std::move(terms));
}

ExpressionPtr make_power(ExpressionPtr const &base,
                         ExpressionPtr const &exponent)
{
	return std::make_shared<PowerExpression>(base, exponent);
}

ExpressionPtr make_function(std::string_view name,
                            ExpressionPtr const &argument)
{
	return std::make_shared<FunctionExpression>(name, argument);
}

} // namespace chalk

#endif
