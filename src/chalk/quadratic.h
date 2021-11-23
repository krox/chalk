#ifndef CHALK_QUADRATIC
#define CHALK_QUADRATIC

namespace chalk {

/**
 * numbers of the form a + b√d with a,b,d ∈ R
 *   - R assumed to be abelian (otherwise, inverse() is broken)
 *   - d should be a non-square, we further normalize it to be be 'square-free'
 *     (though depending on R, it might be tricky to define what that means)
 *   - TODO: this implementation is quite inefficient. For example:
 *         - many of the (simple) operators could actually use an 'unchecked'
 *           version of the constructor that skips the 'removeSquareFactor' step
 *         - for R = rational numbers, it would be more efficient to store this
 *           as four integers instead of three rationals
 *         - the op-assign operators should actually operate inplace
 */
template <typename R> class Quadratic
{
	R a_, b_, d_; // d is undefined if b=0

  public:
	// constructrs
	Quadratic() = default;
	explicit Quadratic(int a) : a_(a), b_(0) {}
	explicit Quadratic(R a) : a_(std::move(a)), b_(0) {}
	Quadratic(R a, R b, R d)
	    : a_(std::move(a)), b_(std::move(b)), d_(std::move(d))
	{
		b_ *= removeSquareFactor(d_);
	}

	// attributes
	R const &real() const { return a_; }
	R const &imag() const { return b_; }
	R const &d() const { return d_; }
};

// unary Quadratic

template <typename R> Quadratic<R> operator-(Quadratic<R> const &a)
{
	return {-a.real(), -a.imag(), a.d()};
}
template <typename R> Quadratic<R> conj(Quadratic<R> const &a)
{
	// 'Galois conjugate'. This might or might not be the same as 'complex conj'
	return {a.real(), -a.imag(), a.d()};
}
template <typename R> R norm2(Quadratic<R> const &a)
{
	// in number theory, this is just called "norm", not "squared norm".
	// but our naming is more useful when using Quadratic<Rational> as an
	// exact representation of (a subset of) real(/complex) numbers.
	return a.real() * a.real() - a.d() * a.imag() * a.imag();
}
template <typename R> Quadratic<R> inverse(Quadratic<R> const &a)
{
	return conj(a) * inverse(norm2(a));
}

template <typename R> std::optional<Quadratic<R>> trySqrt(Quadratic<R> const &a)
{
	if (a.imag() == 0)
		return Quadratic<R>{R(0), R(1), a.real()};

	auto n = trySqrt(norm2(a));
	if (!n)
		return {};

	auto x = trySqrt((a.real() + *n) / 2);
	if (!x)
		return {};
	auto y = a.imag() / (*x * 2);
	return Quadratic<R>{*std::move(x), std::move(y), a.d()};
}
template <typename R> Quadratic<R> sqrt(Quadratic<R> const &a)
{
	return trySqrt(a).value();
}

// binary Quadratic <-> scalar

template <typename R> Quadratic<R> operator+(Quadratic<R> const &a, R const &b)
{
	return {a.real() + b, a.imag(), a.d()};
}
template <typename R> Quadratic<R> operator-(Quadratic<R> const &a, R const &b)
{
	return {a.real() - b, a.imag(), a.d()};
}
template <typename R> Quadratic<R> operator*(Quadratic<R> const &a, R const &b)
{
	return {a.real() * b, a.imag() * b, a.d()};
}
template <typename R> Quadratic<R> operator/(Quadratic<R> const &a, R const &b)
{
	return {a.real() / b, a.imag() / b, a.d()};
}
template <typename R> Quadratic<R> operator+(Quadratic<R> const &a, int b)
{
	return {a.real() + b, a.imag(), a.d()};
}
template <typename R> Quadratic<R> operator-(Quadratic<R> const &a, int b)
{
	return {a.real() - b, a.imag(), a.d()};
}
template <typename R> Quadratic<R> operator*(Quadratic<R> const &a, int b)
{
	return {a.real() * b, a.imag() * b, a.d()};
}
template <typename R> Quadratic<R> operator/(Quadratic<R> const &a, int b)
{
	return {a.real() / b, a.imag() / b, a.d()};
}

// binary Quadratic <-> Quadratic

template <typename R>
Quadratic<R> operator+(Quadratic<R> const &a, Quadratic<R> const &b)
{
	if (a.imag() == 0)
		return {a.real() + b.real(), b.imag(), b.d()};
	if (b.imag() == 0)
		return {a.real() + b.real(), a.imag(), a.d()};
	assert(a.d() == b.d());
	return {a.real() + b.real(), a.imag() + b.imag(), a.d()};
}
template <typename R>
Quadratic<R> operator-(Quadratic<R> const &a, Quadratic<R> const &b)
{
	if (a.imag() == 0)
		return {a.real() - b.real(), -b.imag(), b.d()};
	if (b.imag() == 0)
		return {a.real() - b.real(), a.imag(), a.d()};
	assert(a.d() == b.d());
	return {a.real() - b.real(), a.imag() - b.imag(), a.d()};
}
template <typename R>
Quadratic<R> operator*(Quadratic<R> const &a, Quadratic<R> const &b)
{
	if (a.imag() == 0)
		return {a.real() * b.real(), a.real() * b.imag(), b.d()};
	if (b.imag() == 0)
		return {a.real() * b.real(), a.imag() * b.real(), a.d()};
	assert(a.d() == b.d());
	return {a.real() * b.real() + a.d() * a.imag() * b.imag(),
	        a.real() * b.imag() + a.imag() * b.real(), a.d()};
}
template <typename R>
Quadratic<R> operator/(Quadratic<R> const &a, Quadratic<R> const &b)
{
	return a * inverse(b);
}

// op-assigns

template <typename R, typename U>
Quadratic<R> &operator+=(Quadratic<R> &a, U const &b)
{
	a = a + b;
	return a;
}
template <typename R, typename U>
Quadratic<R> &operator-=(Quadratic<R> &a, U const &b)
{
	a = a - b;
	return a;
}
template <typename R, typename U>
Quadratic<R> &operator*=(Quadratic<R> &a, U const &b)
{
	a = a * b;
	return a;
}
template <typename R, typename U>
Quadratic<R> &operator/=(Quadratic<R> &a, U const &b)
{
	a = a / b;
	return a;
}

// comparisons

template <typename R>
bool operator==(Quadratic<R> const &a, Quadratic<R> const &b)
{
	return a.real() == b.real() && a.imag() == b.imag() &&
	       (a.d() == b.d() || a.imag() == 0);
}
template <typename R> bool operator==(Quadratic<R> const &a, R const &b)
{
	return a.imag() == 0 && a.real() == b;
}
template <typename R> bool operator==(Quadratic<R> const &a, int b)
{
	return a.imag() == 0 && a.real() == b;
}
template <typename R, typename U>
bool operator!=(Quadratic<R> const &a, U const &b)
{
	return !(a == b);
}

} // namespace chalk

template <typename R> struct fmt::formatter<chalk::Quadratic<R>>
{
	constexpr auto parse(format_parse_context &ctx) { return ctx.begin(); }

	template <typename FormatContext>
	auto format(const chalk::Quadratic<R> &x, FormatContext &ctx)
	{
		if (x.imag() == 0)
			return format_to(ctx.out(), "{}", x.real());
		else if (x.imag() < 0)
			return format_to(ctx.out(), "{} - {}√{}", x.real(), -x.imag(),
			                 x.d());
		else
			return format_to(ctx.out(), "{} + {}√{}", x.real(), x.imag(),
			                 x.d());
	}
};

#endif
