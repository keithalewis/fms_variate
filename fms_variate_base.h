// fms_variate.h - Random variates.
#pragma once
#include <cmath>

namespace fms::variate {

	inline const char base_doc[] = R"(
A random variable \(X\) is determined by its cumulative distribution function \(F(x) = P(X <= x)\). 
Its cumulant is \(κ(s) = \log E[\exp(s X)]\) and its Esscher transform \(X_s\) is defined by 
\(P(X_s \le x) = P_s(X \le x) = E[1(X \le x) \exp(s X - κ(s))]\).
A variate must implement the derivatives of the cdf of \(X_s\), cumulant, and
the derivative of the cdf of \(X_s\) with respect to \(s\), called sega.
	)";
	template<class X = double, class S = double>
	struct base {
		typedef X xtype;
		typedef S stype;
		virtual ~base() { }
		X cdf(X x, S s, size_t n = 0) const
		{
			return cdf_(x, s, n);
		}
		X edf(X x, S s) const
		{
			return edf_(x, s);
		}
		S cumulant(S s, size_t n) const
		{
			return cumulant_(s, n);
		}
	private:
		virtual X cdf_(X x, S s, size_t n) const = 0;
		virtual X edf_(X x, S s) const = 0;
		virtual S cumulant_(S s, size_t n) const = 0;
	};

	template<class V, class X = typename V::xtype, class S = typename V::stype>
	class affine : public base<X,S> {
		const V& v;
		X mu, sigma;
	public:
		typedef X xtype; // needed???
		typedef S stype;

		affine(const V& v, X mu = 0, X sigma = 1)
			: v(v), mu(mu), sigma(sigma)
		{ }

		X cdf_(X x, S s = 0, size_t n = 0) const override
		{
			return v.cdf((x - mu) / sigma, s, n) / pow(sigma, X(n));
		}

		S edf_(X x, S s) const override
		{
			return sigma * v.edf((x - mu)/sigma, sigma * x);
		}

		S cumulant_(S s, size_t n = 0) const override
		{
			return v.cumulant(sigma * s, n) * pow(sigma, X(n)) + (n == 0 ? mu * s : n == 1 ? mu : 0);
		}
	};

	template<class X, class S>
	inline X cdf(const base<X,S>& v, X x, S s, size_t n)
	{
		return v.cdf(x, s, n);
	}

}
