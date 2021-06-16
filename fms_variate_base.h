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
/*
		FMS_DOC(cdf) = R"xyzyx(
Returns the \(n\)-th derivative of the Esscher transformed cumulative distrubution.
The <em>Esscher transform</em> of the density function \(f\) of a random variable \(X\) is 
\(f_s(x) = f(x)e^{sX - \kappa(s)}\) whenever the <a href="VARIATE.CUMULANT.html">cumulant</a> 
\(\kappa(s)\) exists.
)xyzyx";
		template<variate_concept V, class X = typename V::xtype, class S = typename V::stype>
		inline X cdf(const V& v, X x, S s = 0, unsigned n = 0)
		{
			return v.cdf(x, s, n);
		}

		FMS_DOC(cumulant) = R"xyzyx(
Return the \(n\)-th derivative of the cumulant \(\kappa(s) = \log(E[e^{sX}]).
)xyzyx";
		template<variate_concept V, class S = typename V::stype>
		inline S cumulant(const V& v, S s, unsigned n = 0)
		{
			return v.cumlant(s, n);
		}

		FMS_DOC(edf) = R"xyzyx(
Return the derivative of the Esscher transformed distribution with respect to \(s\),
\(\partial F_s(x)/\partial s = E[1(X \le x)e^{sX - \kappa(s)}(X - \kappa'(s))].
)xyzyx";
		template<variate_concept V, class X = typename V::xtype, class S = typename V::stype>
		inline X edf(const V& v, S s, X x)
		{
			return v.edf(s, x);
		}

		template<variate_concept V, class X = typename V::xtype>
		inline X mean(const V& v)
		{
			return v.cumulant(0, 1);
		}

		template<variate_concept V, class X = typename V::xtype>
		inline X variance(const V& v)
		{
			return v.cumulant(0, 2);
		}

	}

	template<variate_concept V>
	struct variate_model : public V {
		using V::V;
	};

	//template<class X = double, class S = X>
	//using XXX = variate_model<XXX_impl>;

}
*/
