// fms_variate.h - Random variates.
// 
// A random variable $X$ is determined by its cumulative distribution function $F(x) = P(X <= x)$. 
// Its cumulant is $κ(s) = \log E[\exp(s X)]$ and its Esscher transform $X_s$ is defined by 
// $P(X_s \le x) = P_s(X \le x) = E[1(X \le x) \exp(s X - κ(s))]$.
// A variate must implment the derivatives of the cdf of $X_s$, cumulant, and sega,
// where sega is the derivative of the cdf of $X_s$ with respect to s.

#pragma once
#include <concepts>
#include <cmath>

namespace fms {

	template<typename V, class X = typename V::xtype, class S = typename V::stype>
	concept variate_concept = requires (V v, X x, S s, size_t n) {
		std::semiregular<V>;
		{ v.cdf(x, s, n) } -> std::convertible_to<X>;
		{ v.cumulant(s, n) } -> std::convertible_to<S>;
		{ v.edf(x, s) } -> std::convertible_to<X>;
	};

	namespace variate {

		// affine transformation mu + sigma X
		template<variate_concept V, class X = typename V::xtype, class S = typename V::stype>
		class affine {
			V v;
			X mu, sigma;
		public:
			typedef X xtype;
			typedef S stype;

			affine(const V& v, X mu = 0, X sigma = 1)
				: v(v), mu(mu), sigma(sigma)
			{ }

			X cdf(X x, S s = 0, size_t n = 0) const
			{
				return v.cdf((x - mu) / sigma, s, n) / pow(sigma, X(n));
			}

			S cumulant(S s, size_t n = 0) const
			{
				return v.cumulant(sigma * s, n) * pow(sigma, X(n)) + (n == 0 ? mu * s : n == 1 ? mu : 0);
			}

			S edf(X x, S s) const
			{
				return sigma * v.edf((x - mu)/sigma, sigma * x);
			}
		};

		template<variate_concept V, class X = typename V::xtype, class S = typename V::stype>
		inline X cdf(const V& v, X x, S s = 0, size_t n = 0)
		{
			return v.cdf(x, s, n);
		}

		template<variate_concept V, class S = typename V::stype>
		inline S cumulant(const V& v, S s, size_t n = 0)
		{
			return v.cumlant(s, n);
		}

		template<variate_concept V, class X = typename V::xtype, class S = typename V::stype>
		inline X edf(const V& v, X x, S s)
		{
			return v.edf(x, s);
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
