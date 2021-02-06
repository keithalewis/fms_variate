// fms_variate.h - Random variates.
#define FMS_DOC(name) inline static const char8_t fms_ ## name ## _documentation[]

//inline static const char8_t* fms_variate_documentation 
FMS_DOC(variate) = u8R"xyzyx(
A random variable \(X\) is determined by its cumulative distribution function \(F(x) = P(X <= x)\). 
Its cumulant is \(κ(s) = \log E[\exp(s X)]\) and its Esscher transform \(X_s\) is defined by 
\(F_s(x) = P(X_s \le x) = P_s(X \le x) = E[1(X \le x) \exp(s X - κ(s))]\).
Expected value under \(P_s\) is denoted \(E_s\).
A variate must implment the derivatives of the cdf, cumulant, and edf,
where edf is the derivative of the cdf of \(X_s\) with respect to s.
)xyzyx";

#pragma once
#include <concepts>
#include <cmath>

namespace fms {

	template<typename V, class X = typename V::xtype, class S = typename V::stype>
	concept variate_concept = requires (V v, X x, S s, unsigned n) {
		std::semiregular<V>;
		{ v.cdf(x, s, n) } -> std::convertible_to<X>;
		{ v.cumulant(s, n) } -> std::convertible_to<S>;
		{ v.edf(s, x) } -> std::convertible_to<X>;
	};

	namespace variate {

		// affine transformation mu + sigma X
		template<variate_concept V, class X = typename V::xtype, class S = typename V::stype>
		class affine {
			const V& v;
			X mu, sigma;
		public:
			typedef X xtype;
			typedef S stype;

			FMS_DOC(affine) = u8R"xyzyx(
Given a variate \(X\) construct the variate \(\mu + \sigma X\).
)xyzyx";
			affine(const V& v, X mu = 0, X sigma = 1)
				: v(v), mu(mu), sigma(sigma == 0 ? 1 : sigma)
			{ }

			FMS_DOC(affine_cdf) = u8R"xyzyx(
Return the \(n\)-th derivative of the transformed cumulative distribution function.
)xyzyx";
			X cdf(X x, S s = 0, unsigned n = 0) const
			{
				return v.cdf((x - mu) / sigma, s, n) / pow(sigma, X(n));
			}

			S cumulant(S s, unsigned n = 0) const
			{
				return v.cumulant(sigma * s, n) * pow(sigma, X(n)) + (n == 0 ? mu * s : n == 1 ? mu : 0);
			}

			S edf(S s, X x) const
			{
				return sigma * v.edf(sigma * s, (x - mu)/sigma);
			}
		};

		template<variate_concept V, class X = typename V::xtype, class S = typename V::stype>
		inline X cdf(const V& v, X x, S s = 0, unsigned n = 0)
		{
			return v.cdf(x, s, n);
		}

		template<variate_concept V, class S = typename V::stype>
		inline S cumulant(const V& v, S s, unsigned n = 0)
		{
			return v.cumlant(s, n);
		}

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
