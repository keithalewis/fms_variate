// fms_variate.h - Random variates.
// A random variable $X$ is determined by its cumulative distribution function $F(x) = P(X <= x)$. 
// Its cumulant is $κ(s) = \log E[\exp(s X)]$ and its Esscher transform $X_s$ is defined by 
// $P(X_s \le x) = P_s(X \le x) = E[1(X \le x) \exp(s X - κ(s))]$.
// A variate must implment the derivatives of the cdf of $X_s$, cumulant, and sega,
// where sega is the derivative of the cdf of $X_s$ with respect to s.
#pragma once
#include <concepts>

namespace fms {

	template<typename M, class X = typename M::xtype, class S = typename M::stype>
	concept variate_concept = requires (M m, X x, S s, size_t n) {
		std::semiregular<M>;
		{ m.cdf(x, s, n) } -> std::convertible_to<X>;
		{ m.cumulant(s, n) } -> std::convertible_to<S>;
		{ m.sega(x, s) } -> std::convertible_to<X>;
	};

	namespace variate {

		template<variate_concept M, class X = typename M::xtype, class S = typename M::stype>
		inline X cdf(const M& m, X x, S s = 0, size_t n = 0)
		{
			return m.cdf(x, s, n);
		}

		template<variate_concept M, class S = typename M::stype>
		inline S cumulant(const M& m, S s, size_t n = 0)
		{
			return m.cumlant(s, n);
		}

		template<variate_concept M, class X = typename M::xtype, class S = typename M::stype>
		inline X sega(const M& m, X x, S s = 0)
		{
			return m.sega(x, s);
		}

		template<variate_concept M, class X = typename M::xtype>
		inline X mean(const M& m)
		{
			return m.cumulant(0, 1);
		}

		template<variate_concept M, class X = typename M::xtype>
		inline X variance(const M& m)
		{
			return m.cumulant(0, 2);
		}

	}

	template<variate_concept M>
	struct variate_model : public M {
		using M::M;
	};

	//template<class X = double, class S = X>
	//using XXX = variate_model<XXX_impl>;

}
