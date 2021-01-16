// fms_variate_normal.h - normal distribution
#pragma once
#include <cmath>

namespace fms::variate {

	// Hermite polynomials H_0(x) = 1, H_1(x) = x, H_{n+1}(x) = x H_n(x) - n H_{n-1}(x)
	template<class X>
	inline constexpr X Hermite(unsigned n, X x)
	{
		if (n == 0) {
			return X(1);
		}
		if (n == 1) {
			return x;
		}

		return x * Hermite(n - 1, x) - X(n - 1) * Hermite(n - 2, x);
	}


	// Normal mean 0 variance 1
	template<class X = double, class S = X>
	class standard_normal
	{
		static constexpr X M_SQRT2 = X(1.41421356237309504880);
		static constexpr X M_SQRT2PI = X(2.50662827463100050240);
	public:
		typedef X xtype;
		typedef S stype;

		/*
		standard_normal()
		{ }
		*/

		static X cdf(X x, S s = 0, unsigned n = 0)
		{
			X x_ = x - s;

			if (n == 0) {
				return (1 + erf(x_ / X(M_SQRT2))) / 2;
			}

			X phi = exp(-x_ * x_ / X(2)) / X(M_SQRT2PI);

			if (n == 1) {
				return phi;
			}

			// (d/dx)^n phi(x) = (-1)^n phi(x) H_n(x)
			return phi * Hermite(n - 1, x_) * ((n&1) ? 1 : -1);
		}

		// (d/ds) cdf(x, s, 0)
		static X edf(X x, S s)
		{
			return -cdf(x, s, 1);
		}

		static S cumulant(S s, unsigned n = 0)
		{
			if (n == 0) {
				return s * s / 2;
			}
			if (n == 1) {
				return s;
			}
			if (n == 2) {
				return 1;
			}

			return S(0);
		}
		
		/*
		template<unsigned N>
		static S cumulant(S s)
		{
			epsilon<S,N> s_(s);

			return s_ * s_ / 2;
		}
		*/
	};
}