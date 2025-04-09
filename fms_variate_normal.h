// fms_variate_normal.h - normal distribution
#pragma once
#include <cmath>
#include <numbers>
#include "fms_variate_base.h"

namespace fms::variate {

	template<class X = double>
	constexpr X Hermite(unsigned n, X x)
	{
		if (n == 0) {
			return 1;
		}
		if (n == 1) {
			return x;
		}

		return x * Hermite(n - 1, x) - (n - 1) * Hermite(n - 2, x);
	}


	// Normal mean 0 variance 1
	template<class X = double, class S = X>
	class standard_normal : public interface<X,S> {
#ifndef M_SQRT2
		static constexpr X M_SQRT2 = X(1.41421356237309504880);
#endif
		static constexpr X M_SQRT2PI = X(2.50662827463100050240);
	public:
		typedef X xtype;
		typedef S stype;

		X cdf_(X x, S s, unsigned n) const override
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
			return phi* Hermite(n - 1, x_)* ((n & 1) ? 1 : -1);
		}

		// (d/ds) cdf(x, s, 0)
		X edf_(S s, X x) const override
		{
			return -cdf_(x, s, 1);
		}

		S cumulant_(S s, unsigned n) const override
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
		
	};
}
