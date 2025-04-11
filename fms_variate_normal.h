// fms_variate_normal.h - normal distribution
#pragma once
#include <cmath>
#include <numbers>
#include "fms_variate_interface.h"

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
	class standard_normal : public interface<X, S> {
#ifndef M_SQRT2
		static constexpr X M_SQRT2 = X(1.41421356237309504880);
#endif
		static constexpr X M_SQRT2PI = X(2.50662827463100050240);
	public:
		typedef X xtype;
		typedef S stype;

		X cdf_(X x, S s) const override
		{
			return 0.5 * std::erfc(-(x - s) / M_SQRT2);
		}
		X pdf_(X x, S s) const override
		{
			X x_ = x - s;
			return exp(-x_ * x_ / X(2)) / X(M_SQRT2PI);
		}
		// (d/ds) cdf(x, s, 0)
		X sdf_(S s, X x) const override
		{
			return -cdf_(x, s);
		}
		X mgf_(S s) const override
		{
			return exp(s * s / 2);
		}
		X cgf_(S s) const override
		{
			return s * s / 2;
		}

	};
}
