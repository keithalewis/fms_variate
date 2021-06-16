// fms_variate_constant.h - constant variate
#pragma once
#include <cmath>
#include <limits>
#include "fms_variate_base.h"

namespace fms::variate {

	// constant random variate
	template<class X, class S = X>
	class constant : public base<X,S> {
		X c;
	public:
		typedef X xtype;
		typedef S stype;

		constant(X c)
			: c(c)
		{ }

		// F_s(x) = 1(c <= x) independent of s
		X cdf_(X x, S s = 0, size_t n = 0) const override
		{
			if (n == 0) {
				return X(1) * (c <= x);
			}
			if (n == 1) {
				// really δ_c(x)
				return x == c ? std::numeric_limits<X>::infinity() : 0;
			}

			// really δ_c^{(n - 1)}(x)
			return std::numeric_limits<X>::quiet_NaN();
		}

		// F_s(x) = 1(c <= x) does not depend on s
		X edf_(X, S) const override
		{
			return 0;
		}

		// κ(s) = cs
		S cumulant_(S s, size_t n = 0) const override
		{
			if (n == 0) {
				return c * s;
			}
			if (n == 1) {
				return c;
			}

			return 0;
		}

	};

}
