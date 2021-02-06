// fms_variate_constant.h - constant variate
#pragma once
#include <cmath>
#include <limits>

namespace fms::variate {

	// constant random variate
	template<class X, class S = X>
	class constant {
		X c;
	public:
		typedef typename X xtype;
		typedef typename S stype;

		constant(X c)
			: c(c)
		{ }

		// F_s(x) = 1(c <= x) independent of s
		X cdf(X x, S s = 0, unsigned n = 0) const
		{
			if (n == 0) {
				return 1 * (c <= x);
			}
			if (n == 1) {
				// really δ_c(x)
				return x == c ? std::numeric_limits<X>::infinity() : 0;
			}

			// really δ_c^{(n - 1)}(x)
			return std::numeric_limits<X>::quiet_NaN();
		}

		// κ(s) = cs
		S cumulant(S s, unsigned n = 0) const
		{
			if (n == 0) {
				return c * s;
			}
			if (n == 1) {
				return c;
			}

			return 0;
		}

		// F_s(x) = 1(c <= x) does not depend on s
		X edf(S, X) const
		{
			return 0;
		}
	};

}