// fms_variate_constant.h - constant variate
#pragma once
#include <cmath>

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

		X cdf(X x, S s = 0, size_t n = 0) const
		{
			if (n == 0) {
				return 1 * (c <= x);
			}
			if (n == 1) {
				return x == c ? std::numeric_limits<X>::infinity() : 0;
			}

			return std::numeric_limits<X>::quiet_NaN();
		}

		S cumulant(S s, size_t n = 0) const
		{
			if (n == 0) {
				return c * s;
			}
			if (n == 1) {
				return c;
			}

			return 0;
		}

		X edf(X, S) const
		{
			return 0;
		}
	};

}