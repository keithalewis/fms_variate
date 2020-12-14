// fms_variate.cpp - random variates
#include <cassert>
#include <limits>
#include "fms_variate.h"

using namespace fms::variate;

// constant random variable
template<class X, class S = X>
class constant {
	X c;
public:
	typedef typename X xtype;
	typedef typename S stype;

	constant(X c)
		: c(c)
	{ }

	X cdf(X x, S s, size_t n) const
	{
		if (n == 0) {
			return 1 * (c <= x);
		}
		if (n == 1) {
			return x == c ? std::numeric_limits<X>::infinity() : 0;
		}

		return std::numeric_limits<X>::quiet_NaN();
	}
	S cumulant(S s, size_t n) const
	{
		if (n == 0) {
			return c * s;
		}
		if (n == 1) {
			return c;
		}

		return 0;
	}
	X sega(X, S) const
	{
		return 0;
	}
};

template<class X, class S = X>
int test_constant(X x)
{
	{
		constant one(x);

		assert(mean(one) == x);
		assert(variance(one) == 0);
		assert(cdf(one, x - 1) == 0);
		assert(cdf(one, x + 1) == 1);
		assert(cdf(one, x) == 1);
	}

	return 0;
}
int test_constant_d = test_constant<double>(1);

int main()
{
	return 0;
}