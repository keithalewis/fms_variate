// fms_variate_constant.t.cpp - test constant variate
#include <cassert>
#include <initializer_list>
#include "fms_variate_constant.h"

using namespace fms::variate;

int test_variate_constant()
{
	{
		constant c(1.23);

		for (double s : {-1, 0, 1}) 
		{
			assert(c.cdf(1.22, s) == 0);
			assert(c.cdf(1.23, s) == 1);
			assert(c.cdf(1.24, s) == 1);

			assert(c.cdf(1.22, s, 1) == 0);
			assert(c.cdf(1.23, s, 1) == std::numeric_limits<double>::infinity());
			assert(c.cdf(1.24, s, 1) == 0);

			assert(c.cgf(s, 0) == 1.23*s);
			assert(c.cgf(s, 1) == 1.23);
			assert(c.cgf(s, 2) == 0);
			assert(c.cgf(s, 3) == 0);

			assert(c.sdf(s, 1.22) == 0);
			assert(c.sdf(s, 1.23) == 0);
			assert(c.sdf(s, 1.24) == 0);
		}
	}

	return 0;
}
int test_variate_constant_ = test_variate_constant();
