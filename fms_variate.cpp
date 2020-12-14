// fms_variate.cpp - random variates
#include <cassert>
#include <limits>
#include "fms_variate.h"
#include "fms_variate_constant.h"
#include "fms_variate_normal.h"

using namespace fms;
using namespace fms::variate;

// generic variate tests
template<variate_concept M>
int test_variate(const M& m)
{
	{
		assert(m.cumulant(0) == 0);
		assert(m.cdf(0) == cdf(m, 0.));
	}

	return 0;
}
int test_variate_constant = test_variate(constant(1.23));
int test_variate_normal = test_variate(standard_normal<double>{});

template<variate_concept M>
int test_standard_variate(const M& m)
{
	{
		assert(mean(m) == 0);
		assert(variance(m) == 1);
		assert(m.cumulant(0, 0) == 0);
		assert(m.cumulant(0, 1) == 0);
		assert(m.cumulant(0, 2) == 1);
	}

	return 0;
}
int test_standard_variate_normal = test_variate(standard_normal<double>{});

int main()
{
	return 0;
}