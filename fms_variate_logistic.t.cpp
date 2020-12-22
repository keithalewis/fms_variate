// fms_variate_logistic.t.cpp - test logistic variate
#include <cassert>
#include "fms_test.h"
#include "fms_variate_logistic.h"

using namespace fms::test;
using namespace fms::variate;

template<class X>
int test_variate_logistic()
{
	{
		X s = 0;
		auto k0 = logistic<X>::cumulant(s, 0);
		assert(k0 == 0);
		auto k1 = logistic<X>::cumulant(s, 1);
		assert(k1 == 0);
		auto k2 = logistic<X>::cumulant(s, 2);
		X a = M_PI * M_PI / 3;
		assert(k2 == a);
		auto k3 = logistic<X>::cumulant(s, 1);
		assert(k3 == 0);
	}
	{
		for (size_t n = 0; n < 4; ++n) {
			for (X s : range(-0.5, 0.5, 0.1)) {
				for (int i : range(2, 5, 1)) {
					X h = pow(0.1, i);
					auto f = [n](X s) { return logistic<X>::cumulant(s, n); };
					X df = diff(f, s, h);
					X f1 = logistic<X>::cumulant(s, n + 1);
					check(df, f1, h);
				}
			}
		}
	}
	{
		logistic<X> n;

		auto xs = range<X>(-2, 1, 1);
		auto ss = range(X(-0.1), X(0.2), X(0.1));
		auto ns = std::valarray<size_t>({ 0,1,2,3 });

		check_range(n, xs, ss, ns, 1e-4);
	}


	return 0;
}
int test_variate_logistic_d = test_variate_logistic<double>();