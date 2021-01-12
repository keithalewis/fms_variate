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
		logistic<X> v;

		X s = 0;
		auto k0 = v.cumulant(s, 0);
		assert(k0 == 0);
		auto k1 = v.cumulant(s, 1);
		assert(k1 == 0);
		auto k2 = v.cumulant(s, 2);
		X a = M_PI * M_PI / 3;
		assert(k2 == a);
		auto k3 = v.cumulant(s, 1);
		assert(k3 == 0);
	}
	{
		logistic<X> v;

		for (size_t n = 0; n < 4; ++n) {
			for (X s : range(-0.5, 0.5, 0.1)) {
				for (int i : range(2, 5, 1)) {
					X h = pow(0.1, i);
					auto f = [n,&v](X s) { return v.cumulant(s, n); };
					X df = diff(f, s, h);
					check(df, f, s, h);
				}
			}
		}
	}
	{
		logistic<X> v;

		auto xs = range<X>(-2, 3, 1);
		auto ss = range(X(-0.1), X(0.2), X(0.1));
		std::valarray<X> hs = { 0.01, 0.001, 0.0001 };

		for (auto s : ss) {
			for (auto n : { 0, 1, 2, 3 }) {
				auto f = [s, n, &v](X x) { return v.cdf(x, s, n); };
				auto df = [s, n, &v](X x) { return v.cdf(x, s, n + 1); };
				check(f, df, xs, hs);
			}
		}
	}


	return 0;
}
int test_variate_logistic_d = test_variate_logistic<double>();