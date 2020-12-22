// fms_option.t.cpp - test fms option

#include <cassert>
#include <vector>
#include "fms_test.h"
#include "fms_variate.h"
#include "fms_variate_normal.h"

using namespace fms::test;
using namespace fms::variate;

template<class X>
int test_hermite()
{
	for (X x : {X(-2), X(-1), X(0), X(0.1), X(1), X(2)})
	{
		assert(Hermite(0, x) == 1);
		assert(Hermite(1, x) == x);
		assert(Hermite(2, x) == x*x - 1);
		assert(Hermite(3, x) == x*x*x - 3*x);
	}

	return 0;
}
int test_hermite_d = test_hermite<double>();

template<class X>
int test_variate_normal()
{
	//X eps = std::numeric_limits<X>::epsilon();
	//X dx = X(0.001);

	{
		standard_normal<X> n;

		for (size_t k : {0, 1, 2, 3}) {
			for (X x : range<X>(-2, 3, 1)) {
				for (X s : range(X(-0.1), X(0.2), X(0.1))) {
					for (X h : range<X>(2, 4, 1)) {
						h = pow(X(10), -h);
						
						X df = diff([k, x, s, h, &n](X x) { return n.cdf(x, s, k); }, x, h);
						check(df, n.cdf(x, s, k + 1), h);

						X ds = diff([x, s, h, &n](X s) { return n.cdf(x, s); }, s, h);
						check(ds, n.edf(x, s), h);
					}
				}
			}
		}

	}
	{
		standard_normal<X> n;

		assert(n.cumulant(0, 0) == 0); // true for all cumulants
		assert(n.cumulant(0, 1) == 0); // mean
		assert(n.cumulant(0, 2) == 1); // variance
		assert(n.cumulant(0, 3) == 0);
	}

	{
		X mu = 2;
		X sigma = 3;
		affine n(standard_normal<X>{}, mu, sigma);

		assert(n.cumulant(0, 0) == 0); // true for all cumulants
		assert(n.cumulant(0, 1) == mu); // mean
		assert(n.cumulant(0, 2) == sigma * sigma); // variance
		assert(n.cumulant(0, 3) == 0);
	}
	return 0;
}
//int test_variate_normal_f = test_variate_normal<float>();
int test_variate_normal_d = test_variate_normal<double>();
