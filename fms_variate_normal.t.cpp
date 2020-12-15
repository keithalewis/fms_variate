// fms_option.t.cpp - test fms option

#include <cassert>
#include <vector>
#include "fms_variate.h"
#include "fms_variate_normal.h"

using namespace fms;
using namespace fms::variate;
// f'(x) + f'''(x)h^2/6 + ...
template<class F, class X>
inline X diff(const F& f, X x, X h)
{
	return (f(x + h) - f(x - h)) / (2 * h);
}

template<class X>
inline void check(X df, X f1, X h)
{
	X err = fabs(df - f1);
	X rem = h * h;
	X y;
	y = err / rem;

	assert(fabs(df - f1) < 20 * h * h);
}

template<class X>
inline std::vector<X> range(X a, X b, X h)
{
	std::vector<X> r;

	for (X x = a; x < b; x += h) {
		r.push_back(x);
	}

	return r;
}

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

		for (size_t k : {0, 1, 2}) {
			for (X x : range<X>(-2, 3, 1)) {
				for (X s : range(X(-0.1), X(0.2), X(0.1))) {
					for (X h : range<X>(2, 4, 1)) {
						h = pow(X(10), -h);
						X df = diff([k, x, s, h, &n](X x) { return n.cdf(x, s, k); }, x, h);
						check(df, n.cdf(x, s, k + 1), h);
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
int test_variate_normal_f = test_variate_normal<float>();
int test_variate_normal_d = test_variate_normal<double>();
