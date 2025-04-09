// fms_variate_hypergeometric.t.cpp - General hypergeometric function
#include <cassert>
#include <algorithm>
#include "fms_test.h"
#include "fms_sf_hypergeometric.h"

using namespace fms::sf;

template<class X>
inline X abs1(X x) { return std::max<X>(X(1), abs(x)); }

#define IF_DOUBLE(x) std::is_same_v<double, X> ? x
#define IF_FLOAT(x) std::is_same_v<float, X> ? x

template<class X>
int test_hypergeometric()
{
	constexpr X epsilon = std::numeric_limits<X>::epsilon();

	// _0F_0({},{},x) = exp(x)
	{
		Hypergeometric<X> F_00({}, {});
		X x = 1;
		auto [F, eps, small, iters] = F_00.value(x, epsilon, 40, 40);
		assert(abs(F - exp(x)) <= epsilon * F);
		assert(abs(eps) <= epsilon);
		assert(iters == 40);
	}
	{
		Hypergeometric<X> F_00({}, {});
		X x = 1;
		X eps = epsilon;
		int skip = 1;
		int terms = 40;
		double t0, t1;
		auto f0 = [x, eps, skip, terms, &F_00]() { return F_00.value(x, eps, skip, terms); };
		t0 = fms::test::time(f0);
		auto f1 = [x]() { return exp(x); };
		t1 = fms::test::time(f1);
		//assert(t0 < 100 * t1); // not horrible
		t0 = fms::test::time(fms::test::repeat(100, f0));
		t1 = fms::test::time(fms::test::repeat(100, f1));
		//assert(t0 < 100 * t1); // not horrible
	}
	{
		X xs[] = { X(-1), X(0), X(1) };

		X tol = 2;
		for (X x : xs) {
			Hypergeometric<X> F_00({}, {});
			auto [F, eps, small, iters] = F_00.value(x, epsilon, 1, 40);
			X expx = exp(x);
			/*
			X TOL;
			TOL = (F - expx);
			TOL /= F * epsilon;
			TOL = (eps) / (F * epsilon);
			*/
			assert(abs(F - expx) <= tol * F * epsilon);
			assert(abs(eps) <= tol * F * epsilon);
			assert(small == 1);
			if (x == 1) {
				assert(iters == (IF_DOUBLE(19) : IF_FLOAT(11) : 0));
			}
		}
	}
	// (1 - x)^n = _1F_0(-n;x) = sum_{k=0}^n (-1)^k C(n,k) x^k
	{
		X xs[] = { X(-1.5), X(-1), X(0), X(1), X(2)};
		int ns[] = { 0, 1, 2, 3, 4 };

		X tol = 1;
		for (int n : ns) {
			for (X x : xs) {
				Hypergeometric<X> F_10({X(-n)}, {});
				auto [F, eps, small, iters] = F_10.value(x, epsilon, 1, 40);
				X powx = pow(1 - x, X(n));
				// /*
				X TOL;
				TOL = (F - powx);
				TOL /= F * epsilon;
				TOL = (eps) / (F * epsilon);
				// */
				assert(abs(F - powx) <= tol * abs1(F) * epsilon);
				assert(abs(eps) <= tol * abs1(F) * epsilon);
				//assert(small == 1);
				if (x == 1) {
					//assert(iters == (IF_DOUBLE(19) : IF_FLOAT(11) : 0));
				}

			}
		}
	}
	// arcsin(x) = x F(1/2, 1/2; 3/2; x^2),
	{
		X xs[] = { X(-.5), X(0.1), X(0.5) };

		X tol = 1;
		for (X x : xs) {
			Hypergeometric<X> F_21({ X(.5), X(.5) }, { X(1.5) });
			auto [F, eps, small, iters] = F_21.value(x*x, epsilon, 1, 40);
			F *= x;
			X asinx = asin(x);
			// /*
			X TOL;
			TOL = (F - asinx);
			TOL /= F * epsilon;
			TOL = (eps) / (F * epsilon);
			// */
			assert(abs(F - asinx) <= tol * abs1(F) * epsilon);
			assert(abs(eps) <= tol * abs1(F) * epsilon);
			//assert(small == 1);
			if (x == 1) {
				//assert(iters == (IF_DOUBLE(19) : IF_FLOAT(11) : 0));
			}

		}
	}

	return 0;
}
int test_hypergeometric_d = test_hypergeometric<double>();
int test_hypergeometric_f = test_hypergeometric<float>();