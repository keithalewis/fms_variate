// fms_hypergeometric.h - Hypergeometric function
#pragma once
#include <concepts>
#include <initializer_list>
#include <gsl/gsl_sf_gamma.h>

namespace fms {

	template<class X>
	using list = std::initializer_list<X>;

	// primitive implementation
	// pFq(a,b,x) = sum_n (a)_n/(b)_n x^n/n!
	template<class X> requires std::is_floating_point_v<X>
	static X HypergeometricPFQ(const list<X>& a, const list<X>& b, X x, bool regularized = false)
	{
		X pFq = 0;

		X an = 1; // (a)_n
		X bn = 1; // (b)_n
		X xn = 1; // x^n
		X n_ = 1; // n!

		for (X n = 0; n < 20; ) {
			for (X ai : a) {
				an *= ai + n;
			}
			for (X bi : b) {
				bn *= bi + n;
			}
			pFq += (an / bn) * xn / n_;
			++n;
			xn *= x;
			n_ *= n;
		}

		if (regularized) {
			for (X bi : b) {
				pFq /= gsl_sf_gamma(bi);
			}
		}

		return pFq;
	}

}