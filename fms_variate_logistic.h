// fms_variate_logistic
#pragma once
#include <concepts>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include "fms_ensure.h"

namespace fms::variate {

	namespace {
		inline constexpr int64_t A(uint64_t n, uint64_t k)
		{
			if (n == 0 or k > n) {
				return 0;
			}

			if (k == 1) {
				return (n & 1) ? 1 : -1;
			}

			if (n == 1) {
				return k == 1;
			}

			return k * (A(n - 1, k - 1) - A(n - 1, k));
		}
		// n choose k
		inline constexpr uint64_t C(uint64_t n, uint64_t k)
		{
			if (n == 0) {
				return k == 0;
			}

			if (2 * k > n) {
				k = n - k;
			}

			return (n == 0 or k > n) ? 0 : C(n - 1, k) + C(n - 1, k - 1);
		}
	}

	template<class X = double, class S = X>
		requires std::is_floating_point_v<X> && std::is_floating_point_v<S>
	struct logistic {
		typedef X xtype;
		typedef S stype;

		static X cdf0(X x, size_t n = 0)
		{
			X e = exp(-x);
			X e1 = 1 / (1 + e);

			if (n == 0) {
				return e1;
			}

			X f = 0;
			X e_ = e;
			X e1_ = e1;
			for (size_t k = 1; k <= n; ++k) {
				e1_ *= e1;
				f += A(n, k) * e_ / e1_;
				e_ *= e;
			}

			return f;
		}
		static X cdf(X x, S s = 0, size_t n = 0)
		{
			ensure(-1 < s and s < 1);

			if (s == 0) {
				return cdf0(x, n);
			}

			if (n == 0) {
				return gsl_sf_beta_inc(1 + s, 1 - s, 1/(1 + exp(-x)));
			}

			X es = exp(s * x - cumulant(s, 0));

			X f = 0;
			S sk = 1; // s^k
			for (size_t k = 0; k < n; ++k) {
				f += C(n, k) * cdf0(x, k + 1) * sk;
				sk *= s;
			}

			return es * f;
		}
		static S cumulant(S s, size_t n = 0)
		{
			ensure(-1 < s and s < 1);

			if (n == 0) {
				return gsl_sf_lngamma(1 + s) + gsl_sf_lngamma(1 - s);
			}

			int n_ = static_cast<int>(n - 1);

			return gsl_sf_psi_n(n_, 1 + s) + (n_&1 ? 1 : -1)*gsl_sf_psi_n(n_, 1 - s);
		}
		static X edf(X x, S s)
		{
			return x + s;
		}
	};

}
