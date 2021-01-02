// fms_variate_logistic
#pragma once
#include <concepts>
#include <initializer_list>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_hyperg.h>
#include "fms_ensure.h"
#include "fms_hypergeometric.h"

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

		// (d/dx)^n 1/(1 + e^{-x} = sum_{k=1}^n A_{n,k} e^{-k x}/(1 + e^{-x})^{k + 1}
		static X cdf0(X x, size_t n = 0)
		{
			X e = exp(-x);

			if (n == 0) {
				return 1/(1 + e);
			}

			X f = 0;
			X e_ = e;
			X e1_ = 1 + e;
			for (size_t k = 1; k <= n; ++k) {
				e1_ *= 1 + e;
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
				return beta_inc(1 + s, 1 - s, 1/(1 + exp(-x)));
			}

			X f = 0;
			S sk = 1; // s^k
			for (size_t k = 0; k < n; ++k) {
				f += gsl_sf_choose((unsigned int)n - 1, (unsigned int)k) * cdf0(x, n - k) * sk;
				sk *= s;
			}

			return exp(s * x - cumulant(s, 0)) * f;
		}
		static S cumulant(S s, size_t n = 0)
		{
			ensure(-1 < s and s < 1);

			if (n == 0) {
				return gsl_sf_lngamma(1 + s) + gsl_sf_lngamma(1 - s);
			}

			int n_ = static_cast<int>(n - 1);

			return gsl_sf_psi_n(n_, 1 + s) + ((n_&1) ? 1 : -1) * gsl_sf_psi_n(n_, 1 - s);
		}
		static X edf(X x, S s)
		{
			X u = cdf0(x, 0);

			return beta_inc_1(1 + s, 1 - s, u) - beta_inc_2(1 + s, 1 - s, u);
		}

		static X beta_inc(X a, X b, X u)
		{
			return gsl_sf_beta_inc(a, b, u);
		}

		// d/da I_u(a, b) = (log u - psi(a) + psi(a + b)) I_u(a,b)
		static X beta_inc_1(X a, X b, X u)
		{
			X Iu = beta_inc(a, b, u);

			return (log(u) - gsl_sf_psi_n(0, a) + gsl_sf_psi_n(0, a + b)) * Iu;
		}
		// I_u(a, b) = 1 - I_{1 - u}(b, a)
		// d/db I_u(a, b) = d/db 1 - I_{1 - u}(b, a)
		static X beta_inc_2(X a, X b, X u)
		{
			X Iu_ = beta_inc(b, a, 1 - u);

			return -(log(1 - u) - gsl_sf_psi_n(0, b) + gsl_sf_psi_n(0, b + a)) * Iu_;
		}
	};

}
