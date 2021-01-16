// fms_variate_logistic
#pragma once
#include <concepts>
#include <initializer_list>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_hyperg.h>
#include "fms_ensure.h"
#include "fms_sf_hypergeometric.h"

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

		// A_{n,k} = - (b + k) A_{n-1, k} + (a + b + k - 1) A_{n-1, k-1}, A_{0,0} = 1
		template<class X = double>
		inline X A(X a, X b, unsigned n, unsigned k)
		{
			if (k > n) {
				return 0;
			}
			if (n == 0 and k == 0) {
				return 1;
			}

			return -(b + k) * A(a, b, n - 1, k) + (a + b + k - 1) * A(a, b, n - 1, k - 1);
		}
	}

#ifdef _DEBUG
	template<class X>
	inline void check_A(X a, X b)
	{
		static constexpr X eps = std::numeric_limits<X>::epsilon();
		auto eq = [](X x, X y) { return abs(x - y) <= 10*eps; };
		assert(eq(A(a, b, 0, 0), 1));
		assert(eq(A(a, b, 0, 1), 0));
		assert(eq(A(a, b, 0, -1), 0));

		assert(eq(A(a, b, 1, 0), -b));
		assert(eq(A(a, b, 1, 1), a + b));
		assert(eq(A(a, b, 1, 2), 0));
		assert(eq(A(a, b, 1, -1), 0));

		assert(eq(A(a, b, 2, 0), b * b));
		assert(eq(A(a, b, 2, 1), -(a + b) * (1 + 2 * b)));
		assert(eq(A(a, b, 2, 2), (1 + a + b) * (a + b)));
		assert(eq(A(a, b, 2, 3), 0));
		assert(eq(A(a, b, 2, -1), 0));
	}
#endif // _DEBUT

	// generalized logistic density is f(a,b;x) = e^{-b x}/(1 + e^{-x})^{a + b} / B(a,b)
	template<class X = double, class S = X>
		requires std::is_floating_point_v<X> && std::is_floating_point_v<S>
	struct logistic {
		typedef X xtype;
		typedef S stype;
		X a, b;
		logistic(X a = 1, X b = 1)
			: a(a), b(b)
		{ }

		// (d/dx)^n f(x) = sum_{k=0}^n A_{n,k} e^{-(b + k) x}/(1 + e^{-x})^{a + b + k}
		static X cdf0(X a, X b, X x, unsigned n = 0)
		{
			ensure(a > 0 and b > 0);

			X e_x = exp(-x);
	
			if (n == 0) {
				return beta_inc(a, b, 1 / (1 + e_x));
			}

			unsigned n_ = n - 1;
			X e_ = e_x / (1 + e_x);
			X e_k = 1; // e_^k
			X Ak = 0;
			for (unsigned k = 0; k <= n_; ++k) {
				Ak += A(a, b, n_, k) * e_k;
				e_k *= e_;
			}

			return exp(-b * x) * pow(1 + e_x, -a - b) * Ak / gsl_sf_beta(a, b);
		}
		X cdf(X x, S s = 0, unsigned n = 0) const
		{
			ensure(-a < s and s < b);

			if (n == 0) {
				return beta_inc(a + s, b - s, 1/(1 + exp(-x)));
			}

			return cdf0(a + s, b - s, x, n);
		}
		S cumulant(S s, unsigned n = 0) const
		{
			ensure(-1 < s and s < 1);

			if (n == 0) {
				return gsl_sf_lngamma(a + s) - gsl_sf_lngamma(a) 
				     + gsl_sf_lngamma(b - s) - gsl_sf_lngamma(b);
			}

			int n_ = static_cast<int>(n - 1);

			return gsl_sf_psi_n(n_, a + s) + ((n_&1) ? 1 : -1) * gsl_sf_psi_n(n_, b - s);
		}
		static X edf(X x, S s)
		{
			X u = cdf0(1, 1, x, 0);

			return beta_inc_1(1 + s, 1 - s, u) - beta_inc_2(1 + s, 1 - s, u);
		}

		static X beta(X a, X b)
		{
			return gsl_sf_beta(a, b);
		}
		static X beta_1(X a, X b)
		{
			return gsl_sf_beta(a, b) * (gsl_sf_psi_n(0, a) - gsl_sf_psi_n(0, a + b));
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
