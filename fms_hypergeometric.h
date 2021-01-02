// fms_hypergeometric.h - Hypergeometric function
#pragma once
#include <cmath>
#include <concepts>
#include <initializer_list>
#include <limits>
#include <tuple>

namespace fms {

	template<class X>
	using list = std::initializer_list<X>;

	// primitive implementation of general hypergeometric function
	template<class X>
		requires std::is_floating_point_v<X>
	class Hypergeometric {
		const list<X>& a;
		const list<X>& b;
		X n;  // current n
		X an; // (a)_n
		X bn; // (b)_n
		X xn; // x^n
		X n_; // n!
		X pFq; // running value
	public:
		Hypergeometric(const list<X>& a, const list<X>& b)
			: a(a), b(b), n(0), an(1), bn(1), xn(1), n_(1), pFq(0)
		{ }

		// return next term
		X next(X x)
		{
			for (X ai : a) {
				an *= ai + n;
			}
			for (X bi : b) {
				bn *= bi + n;
			}

			X dF = (an / bn) * xn / n_;

			xn *= x;
			n_ *= ++n;

			return dF;
		}

		// idempotent
		X regularized() const
		{
			X F = pFq;

			for (X bi : b) {
				F /= std::tgamma(bi);
			}

			return F;
		}

		// square root of machine epsilon
		static constexpr X sqrt_eps = X(1) / (1ul << (std::numeric_limits<X>::digits / 2));
		
		// policy based convergence
		std::tuple<X, X, int, int> value(X x, X eps = sqrt_eps, int skip = 40, int terms = 40)
		{
			X dF, maxF = 1;
			int ignore = skip;
			int small = 0;
			int iters = 0;

			// if (a)_n = 0 then all follwing terms are 0
			while (an and ignore and terms - iters) {

				dF = next(x);
				pFq += dF;
				maxF = std::max(maxF, abs(pFq));

				if (abs(dF) < maxF * eps) {
					++small;
					--ignore;
				}
				else {
					ignore = skip;
				}

				++iters;
			}

			return std::tuple(pFq, dF, small, iters);
		}
	};

	// pFq(a,b,x) = sum_n (a_1)_n ... (a_p)_n/((b_1)_n ... (b_q)_n) x^n/n!
	template<class X> requires std::is_floating_point_v<X>
	static X HypergeometricPFQ(const list<X>& a, const list<X>& b, X x, bool regularized = false)
	{
		Hypergeometric pFq(a, b);

		X F = pFq.value(x);

		return regularized ? pFq.regularized() : F;
	}

}