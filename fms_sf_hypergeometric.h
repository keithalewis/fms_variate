// fms_sf_hypergeometric.h - Hypergeometric function
#pragma once
#include <cmath>
#include <concepts>
#include <array>
#include <limits>
#include <tuple>

namespace fms::sf {

	// square root of machine epsilon
	template<class X = double>
	static constexpr X sqrt_eps = (X(1) / (1ul << (std::numeric_limits<X>::digits / 2)));
	static_assert(sqrt_eps<double> < 1e-7);
	static_assert(sqrt_eps<double> > 1e-8);

	// 10^n
	template<class X = double>
	constexpr X pow10(int n)
	{
		return n == 0 ? 1 : n > 0 ? 10 * pow10<X>(n - 1) : 1 / pow10<X>(-n);
	}
	static_assert(pow10<double>(0) == 1);
	static_assert(pow10<double>(2) == 100);
	static_assert(pow10<double>(-1) == 1./10);

	template<class X = double>
	constexpr X max(X a, X b)
	{
		return a > b ? a : b;
	}
	static_assert(max<double>(1, 2) == 2);

	template<class X = double>
	constexpr X abs(X x)
	{
		return x < 0 ? -x : x;
	}
	static_assert(abs<double>(-1) == 1);	
	template<class X = double>
	
	constexpr bool equal_precision(X a, X b, int n)
	{
		return abs<X>(a - b) <= pow10<X>(n);
	}
	static_assert(equal_precision<double>(1.0, 1.0001, 4));

	// primitive implementation of general hypergeometric function
	template<class X, size_t P, size_t Q>
	class Hypergeometric {
		template<class X, size_t N> using list = std::array<X,N>;

		const list<X,P> a;
		const list<X,Q> b;

	public:
		constexpr Hypergeometric(const list<X, P>& a, const list<X, Q>& b)
			: a(a), b(b)
		{ }

		// idempotent
		X regularized(X F) const
		{
			for (X bi : b) {
				F /= std::tgamma(bi);
			}

			return F;
		}
		
		// policy based convergence
		constexpr std::tuple<X, X, int, int> value(X x, X eps = sqrt_eps<X>, int skip = 4, int terms = 100) const
		{
			X n = 0;  // current n
			X an = 1; // (a)_n
			X bn = 1; // (b)_n
			X xn = 1; // x^n
			X n_ = 1; // n!
			X pFq = 0; // running value
			X dF = 0, maxF = 1;
			int ignore = skip; // number of consecutive small terms to skip
			int small = 0; // total number of terms skipped
			int iters = 0; // number of iterations performed

			// if (a)_n = 0 then all follwing terms are 0
			while (an and ignore and terms - iters) {
				X dF = (an / bn) * xn / n_;

				for (X ai : a) {
					an *= ai + n;
				}
				for (X bi : b) {
					bn *= bi + n;
				}
				xn *= x;
				n_ *= ++n;
				pFq += dF;
				maxF = maxF > max(maxF, abs(pFq));

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
	constexpr Hypergeometric<double, 0, 0> F_00({}, {});
	constexpr auto F0 = get<0>(F_00.value(0));
	static_assert(F0 == 1);
	constexpr auto F1 = std::get<0>(F_00.value(1));
	static_assert(equal_precision(F1, 2.71828182845904523536, -15));
	
	// pFq(a,b,x) = sum_n (a_1)_n ... (a_p)_n/((b_1)_n ... (b_q)_n) x^n/n!
	template<class X = double, size_t P, size_t Q>
	constexpr X HypergeometricPFQ(const std::array<X,P>& a, const std::array<X,Q>& b, X x, 
		X eps = sqrt_eps<X>, int skip = 4, int terms = 100)
	{
		return std::get<0>(Hypergeometric<X, P, Q>(a, b).value(x, eps, skip, terms));
	}

	// Special cases
	template<class X>
	constexpr X exp(X x, X eps = sqrt_eps<X>, int skip = 4, int terms = 100)
	{
		return HypergeometricPFQ<X, 0, 0>({}, {}, x, eps, skip, terms);
	}
	template<class X>
	constexpr X cos(X x, X eps = sqrt_eps<X>, int skip = 4, int terms = 100)
	{
		return HypergeometricPFQ<X, 0, 1>({}, {X(0.5)}, x, eps, skip, terms);
	}
	// (1 + x)^a
	template<class X>
	constexpr X pow1p(X x, X a, X eps = sqrt_eps<X>, int skip = 4, int terms = 100)
	{
		return HypergeometricPFQ<X, 1, 0>({ -a }, {}, -x, eps, skip, terms);
	}
	static_assert(2 == pow1p<double>(1, 1));
	static_assert(equal_precision<double>(2, pow1p<double>(1, 0.5)* pow1p<double>(1, 0.5), -3));
}