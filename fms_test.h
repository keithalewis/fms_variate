// fms_test.h - testing utility
#pragma once
#include <cassert>
#include <algorithm>
#include <initializer_list>
#include <valarray>

namespace fms::test {

	// f'(x) + f'''(x)h^2/6 + ...
	template<class F, class X>
	inline X diff(const F& f, X x, X h)
	{
		return (f(x + h) - f(x - h)) / (2 * h);
	}

	// error is O(h^2) where O ~ f'''/6
	template<class X>
	inline void check(X df, X f1, X h, X O = 150)
	{
		O = O * std::max({ X(1), abs(df), abs(f1) });
		X r = fabs(df - f1) / (h * h);
		X r_ = 1 / r;

		assert(fabs(df - f1) < O * h * h);
	}

	// a, a + h, ... , b
	template<class X>
	inline std::valarray<X> range(X a, X b, X h)
	{
		std::valarray<X> r((size_t)(1 + (b - a)/h));

		for (size_t i = 0; i < r.size(); ++i) {
			r[i] = a + X(i) * h;
		}

		return r;
	}

	template<class M, class Xs, class Ss, class Ns, class X>
	inline void check_range(const M& m, const Xs& xs, const Ss& ss, const Ns& ns, X h) 
	{
		for (auto s : ss) {
			for (auto n : ns) {
				auto f = [m, s, n](X x) { return m.cdf(x, s, n); };
				for (auto x : xs) {
					X df = diff(f, x, h);
					check(df, m.cdf(x, s, n + 1), h);
				}
			}
		}
	}

}
