// fms_test.h - testing utility
#pragma once
#include <cassert>
#include <algorithm>
#include <chrono>
#include <functional>
#include <initializer_list>
#include <valarray>

namespace fms::test {

	template<class F>
	inline std::function<void(void)> repeat(size_t n, const F& f)
	{
		return [&n, f]() { while (n--) f(); };
	}

	// time in milliseconds
	template<class F, class C = std::chrono::high_resolution_clock>
	inline double time(const F& f)
	{
		auto start = C::now();
		f();
		auto stop = C::now();

		return 1000*std::chrono::duration<double>(stop - start).count();
	}

	// f'(x) + f'''(x)h^2/6 + ...
	template<class F, class X>
	inline X diff(const F& f, X x, X h)
	{
		return (f(x + h) - f(x - h)) / (2 * h);
	}

	// error is O(h^2) where O ~ f'''/6
	template<class F, class X>
	inline void check(X df, const F& f, X x, X h, X O = 150)
	{
		X f1 = diff(f, x, h);
		O = O * std::max({ X(1), abs(df), abs(f1) });
		X o = fabs(df - f1) / (h * h);
		X o_ = 1 / o;

		assert(fabs(df - f1) < O * (h * h));
	}

	// range of x and dx
	template<class F, class dF, class Xs>
	inline void check(const F& f, const dF& df, const Xs& xs, const Xs& hs)
	{
		for (auto x : xs) {
			auto dfx = df(x);
			for (auto h : hs) {
				check(dfx, f, x, h);
			}
		}
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

}
