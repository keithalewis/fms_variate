// fms_variate_interface.h - Random standare variates.
#pragma once
#include <cmath>

namespace fms::variate {

	inline const char interfac_doc[] = R"(
A random variable \(X\) is determined by its cumulative distribution function \(F(x) = P(X <= x)\). 
Its cumulant is \(kappa(s) = \log E[\exp(s X)]\) and its Esscher transform \(X_s\) is defined by 
\(P(X_s \le x) = P_s(X \le x) = E[1(X \le x) \exp(s X - kappa(s))]\).
A variate must implement the derivatives of the cdf of \(X_s\), cumulant, and
the derivative of the cdf of \(X_s\) with respect to \(s\), called sega.
	)";
	template<class X = double, class S = double>
	struct interface {
		typedef X xtype;
		typedef S stype;
		virtual ~interface() {}
		X cdf(X x, S s, unsigned n = 0) const
		{
			return cdf_(x, s, n);
		}
		X pdf(X x, S s) const
		{
			return pdf_(x, s);
		}
		X mgf(X s) const
		{
			return mgf_(s);
		}
		X cgf(S s) const
		{
			return cgf_(s);
		}
		X sdf(X x, S s) const
		{
			return sdf_(x, s);
		}
		S cumulant(S s, unsigned n) const
		{
			return cumulant_(s, n);
		}
	private:
		virtual X cdf_(X x, S s, unsigned n) const = 0;
		virtual X pdf_(X x, S s) const = 0;
		virtual X mgf_(S s) const = 0;
		virtual X cgf_(S s) const = 0;
		virtual X sdf_(X x, S s) const = 0;
		virtual S cumulant_(S s, unsigned n) const = 0;
	};


} // namespace fms::variate
