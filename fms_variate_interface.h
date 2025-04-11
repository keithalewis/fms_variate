// fms_variate.h - Random variates.
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
		virtual ~interface() 
		{ }
		// Cumulative share measure E[exp(s X - kappa(s)) 1(X <= x)]
		X cdf(X x, S s) const
		{
			return cdf_(x, s);
		}
		// Share density exp(s X - kappa(s)) dE[X <= x]/dx
		X pdf(X x, S s) const
		{
			return pdf_(x, s);
		}
		// Share vega d/ds cdf(x, s)
		X sdf(X x, S s) const
		{
			return sdf_(x, s);
		}
		// Moment generating function E[exp(s X)]
		S mgf(S s) const
		{
			return mgf_(s);
		}
		// Cumulant generating function log E[exp(s X)]
		S cgf(S s) const
		{
			return cgf_(s);
		}
	private:
		virtual X cdf_(X x, S s) const = 0;
		virtual X pdf_(X x, S s) const = 0;
		virtual X sdf_(X x, S s) const = 0;
		virtual S mgf_(S s) const = 0;
		virtual S cgf_(S s) const = 0;
	};

} // namespace fms::variate
