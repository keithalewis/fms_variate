# NOTES

Variates have a cumulative distribution and cumulant and their derivatives.

```C++
normal N;
N.cdf(x, n);
N.cumulant(s, n);

mean(N) -> N.cumulant(0, 1);
variance(N) -> N.cumulant(0, 2);
standard_deviation(N) ???
```

Global operators constant times variate, sum of (independent)  variates.