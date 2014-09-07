// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Bernoulli distribution.

// BernoulliPMF returns the PMF of the Bernoulli distribution. 
func BernoulliPMF(ρ float64) func(k int64) float64 {
	return func(k int64) float64 {
		if k < 0 || k > 1 {
			return NaN
		}
		if k == 1 {
			return ρ
		}
		return 1 - ρ
	}
}

// BernoulliLnPMF returns the natural logarithm of the PMF of the Bernoulli distribution. 
func BernoulliLnPMF(ρ float64) func(k int64) float64 {
	return func(k int64) float64 {
		if k == 1 {
			return log(ρ)
		}
		return log(1 - ρ)
	}
}

// BernoulliPMFAt returns the value of PMF of Bernoulli distribution at x. 
func BernoulliPMFAt(ρ float64, k int64) float64 {
	pmf := BernoulliPMF(ρ)
	return pmf(k)
}

// BernoulliCDF returns the value of CDF of the Bernoulli distribution, at x. 
func BernoulliCDF(ρ float64) func(k int64) float64 {
	return func(k int64) float64 {
		if k < 0 || k > 1 {
			return NaN
		}
		if k == 1 {
			return 1
		}
		return 1 - ρ
	}
}

// BernoulliCDFAt returns the value of CDF of the Bernoulli distribution, at x. 
func BernoulliCDFAt(ρ float64, k int64) float64 {
	cdf := BernoulliCDF(ρ)
	return cdf(k)
}

// BernoulliNext returns random number drawn from the Bernoulli distribution. 
func BernoulliNext(ρ float64) int64 {
	if UniformNext(0, 1) < ρ {
		return 1
	}
	return 0
}

// Bernoulli returns the random number generator with  Bernoulli distribution. 
func Bernoulli(ρ float64) func() int64 { return func() int64 { return BernoulliNext(ρ) } }
