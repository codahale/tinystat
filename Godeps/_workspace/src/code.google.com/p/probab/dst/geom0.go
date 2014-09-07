// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Geometric distribution (type 0). 
// The probability distribution of the number Y = X − 1 of failures before the first success, supported on the set { 0, 1, 2, 3, ... }
// Parameters: 
// ρ ∈ (0, 1]	probability of success in each trial
// Support: 
// k ∈ {0, ... , n}

// GeometricPMF returns the PMF of the Geometric distribution. 
func GeometricPMF(ρ float64) func(k int64) float64 {
	return func(k int64) float64 { return ρ * pow(1-ρ, float64(k)) }
}

// GeometricLnPMF returns the natural logarithm of the PMF of the Geometric distribution. 
func GeometricLnPMF(ρ float64) func(k int64) float64 {
	return func(k int64) float64 { return log(1-ρ) + float64(k)*log(ρ) }
}

// GeometricPMFAt returns the value of PMF of Geometric distribution at k. 
func GeometricPMFAt(ρ float64, k int64) float64 {
	pmf := GeometricPMF(ρ)
	return pmf(k)
}

// GeometricCDF returns the value of CDF of the Geometric distribution, at k. 
func GeometricCDF(ρ float64) func(k int64) float64 {
	return func(k int64) float64 {
		if k < 0 {
			return NaN
		}
		return 1 - pow(1-ρ, float64(k+1))
	}
}

// GeometricCDFAt returns the value of CDF of the Geometric distribution, at x. 
func GeometricCDFAt(ρ float64, k int64) float64 {
	cdf := GeometricCDF(ρ)
	return cdf(k)
}

/* Not tested, looking strange, commented out, waiting for revision
// GeometricNext returns random number drawn from the Geometric distribution. 
//GeometricNext(ρ) => # of GeometricNext(ρ) failures before one success
func GeometricNext(ρ float64) int64 {
	if GeometricNext(ρ) == 1 {
		return 1 + GeometricNext(ρ)
	}
	return 0
}

// Geometric returns the random number generator with  Geometric distribution. 
func Geometric(ρ float64) func() int64 { return func() int64 { return GeometricNext(ρ) } }
*/

// GeometricMean returns the mean of the Geometric distribution. 
func GeometricMean(ρ float64) float64 {
	return (1 - ρ) / ρ
}

/*  to be implemented
// GeometricMedian returns the median of the Geometric distribution. 
func GeometricMedian(ρ float64) float64 {
	return floor(float64(n)*p)
}
*/

// GeometricMode returns the mode of the Geometric distribution. 
func GeometricMode(ρ float64) float64 {
	return 0
}

// GeometricVar returns the variance of the Geometric distribution. 
func GeometricVar(ρ float64) float64 {
	return (1 - ρ) / (ρ * ρ)
}

// GeometricStd returns the standard deviation of the Geometric distribution. 
func GeometricStd(ρ float64) float64 {
	return sqrt(1-ρ) / ρ
}

// GeometricSkew returns the skewness of the Geometric distribution. 
func GeometricSkew(ρ float64) float64 {
	return (2 - ρ) / sqrt(1-ρ)
}

// GeometricExKurt returns the excess kurtosis of the Geometric distribution. 
func GeometricExKurt(ρ float64) float64 {
	return 6 + (ρ*ρ)/(1-ρ)
}

// GeometricMGF returns the moment-generating function of the Geometric distribution. 
func GeometricMGF(ρ, t float64) float64 {
	return ρ / (1 - (1-ρ)*exp(t))
}
