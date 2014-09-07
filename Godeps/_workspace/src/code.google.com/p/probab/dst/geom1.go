// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Geometric distribution (type 1). 
// The probability distribution of the number Y = X − 1 of failures before the first success, supported on the set {1, 2, 3, ... }
// Parameters: 
// ρ ∈ (0, 1]	probability of success in each trial
// Support: 
// k ∈ {1, ... , n}

// Geometric1PMF returns the PMF of the Geometric1 distribution. 
func Geometric1PMF(ρ float64) func(k int64) float64 {
	return func(k int64) float64 { return ρ * pow(1-ρ, float64(k-1)) }
}

// Geometric1LnPMF returns the natural logarithm of the PMF of the Geometric distribution (type 1). 
func Geometric1LnPMF(ρ float64) func(k int64) float64 {
	return func(k int64) float64 { return log(1-ρ) + float64(k-1)*log(ρ) }
}

// Geometric1PMFAt returns the value of PMF of Geometric distribution (type 1) at k. 
func Geometric1PMFAt(ρ float64, k int64) float64 {
	pmf := Geometric1PMF(ρ)
	return pmf(k)
}

// Geometric1CDF returns the value of CDF of the Geometric distribution (type 1). 
func Geometric1CDF(ρ float64) func(k int64) float64 {
	return func(k int64) float64 {
		if k < 0 {
			return NaN
		}
		return 1 - pow(1-ρ, float64(k))
	}
}

// Geometric1CDFAt returns the value of CDF of the Geometric distribution (type 1) at k. 
func Geometric1CDFAt(ρ float64, k int64) float64 {
	cdf := Geometric1CDF(ρ)
	return cdf(k)
}

/* Not tested, looking strange, commented out, waiting for revision
// Geometric1Next returns random number drawn from the Geometric distribution (type 1). 
//Geometric1Next(ρ) => # of Geometric1Next(ρ) failures before one success
func Geometric1Next(ρ float64) int64 {
	if Geometric1Next(ρ) == 1 {
		return 1 + Geometric1Next(ρ)
	}
	return 0
}

// Geometric1 returns the random number generator with  Geometric distribution (type 1). 
func Geometric1(ρ float64) func() int64 { return func() int64 { return Geometric1Next(ρ) } }
*/

// Geometric1Mean returns the mean of the Geometric distribution (type 1). 
func Geometric1Mean(ρ float64) float64 {
	return 1 / ρ
}

/*  to be implemented
// Geometric1Median returns the median of the Geometric distribution (type 1). 
func Geometric1Median(ρ float64) float64 {
	return floor(float64(n)*p)
}
*/

// Geometric1Mode returns the mode of the Geometric distribution (type 1). 
func Geometric1Mode(ρ float64) float64 {
	return 1
}

// Geometric1Var returns the variance of the Geometric distribution (type 1). 
func Geometric1Var(ρ float64) float64 {
	return (1 - ρ) / (ρ * ρ)
}

// Geometric1Std returns the standard deviation of the Geometric distribution (type 1). 
func Geometric1Std(ρ float64) float64 {
	return sqrt(1-ρ) / ρ
}

// Geometric1Skew returns the skewness of the Geometric distribution (type 1). 
func Geometric1Skew(ρ float64) float64 {
	return (2 - ρ) / sqrt(1-ρ)
}

// Geometric1ExKurt returns the excess kurtosis of the Geometric distribution (type 1). 
func Geometric1ExKurt(ρ float64) float64 {
	return 6 + (ρ*ρ)/(1-ρ)
}

// Geometric1MGF returns the moment-generating function of the Geometric distribution (type 1). 
func Geometric1MGF(ρ, t float64) float64 {
	if t >= -log(1-ρ) {
		return NaN
	}
	return ρ * exp(t) / (1 - (1-ρ)*exp(t))
}
