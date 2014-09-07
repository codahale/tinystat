// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Poisson distribution. 
// A discrete probability distribution that expresses the probability of a given number of events occurring in a fixed interval of time and/or space if these events occur with a known average rate and independently of the time since the last event. (The Poisson distribution can also be used for the number of events in other specified intervals such as distance, area or volume.)
// Frank A. Haight (1967). Handbook of the Poisson Distribution. New York: John Wiley & Sons.
//
// Parameters: 
// λ > 0 (real) 	average rate
// p ∈ [0, 1]	probability of success in each trial
//
// Support: 
// k ∈ {0, ... , n}
// x ∈ (0, ∞)

/*
func PoissonLnPMF(λ float64) (foo func(i int64) float64) {
	pmf := PoissonPMF(λ)
	return func(i int64) (p float64) {
		return log(pmf(i))
		//p = -λ +log(λ)*float64(i)
		//x := log(Γ(float64(i)+1))
		// = x
		//p -= LnΓ(float64(i)+1)
		//return p
	}
}

func PoissonPMF(λ float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := ExponentialNext(-λ) * pow(λ, float64(k)) / Γ(float64(k)+1)
		return p
	}
}

func PoissonPMF(λ float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := exp(-λ) * pow(λ, float64(k)) / Γ(float64(k)+1)
		return p
	}
}
*/

// PoissonPMF returns the PMF of the Poisson distribution. 
func PoissonPMF(λ float64) func(k int64) float64 {
	pmf := PoissonLnPMF(λ)
	return func(k int64) float64 {
		p := exp(pmf(k))
		return p
	}
}

// PoissonLnPMF returns the natural logarithm of the PMF of the Poisson distribution. 
func PoissonLnPMF(λ float64) func(k int64) float64 {
	return func(k int64) (p float64) {
		i := float64(k)
		a := log(λ) * i
		b := log(Γ(i + 1))
		p = a - b - λ
		return p
	}
}

// PoissonPMFAt returns the value of PMF of Poisson distribution at k. 
func PoissonPMFAt(λ float64, k int64) float64 {
	pmf := PoissonPMF(λ)
	return pmf(k)
}

// PoissonCDF returns the CDF of the Poisson distribution. 
func PoissonCDF(λ float64) func(k int64) float64 {
	return func(k int64) float64 {
		var p float64 = 0
		var i int64
		pmf := PoissonPMF(λ)
		for i = 0; i <= k; i++ {
			p += pmf(i)
		}
		return p
	}
}

// PoissonCDFAn returns the CDF of the Poisson distribution. Analytic solution, less precision.
func PoissonCDFAn(λ float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := exp(log(iΓint(k+1, λ)) - (logFact(float64(k))))
		return p
	}
}

// PoissonCDFAt returns the value of CDF of the Poisson distribution, at x. 
func PoissonCDFAt(λ float64, k int64) float64 {
	cdf := PoissonCDF(λ)
	return cdf(k)
}

// LnPoissonCDFAn returns the natural logarithm of the CDF of the Poisson distribution. Analytic solution, less precision.
func LnPoissonCDFAn(λ float64) func(k int64) float64 {
	return func(k int64) float64 {
		k1 := float64(k + 1)
		return log(iΓ(k1, λ)) - logFact(float64(k))
	}
}

// PoissonNext2 returns random number drawn from the Poisson distribution (old version). 
func PoissonNext2(λ float64) int64 {
	var k int64
	if λ < 100 { // Knuth algorithm for small λ
		// Donald E. Knuth (1969). Seminumerical Algorithms. The Art of Computer Programming, Volume 2. Addison Wesley.
		// this can be improved upon
		k = iZero
		t := exp(-λ)
		p := fOne
		for ; p > t; p *= UniformNext(0, 1) {
			k++
		}
		k -= 1

	} else { // use Normal approximation
		k = int64(iround(NormalNext(λ, sqrt(λ))))
	}
	return k
}

// Poisson returns the random number generator with  Poisson distribution. 
func Poisson(λ float64) func() int64 {
	return func() int64 {
		return PoissonNext(λ)
	}
}

// PoissonMean returns the mean of the Poisson distribution. 
func PoissonMean(λ float64, k int64) float64 {
	return λ
}

// PoissonMode returns the mode of the Poisson distribution. 
func PoissonMode(λ float64, k int64) float64 {
	return ceil(λ) - 1
}

// PoissonMedian returns the median of the Poisson distribution. Approximation. 
func PoissonMedian(λ float64, k int64) float64 {
	return floor(λ + 1/3 - 0.02*λ)
}

// PoissonVar returns the variance of the Poisson distribution. 
func PoissonVar(λ float64, k int64) float64 {
	return λ
}

// PoissonSkew returns the skewness of the Poisson distribution. 
func PoissonSkew(λ float64, k int64) float64 {
	return pow(λ, -0.5)
}

// PoissonExKurt returns the excess kurtosis of the Poisson distribution. 
func PoissonExKurt(λ float64, k int64) float64 {
	return 1 / λ
}
