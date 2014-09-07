// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Exponential distribution  (a.k.a. negative exponential distribution). 
// Parameters:
// λ > 0: rate, or inverse scale
// Support: x ∈ [0; ∞).

import (
	"math/rand"
)

// ExponentialPDF returns the PDF of the Exponential distribution. 
func ExponentialPDF(λ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		return λ * exp(-λ*x)
	}
}

// ExponentialLnPDF returns the natural logarithm of the PDF of the Exponential distribution. 
func ExponentialLnPDF(λ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return negInf
		}
		return log(λ) - λ*x
	}
}

// ExponentialPDFAt returns the value of PDF of Exponential distribution at x. 
func ExponentialPDFAt(λ, x float64) float64 {
	pdf := ExponentialPDF(λ)
	return pdf(x)
}

// ExponentialCDF returns the CDF of the Exponential distribution. 
func ExponentialCDF(λ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		return 1 - exp(-λ*x)
	}
}

// ExponentialCDFAt returns the value of CDF of the Exponential distribution, at x. 
func ExponentialCDFAt(λ, x float64) float64 {
	cdf := ExponentialCDF(λ)
	return cdf(x)
}

// ExponentialQtl returns the inverse of the CDF (quantile) of the Exponential distribution. 
func ExponentialQtl(λ float64) func(p float64) float64 {
	// p: probability for which the quantile is evaluated
	return func(p float64) float64 {
		return -log(1-p) / λ
	}
}

// ExponentialQtlFor returns the inverse of the CDF (quantile) of the Exponential distribution, for given probability.
func ExponentialQtlFor(λ, p float64) float64 {
	cdf := ExponentialQtl(λ)
	return cdf(p)
}

// ExponentialNext returns random number drawn from the Exponential distribution. 
func ExponentialNext(λ float64) float64 { return rand.ExpFloat64() / λ }

// Exponential returns the random number generator with  Exponential distribution. 
func Exponential(λ float64) func() float64 { return func() float64 { return ExponentialNext(λ) } }

// ExponentialMean returns the mean of the Exponential distribution. 
func ExponentialMean(λ float64) float64 {
	return 1 / λ
}

// Exponential returns the median of the Exponential distribution. 
func ExponentialMedian(λ float64) (med float64) {
	return (1 / λ) * log(2)
}

// ExponentialMode returns the mode of the Exponential distribution. 
func ExponentialMode(λ float64) float64 {
	return 0
}

// ExponentialVar returns the variance of the Exponential distribution. 
func ExponentialVar(λ float64) float64 {
	return 1 / (λ * λ)
}

// ExponentialStd returns the standard deviation of the Exponential distribution. 
func ExponentialStd(λ float64) float64 {
	return 1 / λ
}

// ExponentialSkew returns the skewness of the Exponential distribution. 
func ExponentialSkew(λ float64) (s float64) {
	return 2
}

// ExponentialExKurt returns the excess kurtosis of the Exponential distribution. 
func ExponentialExKurt(λ float64) float64 {
	return 6
}

// ExponentialMGF returns the moment-generating function of the Exponential distribution. 
func ExponentialMGF(λ, p, t float64) float64 {
	return 1 / (1 - t/λ)
}
