// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Cauchy distribution. 
// A continuous probability distribution, also known as Lorentz distribution, Cauchy–Lorentz distribution, Lorentz(ian) function, or Breit–Wigner distribution. The simplest Cauchy distribution is called the standard Cauchy distribution. It has the distribution of a random variable that is the ratio of two independent standard normal random variables. 
//
// Parameters: 
// δ ∈ R		location
// γ > 0		scale
//
// Support: 
// x ∈ R

import (
	"math/rand"
)

// CauchyPDF returns the PDF of the Cauchy distribution. 
func CauchyPDF(δ, γ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(δ) || isNaN(γ) {
			return x + δ + γ
		}
		if γ <= 0 {
			return NaN
		}
		y := (x - δ) / γ
		return 1 / (π * γ * (1 + y*y))
	}
}

// CauchyLnPDF returns the natural logarithm of the PDF of the Cauchy distribution. 
func CauchyLnPDF(δ, γ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(δ) || isNaN(γ) {
			return x + δ + γ
		}
		if γ <= 0 {
			return NaN
		}
		y := (x - δ) / γ
		return -log(π * γ * (1 + y*y))
	}
}

// CauchyPDFAt returns the value of PDF of Cauchy distribution at x. 
func CauchyPDFAt(δ, γ, x float64) float64 {
	pdf := CauchyPDF(δ, γ)
	return pdf(x)
}

// CauchyCDF returns the CDF of the Cauchy distribution. 
func CauchyCDF(δ, γ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(δ) || isNaN(γ) {
			return x + δ + γ
		}
		if γ <= 0 {
			return NaN
		}

		x = (x - δ) / γ
		if isNaN(x) {
			return NaN
		}
		if isInf(x, -1) {
			return 0
		}
		if isInf(x, 1) {
			return 1
		}

		// for large x, the standard formula suffers from cancellation.
		//This is from Morten Welinder thanks to  Ian Smith's  atan(1/x)
		if abs(x) > 1 {
			y := atan(1/x) / π
			if x > 0 {
				return 0.5 - y + 0.5
			} else {
				return -y
			}
		}
		return 0.5 + atan(x)/π
	}
}

// CauchyLnCDF returns the CDF of the Cauchy distribution. 
func CauchyLnCDF(δ, γ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(δ) || isNaN(γ) {
			return x + δ + γ
		}
		if γ <= 0 {
			return NaN
		}

		x = (x - δ) / γ
		if isNaN(x) {
			return NaN
		}
		if isInf(x, -1) {
			return negInf
		}
		if isInf(x, 1) {
			return 0
		}

		// for large x, the standard formula suffers from cancellation.
		//This is from Morten Welinder thanks to  Ian Smith's  atan(1/x)
		if abs(x) > 1 {
			y := atan(1/x) / π
			if x > 0 {
				return log1p(-y)
			} else {
				return log(-y)
			}
		}
		return log(0.5 + atan(x)/π)
	}
}

// CauchyCDFAt returns the value of CDF of the Cauchy distribution, at x. 
func CauchyCDFAt(δ, γ, x float64) float64 {
	cdf := CauchyCDF(δ, γ)
	return cdf(x)
}

// CauchyQtl returns the inverse of the CDF (quantile) of the Cauchy distribution. 
func CauchyQtl(δ, γ float64) func(p float64) float64 {
	return func(p float64) float64 {
		if isNaN(p) || isNaN(δ) || isNaN(γ) {
			return p + δ + γ
		}
		if γ < 0 || isInf(γ, 0) {
			return NaN
		}
		if γ == 0 {
			return δ
		}
		if p == 1 {
			return posInf
		}
		//	-1/tan(π * p) = -cot(π * p) = tan(π * (p - 1/2)) 
		return δ - γ/tan(π*p)
	}
}

// CauchyQtlFor returns the inverse of the CDF (quantile) of the Cauchy distribution, for given probability.
func CauchyQtlFor(δ, γ, p float64) float64 {
	qtl := CauchyQtl(δ, γ)
	return qtl(p)
}

// CauchyNext returns random number drawn from the Cauchy distribution. 
func CauchyNext(δ, γ float64) float64 {
	//	p := UniformNext(0, 1)
	//	return CauchyQtlFor(δ, γ, p)
	return γ*tan(π*(rand.Float64()-0.5)) + δ // Nolan 2009: 21, Eq. 1.11
}

// Cauchy returns the random number generator with  Cauchy distribution. 
func Cauchy(δ, γ float64) func() float64 {
	return func() float64 { return CauchyNext(δ, γ) }
}

// CauchyMean is not defined. 

// CauchyMode returns the mode of the Cauchy distribution. 
func CauchyMode(δ, γ float64) float64 {
	return δ
}

// CauchyMedian returns the median of the Cauchy distribution. 
func CauchyMedian(δ, γ float64) float64 {
	return δ
}

// CauchyVar is not defined. 

// CauchyStd is not defined. 

// CauchySkew is not defined. 

// CauchyExKurt is not defined. 

// CauchyMGF does not exist.
