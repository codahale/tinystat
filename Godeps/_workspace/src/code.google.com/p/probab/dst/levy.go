// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Lévy distribution is a continuous probability distribution for a non-negative random variable. In spectroscopy this distribution, with frequency as the dependent variable, is known as a van der Waals profile. It is a special case of the inverse-gamma distribution.
//
// Parameters: 
// δ > 0		location
// γ > 0		scale
//
// Support: 
// x ∈ [δ, ∞)

// LevyPDF returns the PDF of the Lévy distribution. 
func LevyPDF(δ, γ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(δ) || isNaN(γ) {
			return x + δ + γ
		}
		if γ <= 0 || δ <= 0 {
			return NaN
		}
		return exp(log(γ/(2*π))/2 - 3*log(x-δ)/2 - γ/(2*(x-δ)))
	}
}

// LevyLnPDF returns the natural logarithm of the PDF of the Lévy distribution. 
func LevyLnPDF(δ, γ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(δ) || isNaN(γ) {
			return x + δ + γ
		}
		if γ <= 0 || δ <= 0 {
			return NaN
		}
		return log(γ/(2*π))/2 - 3*log(x-δ)/2 - γ/(2*(x-δ))
	}
}

// LevyPDFAt returns the value of PDF of Lévy distribution at x. 
func LevyPDFAt(δ, γ, x float64) float64 {
	pdf := LevyPDF(δ, γ)
	return pdf(x)
}

// LevyCDF returns the CDF of the Lévy distribution. 
func LevyCDF(δ, γ float64) func(x float64) float64 {
	return func(x float64) float64 {

		if isNaN(x) || isNaN(δ) || isNaN(γ) {
			return x + δ + γ
		}
		if γ <= 0 || δ <= 0 {
			return NaN
		}
		return 2 * (1 - ZCDFAt(1/sqrt((x-δ)/γ)))
	}
}

// LevyCDFAt returns the value of CDF of the Lévy distribution, at x. 
func LevyCDFAt(δ, γ, x float64) float64 {
	cdf := LevyCDF(δ, γ)
	return cdf(x)
}

// LevyQtl returns the inverse of the CDF (quantile) of the Lévy distribution. 
func LevyQtl(δ, γ float64) func(p float64) float64 {
	return func(p float64) float64 {
		if isNaN(p) || isNaN(δ) || isNaN(γ) {
			return p + δ + γ
		}
		if γ <= 0 || δ <= 0 {
			return NaN
		}
		if p == 1 {
			return posInf
		}
		zq := ZQtlFor(1 - p/2)
		return γ/(zq*zq) + δ
	}
}

// LevyQtlFor returns the inverse of the CDF (quantile) of the Lévy distribution, for given probability.
func LevyQtlFor(δ, γ, p float64) float64 {
	qtl := LevyQtl(δ, γ)
	return qtl(p)
}

// LevyNext returns random number drawn from the Lévy distribution. 
func LevyNext(δ, γ float64) float64 {
	z := NormalNext(0, 1)
	return γ/(z*z) + δ // Nolan 2009: 21, Eq. 1.12
}

// Levy returns the random number generator with  Lévy distribution. 
func Levy(δ, γ float64) func() float64 {
	return func() float64 { return LevyNext(δ, γ) }
}

// LevyMean returns the mean of the Lévy distribution. 
func LevyMean(δ, γ float64) float64 {
	return posInf
}

// LevyMode returns the mode of the Lévy distribution. 
func LevyMode(δ, γ float64) float64 {
	if δ == 0 {
		return γ / 3
	}
	return NaN
}

// LevyMedian returns the median of the Lévy distribution. 
func LevyMedian(δ, γ float64) float64 {
	e := 1 / (erfc(0.5))
	return γ * e * e / 2
}

// LevyVar returns the variance of the Lévy distribution. 
func LevyVar(δ, γ float64) float64 {
	return posInf
}

// LevyStd returns the standard deviation of the Lévy distribution. 
func LevyStd(δ, γ float64) float64 {
	return posInf
}

// LevySkew is not defined. 

// LevyExKurt is not defined. 

// LevyMGF does not exist.
