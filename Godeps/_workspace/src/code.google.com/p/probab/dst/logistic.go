// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Logistic distribution. 
// A continuous probability distribution. Its cumulative distribution function is the logistic function, which appears in logistic regression and feedforward neural networks. It resembles the normal distribution in shape but has heavier tails (higher kurtosis).
//
// Parameters: 
// μ ∈ R		location
// σ > 0		scale
//
// Support: 
// x ∈ R

func log1pexp(x float64) float64 {
	if x <= 18 {
		return log1p(exp(x))
	}
	if x > 33.3 {
		return x
	}
	return x + exp(-x)
}

// LogisticPDF returns the PDF of the Logistic distribution. 
func LogisticPDF(μ, σ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(μ) || isNaN(σ) {
			return x + μ + σ
		}
		if σ <= 0 {
			return NaN
		}

		x = abs((x - μ) / σ)
		e := exp(-x)
		f := 1.0 + e
		return e / (σ * f * f)
	}
}

// LogisticLnPDF returns the natural logarithm of the PDF of the Logistic distribution. 
func LogisticLnPDF(μ, σ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(μ) || isNaN(σ) {
			return x + μ + σ
		}
		if σ <= 0 {
			return NaN
		}

		x = abs((x - μ) / σ)
		e := exp(-x)
		f := 1.0 + e
		return -(x + log(σ*f*f))
	}
}

// LogisticPDFAt returns the value of PDF of Logistic distribution at x. 
func LogisticPDFAt(μ, σ, x float64) float64 {
	pdf := LogisticPDF(μ, σ)
	return pdf(x)
}

// LogisticCDF returns the CDF of the Logistic distribution. 
func LogisticCDF(μ, σ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(μ) || isNaN(σ) {
			return x + μ + σ
		}
		if σ <= 0 {
			return NaN
		}
		x = (x - μ) / σ
		if isNaN(x) {
			return NaN
		}
		return 1 / (1 + exp(-x))
	}
}

// LogisticLnCDF returns the CDF of the Logistic distribution. 
func LogisticLnCDF(μ, σ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(μ) || isNaN(σ) {
			return x + μ + σ
		}
		if σ <= 0 {
			return NaN
		}
		x = (x - μ) / σ
		if isNaN(x) {
			return NaN
		}
		return -log1pexp(-x)
	}
}

// LogisticCDFAt returns the value of CDF of the Logistic distribution, at x. 
func LogisticCDFAt(μ, σ, x float64) float64 {
	cdf := LogisticCDF(μ, σ)
	return cdf(x)
}

// LogisticQtl returns the inverse of the CDF (quantile) of the Logistic distribution. 
func LogisticQtl(μ, σ float64) func(p float64) float64 {
	return func(p float64) float64 {
		if isNaN(p) || isNaN(μ) || isNaN(σ) {
			return p + μ + σ
		}
		if σ <= 0 {
			return NaN
		}
		if σ == 0 {
			return μ
		}

		// p := logit(p) = log( p / (1-p) )
		p = log(p / (1 - p))
		return μ + σ*p
	}
}

// LogisticQtlFor returns the inverse of the CDF (quantile) of the Logistic distribution, for given probability.
func LogisticQtlFor(μ, σ, p float64) float64 {
	qtl := LogisticQtl(μ, σ)
	return qtl(p)
}

// LogisticNext returns random number drawn from the Logistic distribution. 
func LogisticNext(μ, σ float64) float64 {
	p := UniformNext(0, 1)
	return LogisticQtlFor(μ, σ, p)
}

// Logistic returns the random number generator with  Logistic distribution. 
func Logistic(μ, σ float64) func() float64 {
	return func() float64 { return LogisticNext(μ, σ) }
}

// LogisticMean returns the mean of the Logistic distribution. 
func LogisticMean(μ, σ float64) float64 {
	return μ
}

// LogisticMode returns the mode of the Logistic distribution. 
func LogisticMode(μ, σ float64) float64 {
	return μ
}

// LogisticMedian returns the median of the Logistic distribution. 
func LogisticMedian(μ, σ float64) float64 {
	return μ
}

// LogisticVar returns the variance of the Logistic distribution. 
func LogisticVar(μ, σ float64) float64 {
	return σ * σ * π * π / 3
}

// LogisticStd returns the standard deviation of the Logistic distribution. 
func LogisticStd(μ, σ float64) float64 {
	return σ * π / sqrt(3)
}

// LogisticSkew returns the skewness of the Logistic distribution. 
func LogisticSkew(μ, σ float64) float64 {
	return 0
}

// LogisticExKurt returns the excess kurtosis of the Logistic distribution. 
func LogisticExKurt(μ, σ float64) float64 {
	return 6.0 / 5.0
}

// LogisticMGF returns the moment-generating function of the Logistic distribution. 
func LogisticMGF(μ, σ, t float64) float64 {
	return exp(μ*t) * B(0, 2) // TO BE CHECKED
}
