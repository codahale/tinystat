// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Log-normal () distribution is a continuous probability distribution of a random variable whose logarithm is normally distributed. If X is a random variable with a normal distribution, then Y = exp(X) has a log-normal distribution; likewise, if Y is log-normally distributed, then X = log(Y) has a normal distribution. 
//
// Parameters: 
// μ ∈ R		(location)
// σ > 0		standard deviation  (scale)
//
// Support: 
// x ∈ R

import (

//	"math/rand"
)

// LogNormalPDF returns the PDF of the LogNormal distribution. 
func LogNormalPDF(μ, σ float64) func(x float64) float64 {
	normalogormalizer := 0.3989422804014327 / σ
	return func(x float64) float64 { return normalogormalizer * exp(-1*(log(x)-μ)*(log(x)-μ)/(2*σ*σ)) / x }
}

// LogNormalPDFAt returns the value of PDF of LogNormal distribution at x. 
func LogNormalPDFAt(μ, σ, x float64) float64 {
	pdf := LogNormalPDF(μ, σ)
	return pdf(x)
}

// LogNormalCDF returns the CDF of the LogNormal distribution. 
func LogNormalCDF(μ, σ float64) func(x float64) float64 {
	return func(x float64) float64 { return ((1.0 / 2.0) * (1 + erf((log(x)-μ)/(σ*sqrt2)))) }
}

// LogNormalCDFAt returns the value of CDF of the LogNormal distribution, at x. 
func LogNormalCDFAt(μ, σ, x float64) float64 {
	cdf := LogNormalCDF(μ, σ)
	return cdf(x)
}

// LogNormalQtl returns the inverse of the CDF (quantile) of the LogNormal distribution. 
func LogNormalQtl(μ, σ float64) func(p float64) float64 {
	return func(p float64) float64 {
		return exp(σ*ZQtlFor(p) + μ)
	}
}

// LogNormalQtlFor returns the inverse of the CDF (quantile) of the LogNormal distribution, for given probability.
func LogNormalQtlFor(μ, σ, p float64) float64 {
	qtl := LogNormalQtl(μ, σ)
	return qtl(p)
}

// LogNormalNext returns random number drawn from the LogNormal distribution. 
func LogNormalNext(μ, σ float64) float64 { return exp(NormalNext(μ, σ)) }

// LogNormal returns the random number generator with  LogNormal distribution. 
func LogNormal(μ, σ float64) func() float64 {
	return func() float64 { return LogNormalNext(μ, σ) }
}

// LogNormalMean returns the mean of the LogNormal distribution. 
func LogNormalMean(μ, σ float64) float64 {
	return exp(μ + σ*σ/2)
}

// LogNormalMode returns the mode of the LogNormal distribution. 
func LogNormalMode(μ, σ float64) float64 {
	return exp(μ - σ*σ)
}

// LogNormalMedian returns the median of the LogNormal distribution. 
func LogNormalMedian(μ, σ float64) float64 {
	return exp(μ)
}

// LogNormalVar returns the variance of the LogNormal distribution. 
func LogNormalVar(μ, σ float64) float64 {
	return exp(σ*σ) - 1*exp(2*μ+σ*σ)
}

// LogNormalStd returns the standard deviation of the LogNormal distribution. 
func LogNormalStd(μ, σ float64) float64 {
	return sqrt(LogNormalVar(μ, σ))
}

// LogNormalSkew returns the skewness of the LogNormal distribution. 
func LogNormalSkew(μ, σ float64) float64 {
	return exp(σ*σ) + 2*sqrt(exp(σ*σ)-1)
}

// LogNormalExKurt returns the excess kurtosis of the LogNormal distribution. 
func LogNormalExKurt(μ, σ float64) float64 {
	return exp(4*σ*σ) + 2*exp(3*σ*σ) + 3*exp(2*σ*σ) - 6
}
