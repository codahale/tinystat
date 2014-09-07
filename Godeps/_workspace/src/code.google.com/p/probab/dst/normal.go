// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Normal  (or Gaussian, or Gauss-Laplace) distribution. 
// A continuous probability distribution, defined on the entire real line, that has a bell-shaped probability density function, known as the Gaussian function. 
//
// Parameters: 
// μ ∈ R		(location)
// σ > 0		standard deviation  (scale)
//
// Support: 
// x ∈ R

import (
	"math/rand"
)

func rateval(a []float64, na int64, b []float64, nb int64, x float64) float64 {
	var (
		i, j    int64
		u, v, r float64
	)

	u = a[na-1]

	for i = na - 1; i > 0; i-- {
		u = x*u + a[i-1]
	}

	v = b[nb-1]

	for j = nb - 1; j > 0; j-- {
		v = x*v + b[j-1]
	}

	r = u / v

	return r
}

func small(q float64) float64 {
	var (
		a = []float64{3.387132872796366608, 133.14166789178437745, 1971.5909503065514427, 13731.693765509461125, 45921.953931549871457, 67265.770927008700853, 33430.575583588128105, 2509.0809287301226727}

		b = []float64{1.0, 42.313330701600911252, 687.1870074920579083, 5394.1960214247511077, 21213.794301586595867, 39307.89580009271061,
			28729.085735721942674, 5226.495278852854561}

		r, x float64
	)
	r = 0.180625 - q*q
	x = q * rateval(a, 8, b, 8, r)

	return x
}

func intermediate(r float64) float64 {
	var (
		a = []float64{1.42343711074968357734, 4.6303378461565452959, 5.7694972214606914055, 3.64784832476320460504,
			1.27045825245236838258, 0.24178072517745061177, 0.0227238449892691845833, 7.7454501427834140764e-4}

		b = []float64{1.0, 2.05319162663775882187, 1.6763848301838038494, 0.68976733498510000455, 0.14810397642748007459, 0.0151986665636164571966,
			5.475938084995344946e-4, 1.05075007164441684324e-9}

		x float64
	)
	x = rateval(a, 8, b, 8, (r - 1.6))
	return x
}

func tail(r float64) float64 {
	var (
		a = []float64{6.6579046435011037772, 5.4637849111641143699, 1.7848265399172913358, 0.29656057182850489123,
			0.026532189526576123093, 0.0012426609473880784386, 2.71155556874348757815e-5, 2.01033439929228813265e-7}

		b = []float64{1.0, 0.59983220655588793769, 0.13692988092273580531, 0.0148753612908506148525, 7.868691311456132591e-4, 1.8463183175100546818e-5, 1.4215117583164458887e-7, 2.04426310338993978564e-15}

		x float64
	)
	x = rateval(a, 8, b, 8, (r - 5.0))
	return x
}

// NormalPDF returns the PDF of the Normal distribution. 
func NormalPDF(μ, σ float64) func(x float64) float64 {
	normalnormalizer := 0.3989422804014327 / σ
	return func(x float64) float64 { return normalnormalizer * exp(-1*(x-μ)*(x-μ)/(2*σ*σ)) }
}

// NormalLnPDF returns the natural logarithm of the PDF of the Normal distribution. 
func NormalLnPDF(μ, σ float64) func(x float64) float64 {
	lnnormalnormalizer := -0.91893853320467267 - log(σ)
	return func(x float64) float64 { return lnnormalnormalizer - (x-μ)*(x-μ)/(2*σ*σ) }
}

// NormalPDFAt returns the value of PDF of Normal distribution at x. 
func NormalPDFAt(μ, σ, x float64) float64 {
	pdf := NormalPDF(μ, σ)
	return pdf(x)
}

// NormalCDF returns the CDF of the Normal distribution. 
func NormalCDF(μ, σ float64) func(x float64) float64 {
	return func(x float64) float64 { return ((1.0 / 2.0) * (1 + erf((x-μ)/(σ*sqrt2)))) }
}

// NormalCDFAt returns the value of CDF of the Normal distribution, at x. 
func NormalCDFAt(μ, σ, x float64) float64 {
	cdf := NormalCDF(μ, σ)
	return cdf(x)
}

// NormalQtl returns the inverse of the CDF (quantile) of the Normal distribution. 
func NormalQtl(μ, σ float64) func(p float64) float64 {
	return func(p float64) float64 {
		return σ*ZQtlFor(p) + μ
	}
}

// NormalQtlFor returns the inverse of the CDF (quantile) of the Normal distribution, for given probability.
func NormalQtlFor(μ, σ, p float64) float64 {
	qtl := NormalQtl(μ, σ)
	return qtl(p)
}

// NormalNext returns random number drawn from the Normal distribution. 
func NormalNext(μ, σ float64) float64 { return rand.NormFloat64()*σ + μ }

// Normal returns the random number generator with  Normal distribution. 
func Normal(μ, σ float64) func() float64 {
	return func() float64 { return NormalNext(μ, σ) }
}

// NormalMean returns the mean of the Normal distribution. 
func NormalMean(μ, σ float64) float64 {
	return μ
}

// NormalMode returns the mode of the Normal distribution. 
func NormalMode(μ, σ float64) float64 {
	return μ
}

// NormalMedian returns the median of the Normal distribution. 
func NormalMedian(μ, σ float64) float64 {
	return μ
}

// NormalVar returns the variance of the Normal distribution. 
func NormalVar(μ, σ float64) float64 {
	return σ * σ
}

// NormalStd returns the standard deviation of the Normal distribution. 
func NormalStd(μ, σ float64) float64 {
	return σ
}

// NormalSkew returns the skewness of the Normal distribution. 
func NormalSkew(μ, σ float64) float64 {
	return 0
}

// NormalExKurt returns the excess kurtosis of the Normal distribution. 
func NormalExKurt(μ, σ float64) float64 {
	return 0
}

// NormalMGF returns the moment-generating function of the Normal distribution. 
func NormalMGF(μ, σ, t float64) float64 {
	return exp(μ*t + σ*σ*t*t/2)
}
