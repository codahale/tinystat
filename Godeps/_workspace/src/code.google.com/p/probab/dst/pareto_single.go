// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Single-parameter  Pareto distribution. 
// Appendix A of Klugman, Panjer & Willmot, Loss Models, Second Edition, Wiley, 2004.
//
// Parameters: 
// α > 0.0		shape (real)
// μ > 0.0		minimum (real)
//
// Support: 
// x >= μ (real)
// inspired by R:actuar

import (
	"math/rand"
)

// ParetoSingChkParams checks parameters of the Single-parameter Pareto distribution. 
func ParetoSingChkParams(α, μ float64) bool {
	ok := true
	if α <= 0 || μ <= 0 {
		ok = false
	}
	return ok
}

// ParetoSingChkSupport checks support of the Single-parameter Pareto distribution. 
func ParetoSingChkSupport(x, μ float64) bool {
	ok := true
	if x < μ {
		ok = false
	}
	return ok
}

// ParetoSingPDF returns the PDF of the Single-parameter Pareto distribution. 
func ParetoSingPDF(α, μ float64) func(x float64) float64 {
	return func(x float64) float64 {
		p := 0.0
		if x >= μ {
			p = exp(log(α) + α*log(μ) - (α+1.0)*log(x))
		}
		return p
	}
}

// ParetoSingPDFAt returns the value of PDF of Single-parameter  Pareto distribution at x. 
func ParetoSingPDFAt(α, μ, x float64) float64 {
	pdf := ParetoSingPDF(α, μ)
	return pdf(x)
}

// ParetoSingCDF returns the CDF of the Single-parameter  Pareto distribution. 
func ParetoSingCDF(α, μ float64) func(x float64) float64 {
	return func(x float64) float64 {
		p := 0.0
		if x > μ {
			p = (0.5 - (pow(μ/x, α)) + 0.5)
		}
		return p
	}
}

// ParetoSingCDFAt returns the value of CDF of the Single-parameter  Pareto distribution, at x. 
func ParetoSingCDFAt(α, μ, x float64) float64 {
	cdf := ParetoSingCDF(α, μ)
	return cdf(x)
}

// ParetoSingQtl returns the inverse of the CDF (quantile) of the Single-parameter  Pareto distribution. 
func ParetoSingQtl(α, μ float64) func(p float64) float64 {
	return func(p float64) float64 {
		return μ / pow((0.5-(p)+0.5), 1.0/α)
	}
}

// ParetoSingQtlFor returns the inverse of the CDF (quantile) of the Single-parameter  Pareto distribution, for given probability.
func ParetoSingQtlFor(α, μ, p float64) float64 {
	cdf := ParetoSingQtl(α, μ)
	return cdf(p)
}

// ParetoSingNext returns random number drawn from the Single-parameter  Pareto distribution. 
func ParetoSingNext(α, μ float64) float64 {
	qtl := ParetoSingQtl(α, μ)
	p := rand.Float64()
	return qtl(p)
}

// ParetoSing returns the random number generator with  Single-parameter  Pareto distribution. 
func ParetoSing(α, μ float64) func() float64 {
	return func() float64 { return ParetoSingNext(α, μ) }
}

// ParetoSingMoment returns the n-th moment of the Single-parameter  Pareto distribution. 
func ParetoSingMoment(α, μ float64, order int) float64 {
	o := float64(order)
	if o >= α {
		return posInf
	}

	return α * pow(μ, o) / (α - o)
}

// ParetoSingMean returns the mean of the Single-parameter  Pareto distribution. 
func ParetoSingMean(α, μ float64) float64 {
	return ParetoSingMoment(α, μ, 1)
}

// ParetoSingVar returns the variance of the Single-parameter  Pareto distribution. 
func ParetoSingVar(α, μ float64) float64 {
	return ParetoSingMoment(α, μ, 2)
}

// ParetoSingSkew returns the skewness of the Single-parameter  Pareto distribution. 
func ParetoSingSkew(α, μ float64) float64 {
	return ParetoSingMoment(α, μ, 3)
}

// ParetoSingExKurt returns the excess kurtosis of the Single-parameter  Pareto distribution. 
func ParetoSingExKurt(α, μ float64) float64 {
	return ParetoSingMoment(α, μ, 4)
}
