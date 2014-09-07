// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Beta distribution reparametrized using mean (μ) and sample size (ν). 
// Kruschke, J. K. (2011). Doing Bayesian data analysis: A tutorial with R and BUGS. p. 83: Academic Press / Elsevier. ISBN 978-0123814852.

// BetaμνPDF returns the PDF of the Beta distribution reparametrized using mean and sample size. 
func BetaμνPDF(μ, ν float64) func(x float64) float64 {
	α := μ * ν
	β := (1 - μ) * ν
	return BetaPDF(α, β)
}

// BetaμνLnPDF returns the natural logarithm of the PDF of the Beta distribution reparametrized using mean and sample size. 
func BetaμνLnPDF(μ, ν float64) func(x float64) float64 {
	α := μ * ν
	β := (1 - μ) * ν
	return BetaLnPDF(α, β)
}

// BetaμνNext returns random number drawn from the  Beta distribution reparametrized using mean and sample size. 
func BetaμνNext(μ, ν float64) float64 {
	α := μ * ν
	β := (1 - μ) * ν
	if ν <= 0 {
		return NaN
	}
	return BetaNext(α, β)
}

// Betaμν returns the random number generator with  Beta distribution reparametrized using mean and sample size. 
func Betaμν(μ, ν float64) func() float64 {
	α := μ * ν
	β := (1 - μ) * ν
	return func() float64 { return BetaNext(α, β) }
}

// BetaμνPDFAt returns the value of PDF of Beta distribution at x. 
func BetaμνPDFAt(μ, ν, x float64) float64 {
	pdf := BetaμνPDF(μ, ν)
	return pdf(x)
}

// BetaμνCDF returns the CDF of the Beta distribution reparametrized using mean and sample size. 
func BetaμνCDF(μ, ν float64) func(x float64) float64 {
	α := μ * ν
	β := (1 - μ) * ν
	return BetaCDF(α, β)
}

// BetaμνCDFAt returns the value of CDF of the Beta distribution reparametrized using mean and sample size, at x. 
func BetaμνCDFAt(μ, ν, x float64) float64 {
	cdf := BetaCDF(μ, ν)
	return cdf(x)
}

// BetaμνQtl returns the inverse of the CDF (quantile) of the Beta distribution reparametrized using mean and sample size. 
func BetaμνQtl(μ, ν float64) func(p float64) float64 {
	// p: probability for which the quantile is evaluated
	α := μ * ν
	β := (1 - μ) * ν
	return BetaQtl(α, β)
}

// BetaμνQtlFor returns the inverse of the CDF (quantile) of the Beta distribution reparametrized using mean and sample size, for a given probability.
func BetaμνQtlFor(μ, ν, p float64) float64 {
	cdf := BetaμνQtl(μ, ν)
	return cdf(p)
}
