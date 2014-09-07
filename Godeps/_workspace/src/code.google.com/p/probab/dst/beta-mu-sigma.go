// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Beta distribution reparametrized using mean and standard deviation. 

// BetaμσPDF returns the PDF of the Beta distribution reparametrized using mean and standard deviation. 
func BetaμσPDF(μ, σ float64) func(x float64) float64 {
	α := μ * (μ*(1-μ)/(σ*σ) - 1)
	β := (1 - μ) * (μ*(1-μ)/(σ*σ) - 1)
	return BetaPDF(α, β)
}

// BetaμσLnPDF returns the natural logarithm of the PDF of the Beta distribution reparametrized using mean and standard deviation. 
func BetaμσLnPDF(μ, σ float64) func(x float64) float64 {
	α := μ * (μ*(1-μ)/(σ*σ) - 1)
	β := (1 - μ) * (μ*(1-μ)/(σ*σ) - 1)
	return BetaLnPDF(α, β)
}

// BetaμσNext returns random number drawn from the  Beta distribution reparametrized using mean and standard deviation. 
func BetaμσNext(μ, σ float64) float64 {
	α := μ * (μ*(1-μ)/(σ*σ) - 1)
	β := (1 - μ) * (μ*(1-μ)/(σ*σ) - 1)
	return BetaNext(α, β)
}

// Betaμσ returns the random number generator with  Beta distribution reparametrized using mean and standard deviation. 
func Betaμσ(μ, σ float64) func() float64 {
	α := μ * (μ*(1-μ)/(σ*σ) - 1)
	β := (1 - μ) * (μ*(1-μ)/(σ*σ) - 1)
	return func() float64 { return BetaNext(α, β) }
}

// BetaμσPDFAt returns the value of PDF of Beta distribution at x. 
func BetaμσPDFAt(μ, σ, x float64) float64 {
	pdf := BetaμσPDF(μ, σ)
	return pdf(x)
}

// BetaμσCDF returns the CDF of the Beta distribution reparametrized using mean and standard deviation. 
func BetaμσCDF(μ, σ float64) func(x float64) float64 {
	α := μ * (μ*(1-μ)/(σ*σ) - 1)
	β := (1 - μ) * (μ*(1-μ)/(σ*σ) - 1)
	return BetaCDF(α, β)
}

// BetaμσCDFAt returns the value of CDF of the Beta distribution reparametrized using mean and standard deviation, at x. 
func BetaμσCDFAt(μ, σ, x float64) float64 {
	cdf := BetaCDF(μ, σ)
	return cdf(x)
}

// BetaμσQtl returns the inverse of the CDF (quantile) of the Beta distribution reparametrized using mean and standard deviation. 
func BetaμσQtl(μ, σ float64) func(p float64) float64 {
	// p: probability for which the quantile is evaluated
	α := μ * (μ*(1-μ)/(σ*σ) - 1)
	β := (1 - μ) * (μ*(1-μ)/(σ*σ) - 1)
	return BetaQtl(α, β)
}

// BetaμσQtlFor returns the inverse of the CDF (quantile) of the Beta distribution reparametrized using mean and standard deviation, for a given probability.
func BetaμσQtlFor(μ, σ, p float64) float64 {
	cdf := BetaμσQtl(μ, σ)
	return cdf(p)
}
