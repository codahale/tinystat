// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Tapered Pareto distribution. 
//
// Parameters: 
// θ > 0.0 (scale) "a"
// α > 0.0 (shape) "lambda"
// taper > 0.0 (tapering) 
//
// Support: 
//  x >= θ 
// inspired by R:PtProcess

import (
	"math/rand"
)

// ParetoTapChkParams checks parameters of the Tapered Pareto distribution. 
func ParetoTapChkParams(θ, α, taper float64) bool {
	ok := true
	if α <= 0 || θ <= 0 || taper <= 0 {
		ok = false
	}
	return ok
}

// ParetoTapChkSupport checks support of the Tapered Pareto distribution. 
func ParetoTapChkSupport(x float64) bool {
	ok := true
	if x < 0 {
		ok = false
	}
	return ok
}

// ParetoTapPDF returns the PDF of the Tapered Pareto distribution. 
func ParetoTapPDF(θ, α, taper float64) func(x float64) float64 {
	return func(x float64) float64 {
		return (α/x + 1/taper) * pow((θ/x), α) * exp((θ-x)/taper)
	}
}

// ParetoTapPDFAt returns the value of PDF of Tapered Pareto distribution at x. 
func ParetoTapPDFAt(θ, α, taper, x float64) float64 {
	pdf := ParetoTapPDF(θ, α, taper)
	return pdf(x)
}

// ParetoTapCDF returns the CDF of the Tapered Pareto distribution. 
func ParetoTapCDF(θ, α, taper float64) func(x float64) float64 {
	return func(x float64) float64 {
		return 1 - pow(θ/x, α)*exp((θ-x)/taper)
	}
}

// ParetoTapCDFAt returns the value of CDF of the Tapered Pareto distribution, at x. 
func ParetoTapCDFAt(θ, α, taper, x float64) float64 {
	cdf := ParetoTapCDF(θ, α, taper)
	return cdf(x)
}

// ParetoTapQtl returns the inverse of the CDF (quantile) of the Tapered Pareto distribution. 
func ParetoTapQtl(θ, α, taper float64) func(p float64) float64 {
	return func(p float64) float64 {
		tol := 1e-8
		x := θ + 1

		// solve using Newton-Raphson method
		for {
			delta := (ParetoTapCDFAt(θ, α, taper, x) - p) / ParetoTapPDFAt(θ, α, taper, x)
			x -= delta
			if abs(delta) < tol {
				break
			}
		}
		return x
	}
}

// ParetoTapQtlFor returns the inverse of the CDF (quantile) of the Tapered Pareto distribution, for given probability.
func ParetoTapQtlFor(θ, α, taper, p float64) float64 {
	cdf := ParetoTapQtl(θ, α, taper)
	return cdf(p)
}

// ParetoTapNext returns random number drawn from the Tapered Pareto distribution. 
func ParetoTapNext(θ, α, taper float64) float64 {
	qtl := ParetoTapQtl(θ, α, taper)
	p := rand.Float64()
	return qtl(p)
}

// ParetoTap returns the random number generator with  Tapered Pareto distribution. 
func ParetoTap(θ, α, taper float64) func() float64 {
	return func() float64 { return ParetoTapNext(θ, α, taper) }
}

/*

func ParetoTapMean(θ, α, taper float64) float64 {
}

func ParetoTapMedian(θ, α, taper float64) float64 {
}

func ParetoTapMode(θ, α, taper float64) float64 {
}

func ParetoTapVar(θ, α, taper float64) float64 {
}

*/
