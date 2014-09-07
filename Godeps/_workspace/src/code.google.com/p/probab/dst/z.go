// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Standard Normal  (or Gaussian, or Gauss-Laplace) distribution. 
// A continuous probability distribution, defined on the entire real line, that has a bell-shaped probability density function, known as the Gaussian function. 
//
// Parameters: none
//
// Support: 
// x ∈ R

import (

//	"math/rand"
)

// ZPDF returns the PDF of the Standard Normal distribution. 
func ZPDF() func(float64) float64 {
	return NormalPDF(0, 1)
}

// ZPDFAt returns the value of PDF of Standard Normal distribution at x. 
func ZPDFAt(x float64) float64 {
	pdf := NormalPDF(0, 1)
	return pdf(x)
}

// ZCDF returns the CDF of the Standard Normal distribution. 
func ZCDF() func(float64) float64 {
	return NormalCDF(0, 1)
}

// ZCDFAt returns the value of CDF of the Standard Normal distribution, at x. 
func ZCDFAt(x float64) float64 {
	cdf := NormalCDF(0, 1)
	return cdf(x)
}

// ZQtl returns the inverse of the CDF (quantile) of the Standard Normal distribution. 
func ZQtl() func(p float64) float64 {
	return func(p float64) float64 {

		var r, x, pp, dp float64

		dp = p - 0.5
		switch {
		case p == 1.0:
			return posInf
		case p == 0.0:
			return negInf
		}
		if abs(dp) <= 0.425 {
			x = small(dp)
			return x
		}
		if p < 0.5 {
			pp = p
		} else {
			pp = 1.0 - p
		}
		r = sqrt(-log(pp))
		if r <= 5.0 {
			x = intermediate(r)
		} else {
			x = tail(r)
		}
		if p < 0.5 {
			return -x
		}
		return x
	}
}

// ZQtlFor returns the inverse of the CDF (quantile) of the Standard Normal distribution, for given probability.
func ZQtlFor(p float64) float64 {
	qtl := ZQtl()
	return qtl(p)
}
