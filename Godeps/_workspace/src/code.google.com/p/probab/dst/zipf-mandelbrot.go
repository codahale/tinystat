// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Zipf-Mandelbrot distribution. 
// For finite n and q == 0 it reduces to Zipf distribution.
//
// Parameters: 
// n ∈ {1, 2, 3,, ...} 	 (integer)
// q ∈ [0, ∞)	(real)
// s ∈ (0, ∞)	(real)
//
// Support: 
// k ∈ {1, 2, ... , n}

// Zipf-Mandelbrot distribution

import (
	"math/rand"
)

//  ZipfMandelbrotChkParams checks parameters of the Zipf-Mandelbrot  distribution.
func ZipfMandelbrotChkParams(n int64, q, s float64) bool {
	v := false
	if n > 0 && q >= 0 && s > 0 {
		v = true
	}
	return v
}

// ZipfMandelbrotPMF returns the PMF of the Zipf-Mandelbrot distribution. 
func ZipfMandelbrotPMF(n int64, q, s float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := 1 / (pow((float64(k)+q), s) * hNumG(n, q, s))
		return p
	}
}

// ZipfMandelbrotPMFAt returns the value of PMF of Zipf-Mandelbrot distribution at k. 
func ZipfMandelbrotPMFAt(n int64, q, s float64, k int64) float64 {
	pmf := ZipfMandelbrotPMF(n, q, s)
	return pmf(k)
}

// ZipfMandelbrotCDF returns the CDF of the Zipf-Mandelbrot distribution. 
func ZipfMandelbrotCDF(n int64, q, s float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := hNumG(k, q, s) / hNumG(n, q, s)
		return p
	}
}

// ZipfMandelbrotCDFAt returns the value of CDF of the Zipf-Mandelbrot distribution, at x. 
func ZipfMandelbrotCDFAt(n int64, q, s float64, k int64) float64 {
	cdf := ZipfMandelbrotCDF(n, q, s)
	return cdf(k)
}

// Quantile Function for the Zipf-Mandelbrot distribution
func ZipfMandelbrotQtl(n int64, q, s float64) func(p float64) int64 {
	return func(p float64) int64 {
		var k int64
		const kMax = 1e16
		cdf := ZipfMandelbrotCDF(n, q, s)
		if cdf(1) >= p {
			k = 1
		} else {
			for k = 1; cdf(k) < p; k++ {
				if k > kMax {
					panic("not found")
				}
			}
		}
		return k
	}
}

// ZipfMandelbrotNext returns random number drawn from the Zipf-Mandelbrot distribution. 
func ZipfMandelbrotNext(n int64, q, s float64) (k int64) {
	qtl := ZipfMandelbrotQtl(n, q, s)
	p := rand.Float64()
	return qtl(p)
}

// ZipfMandelbrot returns the random number generator with  Zipf-Mandelbrot distribution. 
func ZipfMandelbrot(n int64, q, s float64) func() int64 {
	return func() int64 { return ZipfMandelbrotNext(n, q, s) }
}

// ZipfMandelbrotMean returns the mean of the Zipf-Mandelbrot distribution. 
func ZipfMandelbrotMean(n int64, q, s float64) float64 {
	return hNumG(n, q, s-1)/hNumG(n, q, s) - q
}
