// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Zeta distribution. 
// The zeta distribution is equivalent to the Zipf distribution for infinite N. 
//
// Parameters: 
// s > 1.0	 	(real)
//
// Support: 
// k > 0		(integer)

import (
	"math/rand"
)

// ZetaPMF returns the PMF of the Zeta distribution. 
func ZetaPMF(s float64) func(k int64) float64 {
	return func(k int64) float64 {
		t1 := 1 / pow(float64(k), s)
		t2 := ζ(s)
		p := t1 / t2
		return p
	}
}

// ZetaPMFAt returns the value of PMF of Zeta distribution at k. 
func ZetaPMFAt(s float64, k int64) float64 {
	pmf := ZetaPMF(s)
	return pmf(k)
}

// ZetaCDF returns the CDF of the Zeta distribution. 
func ZetaCDF(s float64) func(k int64) float64 {
	return func(k int64) float64 {
		t1 := hNum(k, s)
		t2 := ζ(s)
		p := t1 / t2
		return p
	}
}

// ZetaCDFAt returns the value of CDF of the Zeta distribution, at x. 
func ZetaCDFAt(s float64, k int64) float64 {
	pmf := ZetaCDF(s)
	return pmf(k)
}

// ZetaNext returns random number drawn from the Zeta distribution. 
func ZetaNext(s float64) (k int64) {
	// Devroye 1986: 550. Called "Zipf distribution" there.
	// Devroye, L. 1986: Non-Uniform Random Variate Generation. Springer-Verlag, New York. ISBN 0-387-96305-7.
	var x float64
	b := pow(2.0, s-1.0)
	for {
		u := rand.Float64()
		v := rand.Float64()
		x = floor(pow(u, -1/(s-1)))
		t := pow(1+1.0/x, s-1)
		delta := v * x * (t - 1.0) / (b - 1.0)
		if delta <= (t / b) {
			break
		}
	}
	k = int64(x)
	return
}

// Zeta returns the random number generator with  Zeta distribution. 
func Zeta(s float64) func() int64 {
	return func() int64 { return ZetaNext(s) }
}

// ZetaMean returns the mean of the Zeta distribution. 
func ZetaMean(s float64) float64 {
	if s <= 2 {
		return NaN
	}
	t1 := ζ(s - 1)
	t2 := ζ(s)
	return t1 / t2
}

// ZetaMode returns the mode of the Zeta distribution. 
func ZetaMode() float64 {
	return 1
}

// ZetaVar returns the variance of the Zeta distribution. 
func ZetaVar(s float64) float64 {
	if s <= 3 {
		return NaN
	}
	t1 := ζ(s)
	t2 := ζ(s - 2)
	t3 := ζ((s - 1) * (s - 1))
	t4 := ζ(s * s)
	return (t1*t2 - t3) / t4
}
