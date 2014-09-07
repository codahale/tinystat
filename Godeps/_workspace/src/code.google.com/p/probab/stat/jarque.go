// Copyright 2012 - 2013 The Probab Authors. All rights reserved. See the LICENSE file.

package stat

// Jarque-Bera test for normality.
// Ref.: Jarque & Bera (1980).

import (
	"code.google.com/p/probab/dst"
)

// Jarque performs performs the Jarque-Bera test on the given data sample to determine if the data are 
// sample drawn from a normal population.
func Jarque(x []float64) (jb, pVal float64) {
	// Arguments: 
	// x - vector of observations
	//
	// Details: 
	// Under the hypothesis of normality, data should be symmetrical (i.e. skewness should be equal to zero) 
	// and have skewness chose to three. The Jarque-Bera statistic is chi-square distributed with two degrees of freedom.
	// Alternative hypothesis is "greater".
	//
	// Returns: 
	// jb - the Jarque-Bera statistic
	// pVal - the p-value for the test.

	n := float64(len(x))
	k := Kurt(x)
	s := Skew(x)
	jb = (n / 6) * (s*s + 0.25*((k-3)*(k-3)))
	pVal = 1 - dst.ChiSquareCDFAt(2, jb)
	return
}
