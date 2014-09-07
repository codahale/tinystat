// Copyright 2012 - 2013 The Probab Authors. All rights reserved. See the LICENSE file.

package stat

// Bonett-Seier test of Geary’s kurtosis.
// Ref.: 

import (
	"code.google.com/p/probab/dst"
	"sort"
)

// Bonett performs Bonett-Seier test of Geary’s measure of kurtosis for normally distributed data vector.
func Bonett(x []float64, alternative int) (kurt, z, pVal float64) {
	// Arguments: 
	// x - vector of observations
	// alternative - 0 = "twoSided", 1 = "less", 2 = "greater"
	//
	// Details: 
	// Under the hypothesis of normality, data should have Geary’s kurtosis equal to sqrt(2/pi) (0.7979).
	// This test has such null hypothesis and is useful to detect a significant difference of Geary’s kurtosis
	// in normally distributed data.
	//
	// Returns: 
	// kurt - kurtosis estimator 
	// z - its transformation
	// pVal - the p-value for the test.

	const (
		twoSided = iota
		less
		greater
	)

	sort.Float64s(x)
	n := float64(len(x))
	dm := diffMean(x)

	d2 := make([]float64, len(dm))
	for i, val := range dm {
		d2[i] = val * val
	}

	adm := make([]float64, len(dm))
	for i, val := range dm {
		adm[i] = abs(val)
	}

	rho := sqrt(sum(d2) / n)
	kurt = sum(adm) / n
	omega := 13.29 * (log(rho) - log(kurt))
	z = sqrt(n+2) * (omega - 3) / 3.54
	pVal = 1 - dst.NormalCDFAt(0, 1, z)

	switch alternative {
	case twoSided:
		pVal = 2 * pVal
		if pVal > 1 {
			pVal = 2 - pVal
		}
	case less: // do nothing
	case greater:
		pVal = 1 - pVal
	}

	return kurt, z, pVal
}
