// Copyright 2012 - 2013 The Probab Authors. All rights reserved. See the LICENSE file.

package stat

// Anscombe-Glynn test of kurtosis for normally distributed data.
// Ref.: Anscombe & Glynn (1983).

import (
	"code.google.com/p/probab/dst"
	"sort"
)

// Anscombe performs Anscombe-Glynn test of kurtosis for normally distributed data vector.
func Anscombe(x []float64, alternative int) (kurt, z, pVal float64) {
	// Arguments: 
	// x - vector of observations
	// alternative - 0 = "twoSided", 1 = "less", 2 = "greater"
	//
	// Details: 
	// Under the hypothesis of normality, data should have kurtosis equal to 3. This test has such null
	// hypothesis and is useful to detect a significant difference of kurtosis in normally distributed data.
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
	d4 := make([]float64, len(dm))
	for i, val := range dm {
		d4[i] = val * val * val * val
	}

	d2 := make([]float64, len(dm))
	for i, val := range dm {
		d2[i] = val * val
	}

	//b <- n*sum( (x-mean(x))^4 )/(sum( (x-mean(x))^2 )^2);
	sum2 := sum(d2)
	kurt = n * sum(d4) / (sum2 * sum2)

	eb2 := 3 * (n - 1) / (n + 1)
	vb2 := 24 * n * (n - 2) * (n - 3) / ((n + 1) * (n + 1) * (n + 3) * (n + 5))
	m3 := (6 * (n*n - 5*n + 2) / ((n + 7) * (n + 9))) * sqrt((6*(n+3)*(n+5))/(n*(n-2)*(n-3)))
	a := 6 + (8/m3)*(2/m3+sqrt(1+4/(m3*m3)))
	xx := (kurt - eb2) / sqrt(vb2)
	z0 := (1 - 2/a) / (1 + xx*sqrt(2/(a-4)))
	z = (1 - 2/(9*a) - pow(z0, 1.0/3.0)) / sqrt(2/(9*a))
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
