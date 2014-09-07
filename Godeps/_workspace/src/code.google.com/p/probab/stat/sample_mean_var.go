// Copyright 2012 - 2013 The Probab Authors. All rights reserved. See the LICENSE file.

package stat

// Sample mean and variance (unbiased estimator)

// Sample mean and unbiased (Bessel correction) variance estimates for a data vector.
func SampleMeanVar(x []float64) (μ, σ2 float64) {
	// Arguments: 
	// x - vector of observations
	//
	// Returns: 
	// μ - mean estimator 
	// σ2 - variance estimator 

	var n int
	var m, m2 float64
	μ = 0.0  // sample mean
	σ2 = 0.0 // sample variance unbiased
	m = 0.0
	m2 = 0.0

	for _, val := range x {
		n += 1
		μ += val
		delta := val - m
		m += delta / float64(n)
		m2 += delta * (val - m)
	}

	σ2 = m2 / float64(n-1)
	μ /= float64(len(x))
	return
}
