// Copyright 2012 - 2013 The Probab Authors. All rights reserved. See the LICENSE file.

// Moments for the data vector.

package stat

// Mean returns the mean of the data vector.
func Mean(x []float64) float64 {
	μ := sum(x)
	μ /= float64(len(x))
	return μ
}

// Skew returns skewness of the data vector.
func Skew(x []float64) float64 {
	n := float64(len(x))
	d3 := diffMean(x)
	for i, val := range d3 {
		d3[i] = val * val * val
	}
	d2 := diffMean(x)
	for i, val := range d2 {
		d2[i] = val * val
	}
	return (sum(d3) / n) / pow((sum(d2)/n), 1.5)
}

// Kurt returns the estimator of Pearson’s measure of kurtosis of the data vector. This is NOT the Excess Kustosis!
func Kurt(x []float64) float64 {
	n := float64(len(x))
	d4 := diffMean(x)
	for i, val := range d4 {
		d4[i] = val * val * val * val
	}
	d2 := diffMean(x)
	for i, val := range d2 {
		d2[i] = val * val
	}
	sum2 := sum(d2)
	return n * sum(d4) / (sum2 * sum2)
}

// Moment returns moment of specified order of the data vector.
func moment(x []float64, order int, central, absolute bool) float64 {
	n := float64(len(x))
	if order < 1 {
		panic(" order < 1")
	}
	if central {
		vCent(x)
	}
	if absolute {
		vAbs(x)
	}
	vPow(x, float64(order))
	return sum(x) / n
}

// Geary returns an estimator of Geary’s measure of kurtosis.
func Geary(x []float64) float64 {
	// The Geary’s kurtosis is computed by dividing average difference between observation 
	// and the mean by standard deviation of the sample.
	n := float64(len(x))
	d2 := diffMean(x)
	for i, val := range d2 {
		d2[i] = val * val
	}
	rho := sqrt(sum(d2) / n)
	vCent(x)
	vAbs(x)
	tau := sum(x) / n
	return tau / rho
}
