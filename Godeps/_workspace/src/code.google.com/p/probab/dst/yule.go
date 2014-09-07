// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Yule–Simon distribution. 
// Yule, G. U. (1925). "A Mathematical Theory of Evolution, based on the Conclusions of Dr. J. C. Willis, F.R.S". Philosophical Transactions of the Royal Society of London, Ser. B 213 (402–410): 21–87. doi:10.1098/rstb.1925.0002
// Simon, H. A. (1955). "On a class of skew distribution functions". Biometrika 42 (3–4): 425–440. doi:10.1093/biomet/42.3-4.425
//
// Parameters: 
// a > 0	 	shape (real)
//
// Support: 
// k ∈ {1, 2, ... }

// YulePMF returns the PMF of the Yule–Simon distribution. 
func YulePMF(a float64) func(k int64) float64 {
	return func(k int64) float64 {
		p := a * B(a+1, float64(k))
		return p
	}
}

// YulePMFAt returns the value of PMF of Yule–Simon distribution at k. 
func YulePMFAt(a float64, k int64) float64 {
	pmf := YulePMF(a)
	return pmf(k)
}

// YuleCDF returns the CDF of the Yule–Simon distribution. 
func YuleCDF(a float64) func(k int64) float64 {
	return func(k int64) float64 {
		kk := float64(k)
		p := 1 - kk*B(kk, a+1)
		return p
	}
}

// YuleCDFAt returns the value of CDF of the Yule–Simon distribution, at x. 
func YuleCDFAt(a float64, k int64) float64 {
	cdf := YuleCDF(a)
	return cdf(k)
}

// YuleNext returns random number drawn from the Yule–Simon distribution. 
func YuleNext(a float64) (k int64) {
	// Devroye 1986: 553.
	// Devroye, L. 1986: Non-Uniform Random Variate Generation. Springer-Verlag, New York. ISBN 0-387-96305-7.
	e1 := ExponentialNext(2)
	e2 := ExponentialNext(2)
	k = int64(ceil(-e1 / (log(1 - exp(-e2/(a-1))))))
	return
}

// Yule returns the random number generator with  Yule–Simon distribution. 
func Yule(a float64) func() int64 {
	return func() int64 { return YuleNext(a) }
}

// YuleMean returns the mean of the Yule–Simon distribution. 
func YuleMean(a float64) float64 {
	if a <= 1 {
		return NaN
	}
	return a / (a - 1)
}

// YuleMode returns the mode of the Yule–Simon distribution. 
func YuleMode(a float64) float64 {
	return 1.00
}

// YuleVar returns the variance of the Yule–Simon distribution. 
func YuleVar(a float64) float64 {
	if a <= 2 {
		return NaN
	}
	return a * a / ((a - 1) * (a - 1) * (a - 2))
}

// YuleStd returns the standard deviation of the Yule–Simon distribution. 
func YuleStd(a float64) float64 {
	if a <= 2 {
		return NaN
	}
	return a / ((a - 1) * sqrt(a-2))
}

// YuleSkew returns the skewness of the Yule–Simon distribution. 
func YuleSkew(a float64) float64 {
	if a <= 3 {
		return NaN
	}
	return ((a + 1) * (a + 1) * sqrt(a-2)) / ((a - 3) * a)
}

// YuleExKurt returns the excess kurtosis of the Yule–Simon distribution. 
func YuleExKurt(a float64) float64 {
	if a <= 4 {
		return NaN
	}
	return a + 3 + (11*a*a*a-49*a-22)/((a-4)*(a-3)*a)
}
