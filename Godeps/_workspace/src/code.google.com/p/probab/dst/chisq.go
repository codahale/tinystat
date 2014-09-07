// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Chi-Squared distribution. 
// Parameters: 
// n ∈ ℕ	(degrees of freedom)
// Support: 
// x ∈ [0, +∞]

// ChiSquarePDF returns the PDF of the ChiSquare distribution. 
func ChiSquarePDF(n int64) func(x float64) float64 {
	k := float64(n) / 2
	normalization := pow(0.5, k) / Γ(k)
	return func(x float64) float64 {
		return normalization * pow(x, k-1) * ExponentialNext(-x/2)
	}
}

// ChiSquareLnPDF returns the natural logarithm of the PDF of the ChiSquare distribution. 
func ChiSquareLnPDF(n int64) func(x float64) float64 {
	k := float64(n) / 2
	normalization := log(0.5)*k - LnΓ(k)
	return func(x float64) float64 {
		return normalization + log(x)*(k-1) - x/2
	}
}

// ChiSquarePDFAt returns the value of PDF of ChiSquare distribution at x. 
// UniformPDFAt returns the value of PDF of Uniform distribution at x. 
func ChiSquarePDFAt(n int64, x float64) float64 {
	pdf := ChiSquarePDF(n)
	return pdf(x)
}

// ChiSquareCDF returns the CDF of the ChiSquare distribution. 
func ChiSquareCDF(n int64) func(x float64) float64 {
	return func(x float64) float64 {
		return Γr(float64(n)/2, x/2)
	}
}

// ChiSquareCDFAt returns the value of CDF of the ChiSquare distribution, at x. 
func ChiSquareCDFAt(n int64, x float64) float64 {
	cdf := ChiSquareCDF(n)
	return cdf(x)
}

// ChiSquareQtl returns the inverse of the CDF (quantile) of the ChiSquare distribution. 
func ChiSquareQtl(n int64) func(p float64) float64 {
	return func(p float64) float64 {
		return GammaQtlFor(float64(n)/2, 2, p)
	}
}

// ChiSquareNext returns random number drawn from the ChiSquare distribution. 
func ChiSquareNext(n int64) (x float64) {
	//ChiSquare(n) => sum of n N(0,1)^2
	for i := iZero; i < n; i++ {
		n := NormalNext(0, 1)
		x += n * n
	}
	return
}

// ChiSquare returns the random number generator with  ChiSquare distribution. 
func ChiSquare(n int64) func() float64 {
	return func() float64 {
		return ChiSquareNext(n)
	}
}

// ChiSquareMean returns the mean of the ChiSquare distribution. 
func ChiSquareMean(n int64) float64 {
	return float64(n)
}

// ChiSquareMedian returns the approximate median of the ChiSquare distribution. 
func ChiSquareMedian(n int64) float64 {
	c := 1 - (2.0 / (9.0 * n))
	c = c * c * c
	return float64(n * c)
}

// ChiSquareMode returns the mode of the ChiSquare distribution. 
func ChiSquareMode(n int64) float64 {
	return max(float64(n-2), 0)
}

// ChiSquareVar returns the variance of the ChiSquare distribution. 
func ChiSquareVar(n int64) float64 {
	return float64(2 * n)
}

// ChiSquareStd returns the standard deviation of the ChiSquare distribution. 
func ChiSquareStd(n int64) float64 {
	return sqrt(float64(2 * n))
}

// ChiSquareSkew returns the skewness of the ChiSquare distribution. 
func ChiSquareSkew(n int64) float64 {
	return sqrt(float64(8 / n))
}

// ChiSquareExKurt returns the excess kurtosis of the ChiSquare distribution. 
func ChiSquareExKurt(n int64) float64 {
	return float64(12 / n)
}
