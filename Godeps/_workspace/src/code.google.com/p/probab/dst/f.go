// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// F-distribution, alias Fisher-Snedecor distribution

// FPDF returns the PDF of the F distribution. 
func FPDF(d1, d2 int64) func(x float64) float64 {
	df1 := float64(d1)
	df2 := float64(d2)
	normalization := 1 / B(df1/2, df2/2)
	return func(x float64) float64 {
		return normalization * sqrt(pow(df1*x, df1)*pow(df2, df2)/pow(df1*x+df2, df1+df2)) / x
	}
}

// FLnPDF returns the natural logarithm of the PDF of the F distribution. 
func FLnPDF(d1, d2 int64) func(x float64) float64 {
	df1 := float64(d1)
	df2 := float64(d2)
	normalization := -logB(df1/2, df2/2)
	return func(x float64) float64 {
		return normalization + log(df1*x)*df1/2 + log(df2)*df2/2 - log(df1*x+df2)*(df1+df2)/2 - log(x)
	}
}

// FPDFAt returns the value of PDF of F distribution at x. 
func FPDFAt(d1, d2 int64, x float64) float64 {
	pdf := FPDF(d1, d2)
	return pdf(x)
}

// FCDF returns the CDF of the F distribution. 
func FCDF(d1, d2 int64) func(x float64) float64 {
	return func(x float64) float64 {
		df1 := float64(d1)
		df2 := float64(d2)
		y := df1 * x / (df1*x + df2)
		return iBr(df1/2.0, df2/2.0, y)
	}
}

// FCDFAt returns the value of CDF of the F distribution, at x. 
func FCDFAt(d1, d2 int64, x float64) float64 {
	cdf := FCDF(d1, d2)
	return cdf(x)
}

// FQtl returns the inverse of the CDF (quantile) of the F distribution. 
func FQtl(d1, d2 int64) func(p float64) float64 {
	df1 := float64(d1)
	df2 := float64(d2)
	return func(p float64) float64 {
		if p < 0.0 {
			return NaN
		}
		if p > 1.0 {
			return NaN
		}
		if df1 < 1.0 {
			return NaN
		}
		if df2 < 1.0 {
			return NaN
		}
		return ((1/BetaQtlFor(df2/2, df1/2, 1-p) - 1) * df2 / df1)
	}
}

// FQtlFor returns the inverse of the CDF (quantile) of the F distribution, for given probability.
func FQtlFor(d1, d2 int64, p float64) float64 {
	cdf := FQtl(d1, d2)
	return cdf(p)
}

// FNext returns random number drawn from the F distribution. 
func FNext(d1, d2 int64) float64 {
	df1 := float64(d1)
	df2 := float64(d2)
	return ChiSquareNext(d1) * df2 / (ChiSquareNext(d2) * df1)
}

// F returns the random number generator with  F distribution. 
func F(d1, d2 int64) func() float64 {
	return func() float64 {
		return FNext(d1, d2)
	}
}

// FMean returns the mean of the F distribution. 
func FMean(d1, d2 int64) float64 {
	if d2 <= 2 {
		return NaN
	}
	df2 := float64(d2)
	return df2 / (df2 - 2)
}

// FMode returns the mode of the F distribution. 
func FMode(d1, d2 int64) float64 {
	if d1 <= 2 {
		return NaN
	}
	df1 := float64(d1)
	df2 := float64(d2)
	return ((df1 - 2) / df1) * (df2 / (df2 + 2))
}

// FVar returns the variance of the F distribution. 
func FVar(d1, d2 int64) float64 {
	if d2 <= 4 {
		return NaN
	}
	df1 := float64(d1)
	df2 := float64(d2)
	return 2 * df2 * df2 * (df1 + df2 - 2) / (df1 * (df2 - 2) * (df2 - 2) * (df2 - 4))
}

// FStd returns the standard deviation of the F distribution. 
func FStd(d1, d2 int64) float64 {
	if d2 <= 4 {
		return NaN
	}
	df1 := float64(d1)
	df2 := float64(d2)
	v := 2 * df2 * df2 * (df1 + df2 - 2) / (df1 * (df2 - 2) * (df2 - 2) * (df2 - 4))
	return sqrt(v)
}

// FSkew returns the skewness of the F distribution. 
func FSkew(d1, d2 int64) float64 {
	if d2 <= 6 {
		return NaN
	}
	df1 := float64(d1)
	df2 := float64(d2)
	return (2*df1 + df2 - 2) * sqrt(8*(df2-4)) / (df2 - 6) * sqrt(df1*(df1+df2-2))
}

// FExKurt returns the excess kurtosis of the F distribution. 
func FExKurt(d1, d2 int64) float64 {
	if d2 <= 8 {
		return NaN
	}
	df1 := float64(d1)
	df2 := float64(d2)
	return 12 * (df1*(5*df2-22)*(df1+df2-2) + (df2-4)*(df2-2)*(df2-2)) / (df1 * (df2 - 6) * (df2 - 8) * (df1 + df2 - 2))
}
