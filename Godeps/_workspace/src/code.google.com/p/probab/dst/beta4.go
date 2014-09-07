// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// The four-parameter Beta distribution. 
// http://en.wikipedia.org/wiki/Betadistribution#Fourparameters2
// A beta distribution with the two shape parameters α and β is supported on the range [0,1]. 
// It is possible to alter the location and scale of the distribution by introducing two further parameters representing the minimum, a, 
// and maximum c (c > a), values of the distribution, by a linear transformation substituting the non-dimensional variable x 
// in terms of the new variable y (with support [a,c]) and the parameters a and c:
// 		y=x*(c-a)+a
// therefore
//		x=(y-a)/(c-a)
//

// Beta4PDF returns the PDF of the four-parameter Beta distribution. 
func Beta4PDF(α, β, a, c float64) func(y float64) float64 {
	return func(y float64) float64 {
		if a >= c {
			return NaN
		}
		dα := []float64{α, β}
		dirPDF := DirichletPDF(dα)
		x := (y - a) / (c - a)
		if 0 > x || x > 1 {
			return 0
		}
		dx := []float64{x, 1 - x}
		return dirPDF(dx) / (c - a)
	}
}

// Beta4Next returns random number drawn from the  four-parameter Beta distribution. 
func Beta4Next(α, β, a, c float64) float64 {
	if a >= c {
		return NaN
	}
	x := BetaNext(α, β)
	y := x*(c-a) + a
	return y
}

// Beta4 returns the random number generator with  four-parameter Beta distribution. 
func Beta4(α, β, a, c float64) func() float64 {
	return func() float64 { return Beta4Next(α, β, a, c) }
}

// Beta4PDFAt returns the value of PDF of four-parameter Beta distribution at x. 
func Beta4PDFAt(α, β, a, c, x float64) float64 {
	pdf := Beta4PDF(α, β, a, c)
	return pdf(x)
}

// Beta4CDF returns the CDF of the four-parameter Beta distribution. 
func Beta4CDF(α, β, a, c float64) func(y float64) float64 {
	return func(y float64) float64 {
		var res float64
		if a >= c {
			return NaN
		}
		x := (y - a) / (c - a)
		z := exp(LnΓ(α+β) - LnΓ(α) - LnΓ(β) + α*log(x) + β*log(1.0-x))
		switch {
		case x == 0:
			res = 0.0
		case x == 1.0:
			res = 1.0
		case x < (α+1.0)/(α+β+2.0):
			res = z * betaContinuedFraction(α, β, x) / α
		default:
			res = 1.0 - z*betaContinuedFraction(β, α, 1.0-x)/β

		}
		return res / (c - a)
	}
}

// Beta4CDFAt returns the value of CDF of the four-parameter Beta distribution, at x. 
func Beta4CDFAt(α, β, a, c, x float64) float64 {
	cdf := Beta4CDF(α, β, a, c)
	return cdf(x)
}

// Beta4Qtl returns the inverse of the CDF (quantile) of the four-parameter Beta distribution. 
func Beta4Qtl(α, β, a, c float64) func(p float64) float64 {
	// p: probability for which the quantile is evaluated
	return func(p float64) float64 {
		var x float64 = 0
		var a float64 = 0
		var b float64 = 1
		var precision float64 = 1e-9
		if a >= c {
			return NaN
		}
		if p < 0.0 {
			return NaN
		}
		if p > 1.0 {
			return NaN
		}
		if α < 0.0 {
			return NaN
		}
		if β < 0.0 {
			return NaN
		}

		for (b - a) > precision {
			x = (a + b) / 2
			if iBr(α, β, x) > p {
				b = x
			} else {
				a = x
			}
		}

		return x*(c-a) + a
	}
}

// Beta4QtlFor returns the inverse of the CDF (quantile) of the four-parameter Beta distribution, for a given probability.
func Beta4QtlFor(α, β, a, c, p float64) float64 {
	cdf := Beta4Qtl(α, β, a, c)
	return cdf(p)
}
