// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Beta distribution. 
// Parameters:
// α > 0:		shape
// β > 0:		shape
// Support:	x ∈ [0; 1]

import (
	"fmt"
)

func bisect(x, p, a, b, xtol, ptol float64) float64 {
	var x0, x1, px float64
	cdf := BetaPDF(a, b)
	for abs(x1-x0) > xtol {
		px = cdf(x)
		switch {
		case abs(px-p) < ptol:
			return x
		case px < p:
			x0 = x
		case px > p:
			x1 = x
		}
		x = 0.5 * (x0 + x1)
	}
	return x
}

func betaContinuedFraction(α, β, x float64) float64 {
	var aa, del, res, qab, qap, qam, c, d, m2, m, acc float64
	var i int64
	const eps = 2.2204460492503131e-16
	const maxIter = 1000000000

	acc = 1e-16
	qab = α + β
	qap = α + 1.0
	qam = α - 1.0
	c = 1.0
	d = 1.0 - qab*x/qap

	if abs(d) < eps {
		d = eps
	}
	d = 1.0 / d
	res = d

	for i = 1; i <= maxIter; i++ {
		m = (float64)(i)
		m2 = 2 * m
		aa = m * (β - m) * x / ((qam + m2) * (α + m2))
		d = 1.0 + aa*d
		if abs(d) < eps {
			d = eps
		}
		c = 1.0 + aa/c
		if abs(c) < eps {
			c = eps
		}
		d = 1.0 / d
		res *= d * c
		aa = -(α + m) * (qab + m) * x / ((α + m2) * (qap + m2))
		d = 1.0 + aa*d
		if abs(d) < eps {
			d = eps
		}
		c = 1.0 + aa/c
		if abs(c) < eps {
			c = eps
		}
		d = 1.0 / d
		del = d * c
		res *= del
		if abs(del-1.0) < acc {
			return res
		}
	}

	panic(fmt.Sprintf("betaContinuedFraction(): α or β too big, or maxIter too small"))
	return -1.00
}

// BetaPDF returns the PDF of the Beta distribution. 
func BetaPDF(α, β float64) func(x float64) float64 {
	if α == 1 && β == 1 { // uniform case
		return UniformPDF(0, 1)
	}
	dα := []float64{α, β}
	dirPDF := DirichletPDF(dα)
	return func(x float64) float64 {
		if 0 > x || x > 1 {
			return 0
		}
		dx := []float64{x, 1 - x}
		return dirPDF(dx)
	}
}

// BetaLnPDF returns the natural logarithm of the PDF of the Beta distribution. 
func BetaLnPDF(α, β float64) func(x float64) float64 {
	dα := []float64{α, β}
	dirLnPDF := DirichletLnPDF(dα)
	return func(x float64) float64 {
		if 0 > x || x > 1 {
			return negInf
		}
		dx := []float64{x, 1 - x}
		return dirLnPDF(dx)
	}
}

// BetaPDFAt returns the value of PDF of Beta distribution at x. 
func BetaPDFAt(α, β, x float64) float64 {
	pdf := BetaPDF(α, β)
	return pdf(x)
}

// BetaCDF returns the CDF of the Beta distribution. 
func BetaCDF(α, β float64) func(x float64) float64 {
	if α == 1 && β == 1 { // uniform case
		return UniformCDF(0, 1)
	}
	return func(x float64) float64 {
		var y, res float64
		y = exp(LnΓ(α+β) - LnΓ(α) - LnΓ(β) + α*log(x) + β*log(1.0-x))
		switch {
		case x == 0:
			res = 0.0
		case x == 1.0:
			res = 1.0
		case x < (α+1.0)/(α+β+2.0):
			res = y * betaContinuedFraction(α, β, x) / α
		default:
			res = 1.0 - y*betaContinuedFraction(β, α, 1.0-x)/β

		}
		return res
	}
}

// BetaCDFAt returns the value of CDF of the Beta distribution, at x. 
func BetaCDFAt(α, β, x float64) float64 {
	cdf := BetaCDF(α, β)
	return cdf(x)
}

// BetaQtl returns the inverse of the CDF (quantile) of the Beta distribution. 
func BetaQtl(α, β float64) func(p float64) float64 {
	// p: probability for which the quantile is evaluated
	return func(p float64) float64 {
		var x float64 = 0
		var a float64 = 0
		var b float64 = 1
		var precision float64 = 1e-9
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
		return x
	}
}

// BetaQtlFor returns the inverse of the CDF (quantile) of the Beta distribution, for given probability.
func BetaQtlFor(α, β, p float64) float64 {
	cdf := BetaQtl(α, β)
	return cdf(p)
}

// BetaNext returns random number drawn from the Beta distribution. 
func BetaNext(α, β float64) float64 {
	if α == 1 && β == 1 { // uniform case
		return UniformNext(0, 1)
	}
	dα := []float64{α, β}
	return DirichletNext(dα)[0]
}

// Beta returns the random number generator with  Beta distribution. 
func Beta(α, β float64) func() float64 {
	if α == 1 && β == 1 { // uniform case
		return Uniform(0, 1)
	}
	return func() float64 { return BetaNext(α, β) }
}

// BetaMean returns the mean of the Beta distribution. 
func BetaMean(α, β float64) (μ float64) {
	if α == β { // symmetric case
		μ = 0.5
	} else {
		μ = α / (α + β)
	}
	return
}

// BetaMedian returns the median of the Beta distribution. 
func BetaMedian(α, β float64) (med float64) {
	//The median of the beta distribution is the unique real number 
	// for which the regularized incomplete beta function  = 0.5 . 
	// There is no general closed-form expression for the median of the beta distribution 
	// for arbitrary values of α and β. 
	switch {
	case α == β: // symmetric case
		med = 0.5
	case α == 1 && β > 0:
		med = 1.0 - pow(0.5, 1/β)
	case β == 1 && α > 0:
		med = pow(0.5, 1/α)
	case α == 3 && β == 2:
		med = 0.6142724318676105
	case α == 2 && β == 3:
		med = 0.38572756813238945
	case α <= 1 || β <= 1:
		med = (α - 1/3) / (α + β - 2/3) // approximation
	default:
		panic("no closed form for median, sorry")
	}
	return
}

// BetaMedianApprox returns the approximate median of the Beta distribution. 
func BetaMedianApprox(α, β float64) float64 {
	if α <= 1 || β <= 1 {
		return NaN
	}
	return (α - 1/3) / (α + β - 2/3)
}

// BetaMode returns the mode of the Beta distribution. 
func BetaMode(α, β float64) float64 {
	if α <= 1 || β <= 1 {
		return NaN
	}
	return (α - 1) / (α + β - 2) // if α < 1 and β < 1, this is the anti-mode
}

// BetaVar returns the variance of the Beta distribution. 
func BetaVar(α, β float64) float64 {
	return (α * β) / ((α + β) * (α + β) * (α + β + 1))
}

// BetaStd returns the standard deviation of the Beta distribution. 
func BetaStd(α, β float64) float64 {
	v := (α * β) / ((α + β) * (α + β) * (α + β + 1))
	return sqrt(v)
}

// BetaSkew returns the skewness of the Beta distribution. 
func BetaSkew(α, β float64) (s float64) {

	if α == β { // symmetric case
		s = 0.0
	} else {
		num := 2 * (β - α) * sqrt(α+β+1)
		den := (α + β + 2) * sqrt(α*β)
		s = num / den
	}
	return
}

// BetaExKurt returns the excess kurtosis of the Beta distribution. 
func BetaExKurt(α, β float64) float64 {
	num := 6 * ((α-β)*(α-β)*(α+β+1) - α*β*(α+β+2))
	den := α * β * (α + β + 2) * (α + β + 2)
	return num / den
}

// BetaReparamMeanStd returns the parameters α, β of the Beta distribution calculated from desired mean and standard deviation. 
// To be used to reparametrize the Beta distribution. 
func BetaReparamMeanStd(μ, σ float64) (α, β float64) {
	// http://linkage.rockefeller.edu/pawe3d/help/Beta-distribution.html
	if σ*σ >= μ*(1-μ) {
		return NaN, NaN
	}
	α = (μ*μ - μ*μ*μ - μ*σ*σ) / (σ * σ)
	β = (μ - 2*μ*μ + μ*μ*μ - σ*σ + μ*σ*σ) / (σ * σ)
	return
}

// BetaReparamModStd returns the parameters α, β of the Beta distribution calculated from modus and standard deviation. 
// To be used to reparametrize the Beta distribution. 
// To be implemented
/* func BetaReparamModStd(μ, σ float64) (α, β float64) {
	if σ*σ >= μ*(1-μ) {
			return NaN
	}
	α = 
	β = 
	return
}

μ=(α - 1) / (α + β - 2)
μ*(α+β-2)-(α-1)=0
μ*α+μ*β-2*μ-α+1=0
μ*α-α=-μ*β+2*μ-1
α*(μ-1)=-μ*β+2*μ-1
α=(-μ*β+2*μ-1)  /(μ-1)

(σ*σ)=(α * β) /((α + β) * (α + β) * (α + β + 1))
(σ*σ)*((α + β) * (α + β) * (α + β + 1))  - (α * β) =0
	(α + β) * (α + β) =(α*α+2*α*β+β*β)
		(α*α+2*α*β+β*β) * (α + β + 1) = (α*α*α+2*α*α*β+α*β*β+α*α*β+2*α*β*β+β*β*β+ α*α+2*α*β+β*β) =
			(σ*σ)*... = (α*α*α*σ*σ+2*α*α*β*σ*σ+α*β*β*σ*σ+α*α*β*σ*σ+2*α*β*β*σ*σ+β*β*β*σ*σ+ α*α*σ*σ+2*α*β*σ*σ+β*β*σ*σ)
(α*α*α*σ*σ+2*α*α*β*σ*σ+α*β*β*σ*σ+α*α*β*σ*σ+2*α*β*β*σ*σ+β*β*β*σ*σ+ α*α*σ*σ+2*α*β*σ*σ+β*β*σ*σ -α*β) 
+2*α*α*β*σ*σ  +α*β*β*σ*σ  +α*α*β*σ*σ  +2*α*β*β*σ*σ +β*β*β*σ*σ +2*α*β*σ*σ  +β*β*σ*σ -α*β)+α*α*α*σ*σ+ α*α*σ*σ=0

β*β*β*σ*σ +  β*β*(α*σ*σ+2*α*σ*σ+σ*σ) +   β*(+2*α*α*σ*σ +α*α*σ*σ+2*α*σ*σ -α)+α*α*α*σ*σ+ α*α*σ*σ=0

a=σ*σ
b=(α*σ*σ+2*α*σ*σ+σ*σ)
c=(2*α*α*σ*σ +α*α*σ*σ+2*α*σ*σ -α)
d=(α*α*α*σ*σ+ α*α*σ*σ)
Δ=18*a*b*c*d - 4*b*b*b*d +b*b*c*c -4*a*c*c*c - 27*a*a*d*d
================
mod=(α - 1) / (α + β - 2)
mod*(α + β - 2)=(α - 1)

μ = α / (α + β)
	α := μ*(μ*(1-μ)/(σ*σ)-1)
	β := (1-μ)*(μ*(1-μ)/(σ*σ)-1)

mod=(α - 1)/(α + β - 2)
mod=((μ*(μ*(1-μ)/(σ*σ)-1)) - 1)/((μ*(μ*(1-μ)/(σ*σ)-1)) + (1-μ)*(μ*(1-μ)/(σ*σ)-1) - 2)
mod=((A) - 1)/((μ*(μ*(1-μ)/(σ*σ)-1)) + (1-μ)*(μ*(1-μ)/(σ*σ)-1) - 2)
A = μ*(μ*(1-μ)/(σ*σ)-1) = (μ*μ*(1-μ)/(σ*σ)-μ) = (μ*μ-μ*μ*μ)/(σ*σ-μ) 
B = 


...
*/

// Beta4Transform transforms Beta Distribution with the support [0,1]  to a Beta Distribution with the support [a, b].
func Beta4Transform(a, b, x float64) float64 {
	return (b-a)*x + a
}
