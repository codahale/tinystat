// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Pareto Type II distribution. 
// blah blah. 
//
// Parameters: 
// θ > 0.0 (scale)
// α > 0.0 (shape)
//
// Support: 
// x >= θ
// inspired by R:actuar

import "math/rand"

// ParetoIIChkParams checks parameters of the ParetoII distribution. 
func ParetoIIChkParams(θ, α float64) bool {
	ok := true
	if α <= 0 || θ <= 0 {
		ok = false
	}
	return ok
}

// ParetoIIChkSupport checks support of the ParetoII distribution. 
func ParetoIIChkSupport(x float64) bool {
	ok := true
	if x < 0 {
		ok = false
	}
	return ok
}

// ParetoIIPDF returns the PDF of the ParetoII distribution. 
func ParetoIIPDF(θ, α float64) func(x float64) float64 {
	// We work with the density expressed as
	// α * u^α * (1 - u) / x
	// with u = 1/(1 + v), v = x/θ.
	return func(x float64) float64 {
		var p float64
		if x < 0 {
			p = 0
		} else if x == 0 {
			p = α / θ
		} else {
			tmp := log(x) - log(θ)
			logu := -log1p(exp(tmp))
			log1mu := -log1p(exp(-tmp))
			p = exp(log(α) + α*logu + log1mu - log(x))
		}
		return p
	}
}

// ParetoIIPDFAt returns the value of PDF of Pareto Type II distribution at x. 
func ParetoIIPDFAt(θ, α, x float64) float64 {
	pdf := ParetoIIPDF(θ, α)
	return pdf(x)
}

// ParetoIICDF returns the CDF of the Pareto Type II distribution. 
func ParetoIICDF(θ, α float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		u := exp(-log1p(exp(log(x) - log(θ))))
		return 1 - pow(u, α)
	}
}

// ParetoIICDFAt returns the value of CDF of the Pareto Type II distribution, at x. 
func ParetoIICDFAt(θ, α, q, x float64) float64 {
	cdf := ParetoIICDF(θ, α)
	return cdf(x)
}

// ParetoIIQtl returns the inverse of the CDF (quantile) of the Pareto Type II distribution. 
func ParetoIIQtl(θ, α float64) func(p float64) float64 {
	return func(p float64) float64 {
		if p < 0 || p > 1 {
			return NaN
		}
		return θ * (pow((0.5-(p)+0.5), -1.0/α) - 1.0)
	}
}

// ParetoIIQtlFor returns the inverse of the CDF (quantile) of the Pareto Type II distribution, for given probability.
func ParetoIIQtlFor(θ, α, p float64) float64 {
	cdf := ParetoQtl(θ, α)
	return cdf(p)
}

// ParetoIINext returns random number drawn from the Pareto Type II distribution. 
func ParetoIINext(θ, α float64) float64 {
	qtl := ParetoIIQtl(θ, α)
	p := rand.Float64()
	return qtl(p)
}

// ParetoII returns the random number generator with  Pareto Type II distribution. 
func ParetoII(θ, α float64) func() float64 {
	return func() float64 { return ParetoIINext(θ, α) }
}

// ParetoIIMoment returns the n-th moment of the Pareto Type II distribution. 
func ParetoIIMoment(θ, α float64, order int) float64 {
	o := float64(order)
	if α <= o {
		return NaN
	}
	return pow(θ, o) * Γ(1.0+o) * Γ(α-o) / Γ(α)
}

// ParetoIIMean returns the mean of the Pareto Type II distribution. 
func ParetoIIMean(θ, α float64) float64 {
	return ParetoIIMoment(θ, α, 1)
}

// ParetoIIVar returns the variance of the Pareto Type II distribution. 
func ParetoIIVar(θ, α float64) float64 {
	return ParetoIIMoment(θ, α, 2)
}

// ParetoIISkew returns the skewness of the Pareto Type II distribution. 
func ParetoIISkew(θ, α float64) float64 {
	return ParetoIIMoment(θ, α, 3)
}

// ParetoIIExKurt returns the excess kurtosis of the Pareto Type II distribution. 
func ParetoIIExKurt(θ, α float64) float64 {
	return ParetoIIMoment(θ, α, 4)
}
