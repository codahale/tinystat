// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Gamma distribution. 
// Parameters: 
// α > 0.0		shape parameter, 
// θ (Theta) > 0.0	scale parameter. 
// Alternatively, shape parameter α = k and an inverse scale parameter β = 1⁄θ, is called a rate parameter.
// If α is an integer, then the distribution represents an Erlang distribution; i.e., the sum of k  independent exponentially-distributed random variables, each of which has a mean of θ (which is equivalent to a rate parameter of 1/θ). Equivalently, if α is an integer, then the distribution again represents an Erlang distribution, i.e. the sum of α independent exponentially-distributed random variables, each of which has a mean of 1/β (which is equivalent to a rate parameter of β).
// Support: 
// x ∈ (0, ∞)

// GammaPDF returns the value of CDF of the Gamma distribution, at x. 
func GammaPDF(α float64, θ float64) func(x float64) float64 {
	//  Computes the density of the gamma distribution,
	//
	//                   1/θ (x/θ)^{α-1} exp(-x/θ)
	//       p(x;α,θ) = -----------------------
	//                            (α-1)!
	//
	//  where `α' is the shape parameter, and
	// `θ' is the scale (= 1/λ in other parametrizations).
	return func(x float64) float64 {
		var pr float64
		if isNaN(x) || isNaN(α) || isNaN(θ) {
			return x + α + θ
		}
		if α < 0 || θ <= 0 {
			return NaN
		}
		if x < 0 {
			return negInf
		}
		if α == 0 {
			//	return (x == 0)? ML_POSINF : R_D__0;
			if x == 0 {
				return posInf
			} else {
				return 0
			}
		}
		if x == 0 {
			if α < 1 {
				return posInf
			}
			if α > 1 {
				return 0
			}
			// else
			return 1 / θ
		}
		if α < 1 {
			pr = dpois_raw(α, x/θ)
			return pr * α / x
		}
		// else  α >= 1
		pr = dpois_raw(α-1, x/θ)
		return pr / θ
	}
}

// GammaPDF2 returns the PDF of the Gamma distribution. Another tested implementation.
func GammaPDF2(α, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return 0
		}
		return pow(x, α-1) * exp(-x/θ) / (Γ(α) * pow(θ, α))
	}
}

// GammaPDFAt returns the value of PDF of Gamma distribution at x. 
func GammaPDFAt(k, θ, x float64) float64 {
	pdf := GammaPDF(k, θ)
	return pdf(x)
}

// GammaLnPDF returns the natural logarithm of the value of PDF of the Gamma distribution, at x. 
func GammaLnPDF(α float64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		var pr float64
		if isNaN(x) || isNaN(α) || isNaN(θ) {
			return x + α + θ
		}
		if α < 0 || θ <= 0 {
			return NaN
		}
		if x < 0 {
			return negInf
		}
		if α == 0 {
			if x == 0 {
				return posInf
			} else {
				return negInf
			}
		}
		if x == 0 {
			if α < 1 {
				return posInf
			}
			if α > 1 {
				return negInf
			}
			// else
			return -log(θ)
		}

		if α < 1 {
			pr = dpois_raw_ln(α, x/θ)
			return pr + log(α/x)
		}
		// else  α >= 1
		pr = dpois_raw_ln(α-1, x/θ)
		return pr - log(θ)
	}
}

// GammaLnPDF2 returns the PDF of Gamma distribution. Another tested implementation.
func GammaLnPDF2(α float64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if x < 0 {
			return negInf
		}
		//		return log(pow(x, α-1) * exp(-x/θ) / (Γ(α) * pow(θ, α)))
		//		return log(pow(x, α-1)) + log(exp(-x/θ)) - log(Γ(α)) - log(pow(θ, α))
		//		return log(pow(x, α-1)) -x/θ - log(Γ(α)) - log(pow(θ, α))
		//		return (α-1)*log(x) -x/θ - log(Γ(α)) - α* log(θ)
		return (α-1)*log(x) - x/θ - LnΓ(α) - α*log(θ)
	}
}

// GammaLnPDFAt returns the value of PDF of Gamma distribution at x. 
func GammaLnPDFAt(α, θ, x float64) float64 {
	pdf := GammaLnPDF(α, θ)
	return pdf(x)
}

// GammaCDF returns the CDF of the Gamma distribution.
func GammaCDF(α float64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(α) || isNaN(θ) {
			return x + α + θ
		}
		if α < 0. || θ <= 0. {
			return NaN
		}
		if x == 0 {
			return 0
		}
		x /= θ
		if isNaN(x) { // eg. original x = θ = +Inf 
			return x
		}
		if α == 0. { // limit case; useful e.g. in pnchisq()
			if x <= 0 {
				return 0
			} else {
				return 1
			}
		}
		return pgamma_raw(x, α)
	}
}

// GammaCDFAt returns the value of CDF of the Gamma distribution, at x. 
func GammaCDFAt(α, θ, x float64) float64 {
	cdf := GammaCDF(α, θ)
	return cdf(x)
}

// GammaLnCDF returns the value of CDF of the Gamma distribution, at x. 
func GammaLnCDF(α float64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isNaN(x) || isNaN(α) || isNaN(θ) {
			return x + α + θ
		}
		if α < 0. || θ <= 0. {
			return NaN
		}
		if x == 0 {
			return negInf
		}
		x /= θ
		if isNaN(x) { // eg. original x = θ = +Inf 
			return x
		}
		if α == 0 { // limit case; useful e.g. in pnchisq()
			if x <= 0 {
				return negInf
			} else {
				return 0
			}
		}
		return pgamma_raw_ln(x, α)
	}
}

// GammaLnCDFAt returns the value of CDF of the Gamma distribution, at x. 
func GammaLnCDFAt(α, θ, x float64) float64 {
	cdf := GammaLnCDF(α, θ)
	return cdf(x)
}

// GammaNext returns random number drawn from the Gamma distribution. 
func GammaNext(α float64, θ float64) float64 {
	//if α is a small integer, this way is faster on my laptop
	if α == float64(int64(α)) && α <= 15 {
		x := ExponentialNext(θ)
		for i := 1; i < int(α); i++ {
			x += ExponentialNext(θ)
		}
		return x
	}

	if α < 0.75 {
		return RejectionSample(GammaPDF(α, θ), ExponentialPDF(θ), Exponential(θ), 1)
	}

	//Tadikamalla ACM '73
	a := α - 1
	b := 0.5 + 0.5*sqrt(4*α-3)
	c := a * (1 + b) / b
	d := (b - 1) / (a * b)
	s := a / b
	p := 1.0 / (2 - exp(-s))
	var x, y float64
	for i := 1; ; i++ {
		u := UniformNext(0, 1)
		if u > p {
			var e float64
			for e = -log((1 - u) / (1 - p)); e > s; e = e - a/b {
			}
			x = a - b*e
			y = a - x
		} else {
			x = a - b*log(u/p)
			y = x - a
		}
		u2 := UniformNext(0, 1)
		if log(u2) <= a*log(d*x)-x+y/b+c {
			break
		}
	}
	return x / θ
}

// Gamma returns the random number generator with  Gamma distribution. 
func Gamma(α, θ float64) func() float64 {
	return func() float64 { return GammaNext(α, θ) }
}

// GammaMean returns the mean of the Gamma distribution. 
func GammaMean(α, θ float64) float64 {
	return α * θ
}

// GammaMode returns the mode of the Gamma distribution. 
func GammaMode(α, θ float64) float64 {
	if α <= 1 {
		return NaN
	}
	return (α - 1) * θ
}

// GammaVar returns the variance of the Gamma distribution. 
func GammaVar(α, θ float64) float64 {
	return α * θ * θ
}

// GammaStd returns the standard deviation of the Gamma distribution. 
func GammaStd(α, θ float64) float64 {
	return sqrt(α) * θ
}

// GammaSkew returns the skewness of the Gamma distribution. 
func GammaSkew(α, θ float64) float64 {
	return 2 / sqrt(α)
}

// GammaRateToScale returns the parameter θ (scale) of the Gamma distribution calculated from β = rate.
// α = shape, β = rate
// To be used to reparametrize the Gamma distribution. 
func GammaRateToScale(β float64) (θ float64) {
	θ = 1 / β
	return
}

// GammaReparamModeStd returns the parameters α, θ (shape, scale) of the Gamma distribution calculated from mode and standard deviation. 
// It is more intuitive to start with the mode and standard deviation, instead of the mean and standard deviation as used in the Kruschke (2011) book. 
// The reason is that the gamma distribution is typically very skewed, and therefore the location of the mean is not very intuitive. 
// This function computes the shape and rate parameters of the gamma distribution from a desired mode and standard deviation.
// After http://doingbayesiandataanalysis.blogspot.com/2012/01/parameterizing-gamma-distribution-by.html
func GammaReparamModeStd(mode, sd float64) (α, θ float64) {
	β := (mode + sqrt(mode*mode+4*sd*sd)) / (2 * sd * sd)
	α = 1 + mode*β
	θ = 1 / β
	return
}

// GammaReparamMeanStd returns the parameters α, θ (shape, scale) of the Gamma distribution calculated from mean and standard deviation. 
func GammaReparamMeanStd(mean, sd float64) (α, θ float64) {
	/*
		mean =α*θ
		sd*sd = α*θ*θ
		α=sd*sd /(θ*θ)
		α=mean/θ
		sd*sd /(θ*θ)=mean/θ
		sd*sd *θ/(θ*θ)=mean
		θ/(θ*θ)=mean/(sd*sd )
		θ=(sd*sd )/mean
	*/
	θ = (sd * sd) / mean
	α = mean / θ
	return
}

/************** some non-working code

// GammaCDF returns the CDF of the Gamma distribution. // TO BE REIMPLEMENTED
// Analytic solution, did not pass some tests!
func GammaCDF(k float64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if k < 0 || θ < 0 {
			return NaN
		}
		if x < 0 {
			return 0
		}
		return Iγ(k, x/θ) / Γ(k)
	}
}

// GammaCDFint returns the CDF of the Gamma distribution, for integer k only. 
// Cumulative distribution function, for integer k only
func GammaCDFint(k int64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if k < 0 || θ < 0 {
			return NaN
		}
		if x < 0 {
			return 0
		}
		return Iγint(k, x/θ) / Γ(float64(k))
	}
}

// Cumulative distribution function, using gamma incomplete integral  DOES NOT WORK !!!
func GammaCDF(k float64, θ float64) func(x float64) float64 {
	return func(x float64) float64 {
		if k < 0 || θ < 0 {
			return NaN
		}
		if x < 0 {
			return 0
		}
		return IGam(θ, k*x)
	}
}


// GammaLnPDF returns the natural logarithm of the PDF of the Gamma distribution. 
func GammaLnPDF(α float64, θ float64) func(x float64) float64 {
	expPart := ExponentialLnPDF(θ)
	return func(x float64) float64 {
		if x < 0 {
			return negInf
		}
		return expPart(x) + (α-1)*log(θ*x) - LnΓ(α)
	}
}
// GammaQtl returns the inverse of the CDF (quantile) of the Gamma distribution. 
func GammaQtl(α, θ float64) func(p float64) float64 {
	return func(p float64) float64 {
		var eps, ynew, h float64
		if p == 0 {
			return 0
		}
		if p == 1 {
			return posInf
		}

		eps = 1e-10
		y := α * θ
		yold := y
	L:
		for i := 0; i < 100; i++ {
			h = (GammaCDFAt(α, θ, yold) - p) / GammaPDFAt(α, θ, yold)
			ynew = yold - h
			if ynew <= eps {
				ynew = yold / 10
				h = yold - ynew
			}
			if abs(h) < eps {
				break L
			}
			yold = ynew
		}
		return ynew
	}
}

*/
