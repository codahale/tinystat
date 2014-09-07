// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Inverse-gamma distribution (not to be confused with Inverse CDF of Gamma distribution). 
//
// A two-parameter family of continuous probability distributions on the positive real line, 
// which is the distribution of the reciprocal of a variable distributed according to the gamma distribution. 
// Perhaps the chief use of the inverse gamma distribution is in Bayesian statistics, 
// where it serves as the conjugate prior of the variance of a normal distribution. 
// However, it is common among Bayesians to consider an alternative parametrization 
// of the normal distribution in terms of the precision, defined as the reciprocal of the variance, 
// which allows the gamma distribution to be used directly as a conjugate prior.
//
// Parameters: 
// α > 0:		shape
// β > 0:		scale
// Support:	x ∈ (0, ∞)

// InvGammaPDF returns the PDF of the InvGamma distribution. 
func InvGammaPDF(α, β float64) func(x float64) float64 {
	return func(x float64) float64 {

		//  We work with the density expressed as
		//  u^α * e^(-u) / (x * gamma(α))
		// with u = β/x.

		if isInf(α, 0) || isInf(β, 0) || α <= 0 || β <= 0 {
			return NaN
		}

		if isInf(x, 0) || x <= 0 {
			return 0
		}

		logu := log(β) - log(x)
		return exp(α*logu - exp(logu) - log(x) - LnΓ(α))
	}
}

// InvGammaLnPDF returns the natural logarithm of the PDF of the InvGamma distribution. 
func InvGammaLnPDF(α, β float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isInf(α, 0) || isInf(β, 0) || α <= 0 || β <= 0 {
			return NaN
		}

		if isInf(x, 0) || x <= 0 {
			return negInf
		}

		logu := log(β) - log(x)
		return α*logu - exp(logu) - log(x) - LnΓ(α)
	}
}

// InvGammaPDFAt returns the value of PDF of InvGamma distribution at x. 
func InvGammaPDFAt(α, β, x float64) float64 {
	pdf := InvGammaPDF(α, β)
	return pdf(x)
}

// InvGammaCDF returns the CDF of the InvGamma distribution. 
func InvGammaCDF(α, β float64) func(x float64) float64 {
	return func(x float64) float64 {
		if isInf(α, 0) || isInf(β, 0) || α <= 0 || β <= 0 {
			return NaN
		}

		if x <= 0 {
			return 0
		}
		return iΓ(α, β/x) / Γ(α)
	}
}

// InvGammaCDFAt returns the value of CDF of the InvGamma distribution, at x. 
func InvGammaCDFAt(α, β, x float64) float64 {
	cdf := InvGammaCDF(α, β)
	return cdf(x)
}

// InvGammaQtl returns the inverse of the CDF (quantile) of the InvGamma distribution. 
func InvGammaQtl(α, β float64) func(p float64) float64 {
	return func(p float64) float64 {
		if isInf(α, 0) || isInf(β, 0) || α <= 0 || β <= 0 {
			return NaN
		}
		if p < 0 || p > 1 {
			return NaN
		}
		if p == 0 {
			return 0
		}
		if p == 1 {
			return posInf
		}
		//R_D_qIv
		//    return β / qgamma(p, α, 1.0, !lower_tail, 0);
		return β / GammaQtlFor(α, 1.0, 1-p)
	}
}

// InvGammaQtlFor returns the inverse of the CDF (quantile) of the InvGamma distribution, for given probability.
func InvGammaQtlFor(δ, γ, p float64) float64 {
	qtl := InvGammaQtl(δ, γ)
	return qtl(p)
}

// InvGammaMean returns the mean of the InvGamma distribution. 
func InvGammaMean(α, β float64) float64 {
	if α <= 1 {
		return NaN
	}
	return β / (α - 1)
}

// InvGammaMedian returns the median of the InvGamma distribution. 
// to be implemented ...

// InvGammaMode returns the mode of the InvGamma distribution. 
func InvGammaMode(α, β float64) float64 {
	return β / (α + 1)
}

// InvGammaVar returns the variance of the InvGamma distribution. 
func InvGammaVar(α, β float64) float64 {
	if α <= 2 {
		return NaN
	}
	return (β * β) / ((α - 1) * (α - 1) * (α - 2))
}

// InvGammaStd returns the standard deviation of the InvGamma distribution. 
func InvGammaStd(α, β float64) float64 {
	if α <= 2 {
		return NaN
	}
	return β / sqrt((α-1)*(α-1)*(α-2))
}

// InvGammaSkew returns the skewness of the InvGamma distribution. 
func InvGammaSkew(α, β float64) (s float64) {
	if α <= 3 {
		return NaN
	}
	return 4 * sqrt(α-2) / (α - 3)
}

// InvGammaExKurt returns the excess kurtosis of the InvGamma distribution. 
func InvGammaExKurt(α, β float64) float64 {
	if α <= 4 {
		return NaN
	}
	return (30*α - 66) / ((α - 3) * (α - 4))
}

// InvGammaMGF returns the moment-generating function of the InvGamma distribution. To be implemented ...

/*  some old code...

// InvGammaPDF returns the PDF of the InvGamma distribution. 
func InvGammaPDF(α, β float64) func(x float64) float64 {
	return func(x float64) float64 {
		return exp(α*log(β) - LnΓ(α) - (α+1)*log(x) - β/x)
	}
}

// InvGammaLnPDF returns the natural logarithm of the PDF of the InvGamma distribution. 
func InvGammaLnPDF(α, β float64) func(x float64) float64 {
	return func(x float64) float64 {
		return α*log(β) - LnΓ(α) - (α+1)*log(x) - β*1.0/x
	}
}

// InvGammaCDF returns the CDF of the InvGamma distribution. 
func InvGammaCDF(α, β float64) func(x float64) float64 {
	return func(x float64) float64 {

		if isInf(α, 0) || isInf(β, 0) || α <= 0 || β <= 0 {
			return NaN
		}

		if x <= 0 {
			return 0
		}

		u := exp(log(β) - log(x))

		//    return pgamma(u, α, 1.0, !lower_tail, log_p);
		return 1 - GammaCDFAt(α, 1.0, u)
	}
}
*/
