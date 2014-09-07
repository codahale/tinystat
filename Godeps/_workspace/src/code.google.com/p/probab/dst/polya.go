// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Pólya distribution. 
// Extension of the negative binomial distribution to the case of a positive real parameter r. 
//
// Parameters: 
// r > 0	 	number of failures until the experiment is stopped (integer, but the definition can also be extended to reals)
// ρ ∈ (0,1)	probability of success in each trial
//
// Support: 
// k ∈ { 0, 1, 2, 3, … }		number of successes

func f_search(p, pr, y, n, incr float64, z *float64) float64 {
	if *z >= p {
		/* search to the left */
	L1:
		for {
			*z = PolyaCDFAt(pr, n, int64(floor(y-incr)))
			if y == 0 || *z < p {
				break L1
			}
			y = max(0, y-incr)
		}
	} else { /* search to the right */

	L2:
		for {
			y += incr
			*z = PolyaCDFAt(pr, n, int64(floor(y-incr)))
			if *z >= p {
				break L2
			}
		}
	}
	return y
}

// PolyaPMF returns the PMF of the Pólya distribution. 
func PolyaPMF(ρ, r float64) func(k int64) float64 {
	return func(k int64) float64 {
		kk := float64(k)
		return (Γ(kk+r) / (float64(fact(k)) * Γ(r))) * pow(1-ρ, r) * pow(ρ, float64(k))
	}
}

// PolyaPMFAt returns the value of PMF of Pólya distribution at k. 
func PolyaPMFAt(ρ, r float64, k int64) float64 {
	pmf := PolyaPMF(ρ, r)
	return pmf(k)
}

// PolyaCDF returns the CDF of the Pólya distribution. 
func PolyaCDF(ρ, r float64) func(k int64) float64 {
	return func(k int64) float64 {
		Ip := BetaCDFAt(float64(k+1), r, ρ)
		return 1 - Ip
	}
}

// PolyaCDFAt returns the value of CDF of the Pólya distribution, at k. 
func PolyaCDFAt(ρ, r float64, k int64) float64 {
	cdf := PolyaCDF(ρ, r)
	return cdf(k)
}

// PolyaMean returns the mean of the Pólya distribution. 
func PolyaMean(ρ, r float64) float64 {
	return ρ * r / (1 - ρ)
}

// PolyaMode returns the mode of the Pólya distribution. 
func PolyaMode(ρ, r float64) float64 {
	if r > 1 {
		return floor(ρ * (r - 1) / (1 - ρ))
	}
	return 0
}

// PolyaVar returns the variance of the Pólya distribution. 
func PolyaVar(ρ, r float64) float64 {
	return ρ * r / ((1 - ρ) * (1 - ρ))
}

// PolyaStd returns the standard deviation of the Pólya distribution. 
func PolyaStd(ρ, r float64) float64 {
	return sqrt(ρ*r) / (1 - ρ)
}

// PolyaSkew returns the skewness of the Pólya distribution. 
func PolyaSkew(ρ, r float64) float64 {
	return (1 + ρ) / sqrt(ρ*r)
}

// PolyaExKurt returns the excess kurtosis of the Pólya distribution. 
func PolyaExKurt(ρ, r float64) float64 {
	return 6/r + ((1-ρ)*(1-ρ))/(ρ*r)
}

// PolyaMGF returns the moment-generating function of the Pólya distribution. 
func PolyaMGF(ρ, r float64, t float64) float64 {
	return pow((1-ρ)/(1-ρ*exp(t)), r)
}

// PolyaPGF returns the probability-generating function of the Pólya distribution. 
func PolyaPGF(ρ, r float64, z float64) float64 {
	return pow((1-ρ)/(1-ρ*z), r)
}

// PolyaQtl returns the inverse of the CDF (quantile) of the Pólya distribution.
func PolyaQtl(ρ, r float64) func(p float64) int64 {
	return func(p float64) int64 {
		var eps, pp, qq, mu, sigma, gamma, z, y float64
		fr := float64(r)

		if ρ <= 0 || ρ > 1 || fr <= 0 { // FIXME: fr = 0 is well defined
			return int64(NaN)
		}

		if ρ == 1 {
			return 0
		}

		qq = 1.0 / ρ
		pp = (1.0 - ρ) * qq
		mu = fr * pp
		sigma = sqrt(fr * pp * qq)
		gamma = (qq + pp) / sigma

		// temporary hack --- FIXME ---
		if p+1.01*eps >= 1. {
			return int64(NaN)
		}

		// approximate by Cornish-Fisher expansion
		z = NormalQtlFor(0, 1, p)
		y = floor(mu + sigma*(z+gamma*(z*z-1)/6) + 0.5)
		z = PolyaCDFAt(ρ, r, int64(y))

		// fuzz to ensure left continuity
		p *= 1 - 64*eps

		// If the C-F value is not too large a simple search is OK
		if y < 1e5 {
			return int64(floor(f_search(p, ρ, y, r, 1, &z)))
		}

		// Otherwise be a bit cleverer in the search
		{
			incr := floor(y / 1000)
			oldincr := incr
			for oldincr > 1 && incr > floor(y*1e-15) {
				//	    y = do_search(y, &z, p, r, ρ, incr)
				y = f_search(p, ρ, y, r, incr, &z)
				incr = max(1, incr/100)
				oldincr = incr
			}
			return int64(floor(y))
		}
	}
}

// PolyaQtlFor returns the inverse of the CDF (quantile) of the Pólya distribution, for given probability.
func PolyaQtlFor(ρ, r float64, p float64) int64 {
	qtl := PolyaQtl(ρ, r)
	return qtl(p)
}
