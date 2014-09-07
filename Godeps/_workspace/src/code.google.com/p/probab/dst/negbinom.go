// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Negative binomial distribution. 
// A discrete probability distribution of the number of successes in a sequence of Bernoulli trials before a specified (non-random) number of failures (denoted r) occur. For example, if one throws a die repeatedly until the third time “1” appears, then the probability distribution of the number of non-“1”s that had appeared will be negative binomial.
//
// Parameters: 
// r > 0	 	number of failures until the experiment is stopped (integer, but the definition can also be extended to reals)
// ρ ∈ (0,1)	probability of success in each trial
//
// Support: 
// k ∈ { 0, 1, 2, 3, … }		number of successes

func do_search(p, pr float64, y, n, incr int64, z *float64) int64 {
	if *z >= p {
		// search to the left
	L1:
		for {
			*z = NegBinomialCDFAt(pr, n, y-incr)
			if y == 0 || *z < p {
				break L1
			}
			y = imax(0, y-incr)
		}
	} else { // search to the right

	L2:
		for {
			y += incr
			*z = NegBinomialCDFAt(pr, n, y-incr)
			if *z >= p {
				break L2
			}
		}
	}
	return y
}

// NegBinomialPMF returns the PMF of the Negative binomial distribution. 
func NegBinomialPMF(ρ float64, r int64) func(k int64) float64 {
	return func(k int64) float64 {
		return BinomCoeff(k+r-1, k) * pow(1-ρ, float64(r)) * pow(ρ, float64(k))
	}
}

// NegBinomialLnPMF returns the natural logarithm of the PMF of the Negative binomial distribution. 
func NegBinomialLnPMF(ρ float64, r int64) func(i int64) float64 {
	return func(k int64) float64 {
		rr := float64(r)
		return logChoose(k+r-1, r-1) + log(ρ)*rr + log(1-ρ)*float64(k)
	}
}

// NegBinomialPMFAt returns the value of PMF of Negative binomial distribution at k. 
func NegBinomialPMFAt(ρ float64, r, k int64) float64 {
	pmf := NegBinomialPMF(ρ, r)
	return pmf(k)
}

// NegBinomialCDF returns the CDF of the Negative binomial distribution. 
func NegBinomialCDF(ρ float64, r int64) func(k int64) float64 {
	return func(k int64) float64 {
		Ip := BetaCDFAt(float64(k+1), float64(r), ρ)
		return 1 - Ip
	}
}

// NegBinomialCDFAt returns the value of CDF of the Negative binomial distribution at k. 
func NegBinomialCDFAt(ρ float64, r, k int64) float64 {
	cdf := NegBinomialCDF(ρ, r)
	return cdf(k)
}

// NegBinomialQtl returns the inverse of the CDF (qquantile) of the Negative binomial distribution.
func NegBinomialQtl(ρ float64, r int64) func(p float64) int64 {
	return func(p float64) int64 {
		var pp, qq, mu, sigma, gamma, z float64
		var y int64
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
		if p+1.01*eps64 >= 1 {
			return int64(NaN)
		}

		// y := approx.value (Cornish-Fisher expansion)

		z = NormalQtlFor(0, 1, p)
		y = int64(floor(mu + sigma*(z+gamma*(z*z-1)/6) + 0.5))
		z = NegBinomialCDFAt(ρ, r, y)

		// fuzz to ensure left continuity
		p *= 1 - 64*eps64

		// If the C-F value is not too large a simple search is OK
		if y < 1e5 {
			return do_search(p, ρ, y, r, 1, &z)
		}
		// Otherwise be a bit cleverer in the search
		{
			incr := int64(floor(float64(y) / 1000))
			oldincr := incr
			for oldincr > 1 && incr > int64(floor(float64(y)*1e-15)) {
				//	    y = do_search(y, &z, p, r, ρ, incr)
				y = do_search(p, ρ, y, r, incr, &z)
				incr = imax(1, incr/100)
				oldincr = incr
			}
			return y
		}
	}
}

// NegBinomialQtlFor returns the inverse of the CDF (quantile) of the Negative binomial distribution, for given probability.
func NegBinomialQtlFor(ρ float64, r int64, p float64) int64 {
	qtl := NegBinomialQtl(ρ, r)
	return qtl(p)
}

// NegBinomialNext returns random number drawn from the Negative binomial distribution. 
func NegBinomialNext(ρ float64, r int64) int64 {
	k := iZero
	for r >= 0 {
		i := BernoulliNext(ρ)
		r -= i
		k += (1 - i)
	}
	return k
}

// NegBinomial returns the random number generator with  Negative binomial distribution. 
func NegBinomial(ρ float64, r int64) func() int64 {
	return func() int64 {
		return NegBinomialNext(ρ, r)
	}
}

// NegBinomialMean returns the mean of the Negative binomial distribution. 
func NegBinomialMean(ρ float64, r int64) float64 {
	return ρ * float64(r) / (1 - ρ)
}

// NegBinomialMode returns the mode of the Negative binomial distribution. 
func NegBinomialMode(ρ float64, r int64) float64 {
	if r > 1 {
		return floor(ρ * float64(r-1) / (1 - ρ))
	}
	return 0
}

// NegBinomialVar returns the variance of the Negative binomial distribution. 
func NegBinomialVar(ρ float64, r int64) float64 {
	return ρ * float64(r) / ((1 - ρ) * (1 - ρ))
}

// NegBinomialStd returns the standard deviation of the Negative binomial distribution. 
func NegBinomialStd(ρ float64, r int64) float64 {
	return sqrt(ρ*float64(r)) / (1 - ρ)
}

// NegBinomialSkew returns the skewness of the Negative binomial distribution. 
func NegBinomialSkew(ρ float64, r int64) float64 {
	return (1 + ρ) / sqrt(ρ*float64(r))
}

// NegBinomialExKurt returns the excess kurtosis of the Negative binomial distribution. 
func NegBinomialExKurt(ρ float64, r int64) float64 {
	rr := float64(r)
	return 6/rr + ((1-ρ)*(1-ρ))/(ρ*rr)
}

// NegBinomialMGF returns the moment-generating function of the Negative binomial distribution. 
func NegBinomialMGF(ρ float64, r int64, t float64) float64 {
	return pow((1-ρ)/(1-ρ*exp(t)), float64(r))
}

// NegBinomialPGF returns the probability-generating function of the Negative binomial distribution. 
func NegBinomialPGF(ρ float64, r int64, z float64) float64 {
	return pow((1-ρ)/(1-ρ*z), float64(r))
}
