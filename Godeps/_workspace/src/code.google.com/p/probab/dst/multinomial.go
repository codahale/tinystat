// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Multinomial distribution 
// is a generalization of the binomial distribution.
// The binomial distribution is the probability distribution of the number of "successes" in n independent Bernoulli trials, with the same probability of "success" on each trial. 
// In a multinomial distribution, the analog of the Bernoulli distribution is the categorical distribution, 
// where each trial results in exactly one of some fixed finite number k of possible outcomes, 
// with probabilities p1, ..., pk (so that pi ≥ 0 for i = 1, ..., k and \sum{i=1}^k pi = 1), and there are n independent trials. 
//
// Parameters: 
// n ∈ {1, 2 ... }	 	number of trials
// θ1, ..., θk  ∈ [0, 1]	event probabilities (must sum to1)
//
// Support: 
// xi ∈ {0, ... , n}
// Σxi = n

// MultinomialPMF returns the PMF of the Multinomial distribution. 
func MultinomialPMF(θ []float64, n int64) func(x []int64) float64 {
	return func(x []int64) float64 {
		if len(x) != len(θ) {
			return 0
		}
		l := fOne
		totalx := iZero
		for i := 0; i < len(x); i++ {
			l *= pow(θ[i], float64(x[i]))
			l /= Γ(float64(x[i] + 1))
			totalx += x[i]
		}
		if totalx != n {
			return 0
		}
		l *= Γ(float64(totalx + 1))
		return l
	}
}

// MultinomialLnPMF returns the natural logarithm of the PMF of the Multinomial distribution. 
func MultinomialLnPMF(θ []float64, n int64) func(x []int64) float64 {
	return func(x []int64) float64 {
		if len(x) != len(θ) {
			return negInf
		}
		l := fZero
		totalx := iZero
		for i := 0; i < len(x); i++ {
			l += log(θ[i]) * float64(x[i])
			l -= LnΓ(float64(x[i] + 1))
			totalx += x[i]
		}
		if totalx != n {
			return negInf
		}
		l += LnΓ(float64(totalx + 1))
		return l
	}
}

// MultinomialPMFAt returns the value of PMF of Multinomial distribution(μ, σ) at k. 
func MultinomialPMFAt(θ []float64, n int64, x []int64) float64 {
	pmf := MultinomialPMF(θ, n)
	return pmf(x)
}

// MultinomialNext returns random number drawn from the Multinomial distribution. 
func MultinomialNext(θ []float64, n int64) []int64 {
	x := make([]int64, len(θ))
	chooser := Choice(θ)
	for i := iZero; i < n; i++ {
		x[chooser()]++
	}
	return x
}

// Multinomial returns the random number generator with  Multinomial distribution. 
func Multinomial(θ []float64, n int64) func() []int64 {
	return func() []int64 {
		return MultinomialNext(θ, n)
	}
}

// MultinomialMean returns the mean of the Multinomial distribution. 
func MultinomialMean(θ []float64, n int64) []float64 {
	k := len(θ)
	x := make([]float64, k)
	for i := 0; i < k; i++ {
		x[i] = float64(n) * θ[i]
	}
	return x
}

// MultinomialVar returns the variance of the Multinomial distribution. 
func MultinomialVar(θ []float64, n int64) []float64 {
	k := len(θ)
	x := make([]float64, k)
	for i := 0; i < k; i++ {
		x[i] = float64(n) * θ[i] * (1 - θ[i])
	}
	return x
}

// MultinomialStd returns the standard deviation of the Multinomial distribution. 
func MultinomialStd(θ []float64, n int64) []float64 {
	k := len(θ)
	x := make([]float64, k)
	for i := 0; i < k; i++ {
		x[i] = sqrt(float64(n) * θ[i] * (1 - θ[i]))
	}
	return x
}

// MultinomialMGF returns the moment-generating function of the Multinomial distribution. 
func MultinomialMGF(θ []float64, n int64, t []float64) float64 {
	k := len(θ)
	sum := 0.0
	for i := 0; i < k; i++ {
		sum += θ[i] * exp(t[i])
	}
	return pow(sum, float64(n))
}

// MultinomialPGF returns the probability-generating function of the Multinomial distribution. 
func MultinomialPGF(θ []float64, n int64, z []float64) float64 {
	k := len(θ)
	sum := 0.0
	for i := 0; i < k; i++ {
		sum += θ[i] * z[i]
	}
	return pow(sum, float64(n))
}
