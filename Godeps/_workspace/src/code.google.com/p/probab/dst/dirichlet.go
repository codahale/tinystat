// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Dirichlet distribution. It is the multivariate generalization of the beta distribution.
// Parameters: 
// αi > 0	 	concentration parameters
// Support: 
// θi ∈ [0, 1] and Σθi = 1

// DirichletPDF returns the PDF of the Dirichlet distribution. 
func DirichletPDF(α []float64) func(θ []float64) float64 {
	return func(θ []float64) float64 {
		k := len(α)
		if len(θ) != k {
			return 0
		}
		l := float64(1.0)
		totalα := float64(0)
		for i := 0; i < k; i++ {
			if θ[i] < 0 || θ[i] > 1 {
				return 0
			}
			l *= pow(θ[i], α[i]-1)
			l /= Γ(α[i])
			totalα += α[i]
		}
		l *= Γ(totalα)
		return l
	}
}

// DirichletLnPDF returns the natural logarithm of the PDF of the Dirichlet distribution. 
func DirichletLnPDF(α []float64) func(x []float64) float64 {
	return func(x []float64) float64 {
		k := len(α)
		if len(x) != k {
			return negInf
		}
		l := fZero
		totalα := float64(0)
		for i := 0; i < k; i++ {
			if x[i] < 0 || x[i] > 1 {
				return negInf
			}
			l += (α[i] - 1) * log(x[i])
			l -= LnΓ(α[i])
			totalα += α[i]
		}
		l += LnΓ(totalα)
		return l
	}
}

// DirichletPDFAt returns the value of PDF of Dirichlet distribution at x. 
func DirichletPDFAt(α, θ []float64) float64 {
	pdf := DirichletPDF(α)
	return pdf(θ)
}

// DirichletNext returns random number drawn from the Dirichlet distribution. 
func DirichletNext(α []float64) []float64 {
	k := len(α)
	x := make([]float64, k)
	sum := fZero
	for i := 0; i < len(α); i++ {
		x[i] = GammaNext(α[i], 1.0)
		sum += x[i]
	}
	for i := 0; i < len(α); i++ {
		x[i] /= sum
	}
	return x
}

// Dirichlet returns the random number generator with  Dirichlet distribution. 
func Dirichlet(α []float64) func() []float64 {
	return func() []float64 { return DirichletNext(α) }
}

// DirichletMean returns the mean of the Dirichlet distribution. 
func DirichletMean(α []float64) []float64 {
	k := len(α)
	x := make([]float64, k)
	sum := fZero
	for i := 0; i < k; i++ {
		sum += α[i]
	}

	for i := 0; i < k; i++ {
		x[i] = α[i] / sum
	}
	return x
}

// DirichletMode returns the mode of the Dirichlet distribution. 
func DirichletMode(α []float64) []float64 {
	k := len(α)
	x := make([]float64, k)
	sum := fZero
	for i := 0; i < k; i++ {
		if α[i] <= 1 { // REVISION and citation NEEDED!
			panic("mode not defined")
		}
		sum += α[i]
	}

	for i := 0; i < k; i++ {
		x[i] = (α[i] - 1) / (sum - float64(k))
	}
	return x
}

// DirichletVar returns the variance of the Dirichlet distribution. 
func DirichletVar(α []float64) []float64 {
	k := len(α)
	x := make([]float64, k)
	sum := fZero
	for i := 0; i < k; i++ {
		sum += α[i]
	}

	for i := 0; i < k; i++ {
		x[i] = (α[i] * (sum - α[i])) / (sum * sum * (sum + 1))

	}
	return x
}
