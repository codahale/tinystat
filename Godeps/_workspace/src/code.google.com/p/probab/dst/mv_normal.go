// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Multivariate normal distribution. 
// The multivariate normal distribution or multivariate Gaussian distribution, is a generalization of the one-dimensional (univariate) normal distribution to higher dimensions. One possible definition is that a random vector is said to be k-variate normally distributed if every linear combination of its k components has a univariate normal distribution. However, its importance derives mainly from the multivariate central limit theorem. The multivariate normal distribution is often used to describe, at least approximately, any set of (possibly) correlated real-valued random variables each of which clusters around a mean value.
//
// Parameters: 
// μ ∈ ℝk		(Rk)	location
// Σ ∈ ℝk✕k	(Rkxk)	covariance (nonnegative-definite matrix)
//
// Support: 
// x ∈ μ+span(Σ) ⊆ ℝk

import (
	. "github.com/skelterjohn/go.matrix"
)

// MVNormalPDF returns the PDF of the Multivariate normal distribution. 
func MVNormalPDF(μ *DenseMatrix, Σ *DenseMatrix) func(x *DenseMatrix) float64 {
	p := μ.Rows()
	backμ := μ.DenseMatrix()
	backμ.Scale(-1)

	Σdet := Σ.Det()
	ΣdetRt := sqrt(Σdet)
	Σinv, _ := Σ.Inverse()

	normalization := pow(2*π, -float64(p)/2) / ΣdetRt

	return func(x *DenseMatrix) float64 {
		δ, _ := x.PlusDense(backμ)
		tmp := δ.Transpose()
		tmp, _ = tmp.TimesDense(Σinv)
		tmp, _ = tmp.TimesDense(δ)
		f := tmp.Get(0, 0)
		return normalization * exp(-f/2)
	}
}

// MVNormalNext returns random number drawn from the Multivariate normal distribution. 
func MVNormalNext(μ *DenseMatrix, Σ *DenseMatrix) *DenseMatrix {
	n := μ.Rows()
	x := Zeros(n, 1)
	for i := 0; i < n; i++ {
		x.Set(i, 0, NormalNext(0, 1))
	}
	C, err := Σ.Cholesky()
	Cx, err := C.TimesDense(x)
	μCx, err := μ.PlusDense(Cx)
	if err != nil {
		panic(err)
	}
	return μCx
}

// MVNormal returns the random number generator with  Multivariate normal distribution. 
func MVNormal(μ *DenseMatrix, Σ *DenseMatrix) func() *DenseMatrix {
	C, _ := Σ.Cholesky()
	n := μ.Rows()
	return func() *DenseMatrix {
		x := Zeros(n, 1)
		for i := 0; i < n; i++ {
			x.Set(i, 0, NormalNext(0, 1))
		}
		Cx, _ := C.TimesDense(x)
		MCx, _ := μ.PlusDense(Cx)
		return MCx
	}
}

// MVNormalMean returns the mean of the Multivariate normal distribution. 
func MVNormalMean(μ *DenseMatrix, Σ *DenseMatrix) *DenseMatrix {
	return μ
}

// MVNormalMode returns the mode of the Multivariate normal distribution. 
func MVNormalMode(μ *DenseMatrix, Σ *DenseMatrix) *DenseMatrix {
	return μ
}

// MVNormalVar returns the variance of the Multivariate normal distribution. 
func MVNormalVar(μ *DenseMatrix, Σ *DenseMatrix) *DenseMatrix {
	return Σ
}
