// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Inverse-Wishart distribution (not to be confused with Inverse CDF of Wishart distribution). 
// The inverse Wishart distribution, also called the inverted Wishart distribution, is a probability distribution defined on real-valued positive-definite matrices. In Bayesian statistics it is used as the conjugate prior for the covariance matrix of a multivariate normal distribution.
//
// Parameters: 
// n > p-1	 	degrees of freedom (real)
// Ψ >0			inverse scale matrix (positive definite)
//
// Support: 
// X is positive definite

import (
	m "github.com/skelterjohn/go.matrix"
)

// InverseWishartPDF returns the PDF of the Inverse-Wishart distribution. 
func InverseWishartPDF(n int, Ψ *m.DenseMatrix) func(B *m.DenseMatrix) float64 {
	p := Ψ.Rows()
	Ψdet := Ψ.Det()
	normalization := pow(Ψdet, -0.5*float64(n)) *
		pow(2, -0.5*float64(n*p)) /
		Γ(float64(n)/2)
	return func(B *m.DenseMatrix) float64 {
		Bdet := B.Det()
		Binv, _ := B.Inverse()
		ΨBinv, _ := Ψ.Times(Binv)
		return normalization *
			pow(Bdet, -.5*float64(n+p+1)) *
			exp(-0.5*ΨBinv.Trace())
	}
}

// InverseWishartLnPDF returns the natural logarithm of the PDF of the Inverse-Wishart distribution. 
func InverseWishartLnPDF(n int, Ψ *m.DenseMatrix) func(W *m.DenseMatrix) float64 {
	p := Ψ.Rows()
	Ψdet := Ψ.Det()
	normalization := log(Ψdet)*-0.5*float64(n) +
		log(2)*-0.5*float64(n*p) -
		LnΓ(float64(n)/2)
	return func(B *m.DenseMatrix) float64 {
		Bdet := B.Det()
		Binv, _ := B.Inverse()
		ΨBinv, _ := Ψ.Times(Binv)
		return normalization +
			log(Bdet)*-.5*float64(n+p+1) +
			-0.5*ΨBinv.Trace()
	}
}

// InverseWishartNext returns random number drawn from the Inverse-Wishart distribution. 
func InverseWishartNext(n int, V *m.DenseMatrix) *m.DenseMatrix {
	return InverseWishart(n, V)()
}

// InverseWishart returns the random number generator with  Inverse-Wishart distribution. 
func InverseWishart(n int, V *m.DenseMatrix) func() *m.DenseMatrix {
	p := V.Rows()
	zeros := m.Zeros(p, 1)
	rowGen := MVNormal(zeros, V)
	return func() *m.DenseMatrix {
		x := make([][]float64, n)
		for i := 0; i < n; i++ {
			x[i] = rowGen().Array()
		}
		X := m.MakeDenseMatrixStacked(x)
		S, _ := X.Transpose().TimesDense(X)
		Sinv, _ := S.Inverse()
		return Sinv
	}
}
