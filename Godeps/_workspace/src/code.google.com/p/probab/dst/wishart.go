// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Wishart distribution. 
// A generalization to multiple dimensions of the chi-squared distribution, or, in the case of non-integer degrees of freedom, of the gamma distribution. 
// Wishart, J. (1928). "The generalised product moment distribution in samples from a normal multivariate population". Biometrika 20A (1-2): 32–52. doi:10.1093/biomet/20A.1-2.32. 
// Parameters: 
// n > p-1	degrees of freedom (real)
// V > 0		pxp scale matrix	(positive definite, real)
//
// Support: 
// X 	pxp positive definite, real

import (
	m "github.com/skelterjohn/go.matrix"
)

// WishartPDF returns the PDF of the Wishart distribution. 
func WishartPDF(n int, V *m.DenseMatrix) func(W *m.DenseMatrix) float64 {
	p := V.Rows()
	Vdet := V.Det()
	Vinv, _ := V.Inverse()
	normalization := pow(2, -0.5*float64(n*p)) *
		pow(Vdet, -0.5*float64(n)) /
		Γ(0.5*float64(n))
	return func(W *m.DenseMatrix) float64 {
		VinvW, _ := Vinv.Times(W)
		return normalization * pow(W.Det(), 0.5*float64(n-p-1)) *
			exp(-0.5*VinvW.Trace())
	}
}

// WishartLnPDF returns the natural logarithm of the PDF of the Wishart distribution. 
func WishartLnPDF(n int, V *m.DenseMatrix) func(W *m.DenseMatrix) float64 {

	p := V.Rows()
	Vdet := V.Det()
	Vinv, _ := V.Inverse()
	normalization := log(2)*(-0.5*float64(n*p)) +
		log(Vdet)*(-0.5*float64(n)) -
		LnΓ(0.5*float64(n))
	return func(W *m.DenseMatrix) float64 {
		VinvW, _ := Vinv.Times(W)
		return normalization +
			log(W.Det())*0.5*float64(n-p-1) -
			0.5*VinvW.Trace()
	}
}

// WishartNext returns random number drawn from the Wishart distribution. 
func WishartNext(n int, V *m.DenseMatrix) *m.DenseMatrix {
	return Wishart(n, V)()
}

// Wishart returns the random number generator with  Wishart distribution. 
func Wishart(n int, V *m.DenseMatrix) func() *m.DenseMatrix {
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
		return S
	}
}
