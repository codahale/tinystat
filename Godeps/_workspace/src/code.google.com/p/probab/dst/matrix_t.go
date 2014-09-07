package dst

import (
	"fmt"
	mx "github.com/skelterjohn/go.matrix"
)

func checkMatrixT(M, Omega, Sigma *mx.DenseMatrix, n int) {
	p := M.Rows()
	m := M.Cols()
	if Omega.Rows() != p {
		panic(fmt.Sprintf("Omega.Rows != M.Rows, %d != %d", Omega.Rows(), p))
	}
	if Omega.Cols() != Omega.Rows() {
		panic("Omega is not square")
	}
	if Sigma.Rows() != m {
		panic(fmt.Sprintf("Sigma.Cols != M.Cols, %d != %d", Sigma.Cols(), m))
	}
	if Sigma.Cols() != Sigma.Rows() {
		panic("Sigma is not square")
	}
	if n <= 0 {
		panic("n <= 0")
	}
}

func MatrixTPDF(M, Omega, Sigma *mx.DenseMatrix, n int) func(T *mx.DenseMatrix) (l float64) {
	checkMatrixT(M, Omega, Sigma, n)

	nf := float64(n)
	p := M.Rows()
	pf := float64(p)
	m := M.Cols()
	mf := float64(m)

	var norm float64 = 1

	norm *= Γpr(p, 0.5*(nf+mf+pf-1), 0.5*(nf+pf-1))
	norm *= pow(π, -0.5*mf*pf)
	norm *= pow(Omega.Det(), -0.5*mf)
	norm *= pow(Sigma.Det(), -0.5*pf)

	SigmaInv, err := Sigma.Inverse()
	if err != nil {
		panic(err)
	}
	OmegaInv, err := Omega.Inverse()
	if err != nil {
		panic(err)
	}

	return func(T *mx.DenseMatrix) (l float64) {
		l = norm

		diff, err := T.MinusDense(M)
		if err != nil {
			panic(err)
		}
		inner := OmegaInv.Copy()
		inner, _ = inner.TimesDense(diff)
		inner, _ = inner.TimesDense(SigmaInv)
		inner, _ = inner.TimesDense(diff.Transpose())

		l *= pow(inner.Det(), -0.5*(nf+mf+pf-1))

		return
	}
}

func MatrixTLnPDF(M, Omega, Sigma *mx.DenseMatrix, n int) func(T *mx.DenseMatrix) (ll float64) {
	checkMatrixT(M, Omega, Sigma, n)

	nf := float64(n)
	p := M.Rows()
	pf := float64(p)
	m := M.Cols()
	mf := float64(m)

	var norm float64 = 0

	norm += logΓpr(p, 0.5*(nf+mf+pf-1), 0.5*(nf+pf-1))
	norm += pow(π, -0.5*mf*pf)
	norm += pow(Omega.Det(), -0.5*mf)
	norm += pow(Sigma.Det(), -0.5*pf)

	SigmaInv, err := Sigma.Inverse()
	if err != nil {
		panic(err)
	}
	OmegaInv, err := Omega.Inverse()
	if err != nil {
		panic(err)
	}

	return func(T *mx.DenseMatrix) (ll float64) {
		ll = norm

		diff, err := T.MinusDense(M)
		if err != nil {
			panic(err)
		}
		inner := OmegaInv.Copy()
		inner, _ = inner.TimesDense(diff)
		inner, _ = inner.TimesDense(SigmaInv)
		inner, _ = inner.TimesDense(diff.Transpose())

		ll += log(inner.Det()) * -0.5 * (nf + mf + pf - 1)

		return
	}
}

func MatrixT(M, Omega, Sigma *mx.DenseMatrix, n int) func() (T *mx.DenseMatrix) {
	checkMatrixT(M, Omega, Sigma, n)

	fmt.Println("M:", M)
	fmt.Println("Sigma:", Sigma)
	fmt.Println("Omega:", Omega)

	p := M.Rows()
	m := M.Cols()

	OmegaInv, err := Omega.Inverse()
	if err != nil {
		panic(err)
	}

	Sdist := Wishart(n+p-1, OmegaInv)

	Xdist := MatrixNormal(mx.Zeros(p, m), mx.Eye(p), Sigma)

	return func() (T *mx.DenseMatrix) {
		S := Sdist()
		Sinv, err := S.Inverse()
		if err != nil {
			panic(err)
		}
		Sinvc, err := Sinv.Cholesky()
		if err != nil {
			panic(err)
		}
		X := Xdist()
		fmt.Println("Sinvc:", Sinvc)
		fmt.Println("X:", X)
		T, err = Sinvc.Transpose().TimesDense(X)
		if err != nil {
			panic(err)
		}
		err = T.AddDense(M)
		if err != nil {
			panic(err)
		}
		return
	}
}

func MatrixTNext(M, Omega, Sigma *mx.DenseMatrix, n int) (T *mx.DenseMatrix) {
	return MatrixT(M, Omega, Sigma, n)()
}
