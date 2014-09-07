// Copyright 2012 - 2013 The Probab Authors. All rights reserved. See the LICENSE file.

// Variance vector for data matrix.

package stat

import (
	. "github.com/skelterjohn/go.matrix"
)

// Population variance vector of columns of data matrix, one-pass algorithm
// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
func Var(data *DenseMatrix) *Vector {
	rows := data.Rows()
	cols := data.Cols()
	out := NewVector(cols)

	for i := 0; i < cols; i++ {
		n := 0.0
		mean := 0.0
		m2 := 0.0
		for j := 0; j < rows; j++ {
			n++
			x := data.Get(j, i)
			delta := x - mean
			mean += delta / n
			if n > 1 {
				m2 += delta * (x - mean)
			}
		}
		v := m2 / n
		out.Set(i, v)
	}
	return out
}

// Sample variance vector of columns of data matrix, one-pass algorithm
// This is R:var()
func SVar(data *DenseMatrix) *Vector {
	rows := data.Rows()
	cols := data.Cols()
	n := float64(rows)
	out := Var(data)

	for i := 0; i < cols; i++ {
		v := out.Get(i) * n / (n - 1)
		out.Set(i, v)
	}
	return out
}
