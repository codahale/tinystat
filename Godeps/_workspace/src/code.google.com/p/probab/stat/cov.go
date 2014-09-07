// Covariance matrix 

package stat

import (
	. "github.com/skelterjohn/go.matrix"
)

// Covariance matrix between columns of data matrix, two-pass algorithm
// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Covariance
func Cov(data *DenseMatrix) *DenseMatrix {
	rows := data.Rows()
	cols := data.Cols()
	out := Zeros(cols, cols)

	for i := 0; i < cols; i++ {
		for j := i; j < cols; j++ {
			meanX := 0.0
			meanY := 0.0
			// calculate column means
			for k := 0; k < rows; k++ {
				meanX += data.Get(k, i)
				meanY += data.Get(k, j)
			}
			meanX /= float64(rows)
			meanY /= float64(rows)
			// calculate covariance
			cov := 0.0
			for k := 0; k < rows; k++ {
				x := data.Get(k, i)
				y := data.Get(k, j)
				cov += (x - meanX) * (y - meanY)
			}
			cov /= float64(rows)
			out.Set(i, j, cov)
			out.Set(j, i, cov)
		}
	}
	return out
}

// Sample covariance matrix between columns of data matrix, for samples
// This is R:cov()
func SCov(data *DenseMatrix) *DenseMatrix {
	rows := data.Rows()
	cols := data.Cols()
	out := Cov(data)

	for i := 0; i < cols; i++ {
		for j := i; j < cols; j++ {
			v := out.Get(i, j) * float64(rows) / (float64(rows) - 1)
			out.Set(i, j, v)
			out.Set(j, i, v)
		}
	}
	return out
}
