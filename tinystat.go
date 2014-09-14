// Package tinystat provides the ability to compare data sets using Student's
// t-test at various levels of confidence.
package tinystat

import (
	"math"

	"code.google.com/p/probab/dst"
)

// A Summary is a statistical summary of a normally distributed data set.
type Summary struct {
	N        float64 // N is the number of measurements in the set.
	Mean     float64 // Mean is the arithmetic mean of the measurements.
	Variance float64 // Variance is the sample variance of the data set.
}

// Summarize analyzes the given data set and returns a Summary.
func Summarize(data []float64) Summary {
	var m, m2 float64
	for n, x := range data {
		// Welford algorithm for corrected variance
		delta := x - m
		m += delta / float64(n+1)
		m2 += delta * (x - m)
	}

	return Summary{
		Mean:     m,
		Variance: m2 / float64(len(data)-1), // Bessel's correction
		N:        float64(len(data)),
	}
}

// Difference represents the statistical difference between two samples.
type Difference struct {
	Delta  float64 // Delta is the difference between the samples' means.
	Error  float64 // Error is the margin of error at the given confidence level.
	StdDev float64 // StdDev is the standard deviation of the two samples.
}

// Significant returns true if the difference is statistically significant.
func (d Difference) Significant() bool {
	return d.Delta > d.Error
}

// Compare returns the statistical difference between the two summaries using a
// two-tailed Student's t-test. The confidence level must be in the range (0,
// 100).
func Compare(a, b Summary, confidence float64) Difference {
	// calculate the quantile for two-tailed Student's t
	t := dst.StudentsTQtlFor(a.N+b.N-2, 1-((1-(confidence/100))/2))

	s := math.Sqrt(
		((a.N-1)*a.Variance + (b.N-1)*b.Variance) /
			(a.N + b.N - 2),
	)
	d := math.Abs(a.Mean - b.Mean)
	e := t * s * math.Sqrt(1.0/a.N+1.0/b.N)

	return Difference{
		Delta:  d,
		Error:  e,
		StdDev: s,
	}
}
