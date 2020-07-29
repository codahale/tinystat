// Package tinystat provides the ability to compare data sets using Student's t-test at various
// levels of confidence.
package tinystat

import (
	"math"

	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

// A Summary is a statistical summary of a normally distributed data set.
type Summary struct {
	N        float64 // N is the number of measurements in the set.
	Mean     float64 // Mean is the arithmetic mean of the measurements.
	Variance float64 // Variance is the sample variance of the data set.
}

// Summarize analyzes the given data set and returns a Summary.
func Summarize(data []float64) Summary {
	m, v := stat.MeanVariance(data, nil)

	return Summary{Mean: m, Variance: v, N: float64(len(data))}
}

// Difference represents the statistical difference between two Summary values.
type Difference struct {
	Delta    float64 // Delta is the difference between the samples' means.
	Error    float64 // Error is the margin of error at the given confidence level.
	RelDelta float64 // RelDelta is the ratio of Delta to the control mean.
	RelError float64 // RelError is the ratio of Error to the control mean.
	StdDev   float64 // StdDev is the pooled standard deviation of the two samples.
}

// Significant returns true if the difference is statistically significant.
func (d Difference) Significant() bool {
	return d.Delta > d.Error
}

// Compare returns the statistical difference between the two summaries using a two-tailed Student's
// t-test. The confidence level must be in the range (0, 100).
func Compare(control, experiment Summary, confidence float64) Difference {
	a, b := control, experiment

	// Calculate the combined degrees of freedom.
	nu := a.N + b.N - 2

	// Calculate the t-value using Student's T. gonum's implementation requires location and scale
	// parameters (mu and sigma), in addition to the degrees of freedom for the distribution. The
	// quantile is relaxed by half because this is a two-tailed test--that is, we're looking for
	// whether experimental measurements are either higher or lower than the control measurements.
	t := distuv.StudentsT{Mu: 0, Sigma: 1, Nu: nu}.
		Quantile(1 - ((1 - (confidence / 100)) / 2))

	// Calculate the pooled estimate of standard deviation.
	s := math.Sqrt(((a.N-1)*a.Variance + (b.N-1)*b.Variance) / nu)

	// Calculate the difference between the means of the two samples.
	d := math.Abs(a.Mean - b.Mean)

	// Calculate the standard error, given the t-value and standard deviation.
	e := t * s * math.Sqrt(1.0/a.N+1.0/b.N)

	return Difference{
		Delta:    d,
		Error:    e,
		RelDelta: d / control.Mean,
		RelError: e / control.Mean,
		StdDev:   s,
	}
}
