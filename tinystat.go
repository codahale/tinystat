// Package tinystat provides the ability to compare data sets using Welch's t-test at various levels
// of confidence.
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

// StdDev returns the standard deviation of the sample.
func (s *Summary) StdDev() float64 {
	return math.Sqrt(s.Variance)
}

// StdErr returns the standard error of the sample.
func (s *Summary) StdErr() float64 {
	return stat.StdErr(s.StdDev(), s.N)
}

// Summarize analyzes the given data set and returns a Summary.
func Summarize(data []float64) Summary {
	m, v := stat.MeanVariance(data, nil)

	return Summary{Mean: m, Variance: v, N: float64(len(data))}
}

// Difference represents the statistical difference between two Summary values.
type Difference struct {
	Delta            float64 // Delta is the absolute difference between the samples' means.
	CriticalValue    float64 // CriticalValue is the maximum allowed Delta at the given confidence level.
	RelDelta         float64 // RelDelta is the ratio of Delta to the control mean.
	RelCriticalValue float64 // RelCriticalValue is the ratio of CriticalValue to the control mean.

	// PValue is the p-value for the test. This is, effectively, the probability that accepting the
	// results of this test will be a Type 1 error, in which the null hypothesis (i.e. there is no
	// difference between the means of the two samples) will be rejected when it is in fact true.
	PValue float64

	// EffectSize is the difference in means between the two samples, normalized for variance.
	// Technically, this is Cohen's d.
	EffectSize float64
}

// Significant returns true if the difference is statistically significant.
func (d Difference) Significant() bool {
	return d.Delta > d.CriticalValue
}

// Compare returns the statistical difference between the two summaries using a two-tailed Welch's
// t-test. The confidence level must be in the range (0, 100).
func Compare(control, experiment Summary, confidence float64) Difference {
	a, b := control, experiment

	// Calculate the degrees of freedom.
	nu := math.Pow(a.Variance/a.N+b.Variance/b.N, 2) /
		(math.Pow(a.Variance, 2)/(math.Pow(a.N, 2)*(a.N-1)) +
			math.Pow(b.Variance, 2)/(math.Pow(b.N, 2)*(b.N-1)))

	// Create a Student's T distribution with location of 0, a scale of 1, and a shape of the number
	// of degrees of freedom in the test.
	dist := distuv.StudentsT{Mu: 0, Sigma: 1, Nu: nu}

	// Calculate the t-value for the given confidence level.
	t := dist.Quantile(1 - ((1 - (confidence / 100)) / tails))

	// Calculate the p-value given the t-value.
	p := dist.CDF(-t) * tails

	// Calculate the difference between the means of the two samples.
	d := math.Abs(a.Mean - b.Mean)

	// Calculate the standard error.
	s := math.Sqrt(a.Variance/a.N + b.Variance/b.N)

	// Calculate the critical value.
	cv := t * s

	// Calculate Cohen's d for the effect size, using the mean variance instead of the pooled
	// variance.
	cd := d / math.Sqrt((a.Variance+b.Variance)/2)

	return Difference{
		Delta:            d,
		CriticalValue:    cv,
		RelDelta:         d / control.Mean,
		RelCriticalValue: cv / control.Mean,
		PValue:           p,
		EffectSize:       cd,
	}
}

// tails is the number of distribution tails used to determine significance. In this case, we always
// use a two-tailed test because our null hypothesis is that the samples are not different.
const tails = 2
