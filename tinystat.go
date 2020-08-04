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
	// Effect is the absolute difference between the samples' means.
	Effect float64

	// EffectSize is the difference in means between the two samples, normalized for variance.
	// Technically, this is Cohen's d.
	EffectSize float64

	// CriticalValue is the maximum allowed Effect at the given confidence level.
	CriticalValue float64

	// PValue is the p-value for the test. This is, effectively, the probability that accepting the
	// results of this test will be a Type 1 error, in which the null hypothesis (i.e. there is no
	// difference between the means of the two samples) will be rejected when it is in fact true.
	PValue float64

	// Alpha is the significance level of the test. It is the maximum allowed value of the p-value.
	Alpha float64

	// Beta is the probability of a Type 2 error. That is, the probability that the null hypothesis
	// will be retained despite it not being true.
	Beta float64
}

// Significant returns true if the difference is statistically significant.
func (d Difference) Significant() bool {
	return d.Effect > d.CriticalValue
}

// Compare returns the statistical difference between the two summaries using a two-tailed Welch's
// t-test. The confidence level must be in the range (0, 100).
func Compare(control, experiment Summary, confidence float64) Difference {
	a, b := control, experiment

	// Calculate the significance level.
	alpha := 1 - (confidence / 100)

	// Calculate the degrees of freedom.
	nu := math.Pow(a.Variance/a.N+b.Variance/b.N, 2) /
		(math.Pow(a.Variance, 2)/(math.Pow(a.N, 2)*(a.N-1)) +
			math.Pow(b.Variance, 2)/(math.Pow(b.N, 2)*(b.N-1)))

	// Create a Student's T distribution with location of 0, a scale of 1, and a shape of the number
	// of degrees of freedom in the test.
	studentsT := distuv.StudentsT{Mu: 0, Sigma: 1, Nu: nu}

	// Calculate the hypothetical two-tailed t-value for the given significance level.
	tHyp := studentsT.Quantile(1 - (alpha / tails))

	// Calculate the absolute difference between the means of the two samples.
	d := math.Abs(a.Mean - b.Mean)

	// Calculate the standard error.
	s := math.Sqrt(a.Variance/a.N + b.Variance/b.N)

	// Calculate the experimental t-value.
	tExp := d / s

	// Calculate the p-value given the experimental t-value.
	p := studentsT.CDF(-tExp) * tails

	// Calculate the critical value.
	cv := tHyp * s

	// Calculate the standard deviating using mean variance.
	sd := math.Sqrt((a.Variance + b.Variance) / 2)

	// Calculate Cohen's d for the effect size.
	cd := d / sd

	// Create a standard normal distribution.
	stdNormal := distuv.UnitNormal

	// Calculate the statistical power.
	z := d / (sd * math.Sqrt(1/a.N+1/b.N))
	za := stdNormal.Quantile(1 - alpha/tails)
	beta := stdNormal.CDF(z-za) - stdNormal.CDF(-z-za)

	return Difference{
		Effect:        d,
		CriticalValue: cv,
		EffectSize:    cd,
		PValue:        p,
		Alpha:         alpha,
		Beta:          beta,
	}
}

// tails is the number of distribution tails used to determine significance. In this case, we always
// use a two-tailed test because our null hypothesis is that the samples are not different.
const tails = 2
