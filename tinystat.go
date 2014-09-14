// Package tinystat provides the ability to compare data sets using Student's
// t-test at various levels of confidence.
package tinystat

import (
	"math"
	"sort"

	"code.google.com/p/probab/dst"
)

// A Summary is a statistical summary of a normally distributed data set.
type Summary struct {
	N        float64 // N is the number of measurements in the set.
	Min      float64 // Min is the smallest measurement.
	Max      float64 // Max is the largest measurement.
	Median   float64 // Median is the median measurement.
	Mean     float64 // Mean is the arithmetic mean of the measurements.
	Variance float64 // Variance is the sample variance of the data set.
	StdDev   float64 // StdDev is the sample standard deviation of the data set.
}

// Summarize analyzes the given data set and returns a Summary.
func Summarize(data []float64) Summary {
	min, max, mean, variance := minMaxMeanVar(data)
	median := median(data)

	return Summary{
		Min:      min,
		Max:      max,
		Mean:     mean,
		Median:   median,
		Variance: variance,
		StdDev:   math.Sqrt(variance),
		N:        float64(len(data)),
	}
}

// Difference represents the statistical difference between two samples.
type Difference struct {
	Delta        float64
	Error        float64
	PctDelta     float64
	PctError     float64
	PooledStdDev float64
}

// Significant returns true if the difference is statistically significant.
func (d Difference) Significant() bool {
	return d.Delta > d.Error
}

// Compare returns the statistical difference between the two summaries using
// Student's t-test. The confidence level must be in the range (0, 99).
func Compare(a, b Summary, confidence float64) Difference {
	// calculate the quantile for two-sided Student's t
	t := dst.StudentsTQtlFor(a.N+b.N-2, 1-((1-(confidence/100))/2))

	s := math.Sqrt(
		((a.N-1)*a.Variance + (b.N-1)*b.Variance) /
			(a.N + b.N - 2),
	)
	d := math.Abs(a.Mean - b.Mean)
	e := t * s * math.Sqrt(1.0/a.N+1.0/b.N)

	return Difference{
		Delta:        d,
		Error:        e,
		PctDelta:     d * 100 / b.Mean,
		PctError:     e * 100 / b.Mean,
		PooledStdDev: s,
	}
}

func minMaxMeanVar(data []float64) (min float64, max float64, mean float64, variance float64) {
	min, max = math.Inf(1), math.Inf(-1)

	for n, x := range data {
		if x < min {
			min = x
		}

		if x > max {
			max = x
		}

		// Welford algorithm for corrected variance
		delta := x - mean
		mean += delta / float64(n+1)
		variance += delta * (x - mean)
	}

	variance /= float64(len(data) - 1) // Bessel's correction

	return
}

func median(data []float64) float64 {
	// don't mutate the argument
	d := make([]float64, len(data))
	copy(d, data)

	sort.Float64s(d)

	if len(d)%2 == 1 {
		return d[len(d)/2]
	}
	return (d[len(d)/2-1] + d[len(d)/2]) / 2

}
