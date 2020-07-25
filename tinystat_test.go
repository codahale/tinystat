package tinystat_test

import (
	"testing"

	"github.com/codahale/tinystat"
	"github.com/codahale/tinystat/internal/assert"
	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestSummarizeOdd(t *testing.T) {
	s := tinystat.Summarize([]float64{1, 2, 3})

	assert.Equal(t, "N", 3.0, s.N, epsilon)
	assert.Equal(t, "Mean", 2.0, s.Mean, epsilon)
	assert.Equal(t, "Variance", 1.0, s.Variance, epsilon)
}

func TestSummarizeEven(t *testing.T) {
	s := tinystat.Summarize([]float64{1, 2, 3, 4})

	assert.Equal(t, "N", 4.0, s.N, epsilon)
	assert.Equal(t, "Mean", 2.5, s.Mean, epsilon)
	assert.Equal(t, "Variance", 1.6666, s.Variance, epsilon)
}

func TestCompareSimilarData(t *testing.T) {
	a := tinystat.Summarize([]float64{1, 2, 3, 4})
	b := tinystat.Summarize([]float64{1, 2, 3, 4})
	d := tinystat.Compare(a, b, 80)

	assert.Equal(t, "Delta", 0.0, d.Delta, epsilon)
	assert.Equal(t, "Error", 1.314, d.Error, epsilon)
	assert.Equal(t, "RelDelta", 0.0, d.RelDelta, epsilon)
	assert.Equal(t, "RelError", 0.5257, d.RelError, epsilon)
	assert.Equal(t, "StdDev", 1.29099, d.StdDev, epsilon)
	assert.Equal(t, "Significant", false, d.Significant())
}

func TestCompareDifferentData(t *testing.T) {
	a := tinystat.Summarize([]float64{1, 2, 3, 4})
	b := tinystat.Summarize([]float64{10, 20, 30, 40})
	d := tinystat.Compare(a, b, 80)

	assert.Equal(t, "Delta", 22.5, d.Delta, epsilon)
	assert.Equal(t, "Error", 9.33993, d.Error, epsilon)
	assert.Equal(t, "RelDelta", 9.0, d.RelDelta, epsilon)
	assert.Equal(t, "RelError", 3.7359, d.RelError, epsilon)
	assert.Equal(t, "StdDev", 9.1742, d.StdDev, epsilon)
	assert.Equal(t, "Significant", true, d.Significant())
}

var epsilon = cmpopts.EquateApprox(0.001, 0.001)
