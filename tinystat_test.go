package tinystat_test

import (
	"testing"

	"github.com/codahale/tinystat"
	"github.com/codahale/tinystat/internal/assert"
	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestSummarizeOdd(t *testing.T) {
	s := tinystat.Summarize([]float64{1, 2, 3})

	assert.Equal(t, "Summarize",
		tinystat.Summary{
			N:        3,
			Mean:     2,
			Variance: 1,
		},
		s, epsilon)
}

func TestSummarizeEven(t *testing.T) {
	s := tinystat.Summarize([]float64{1, 2, 3, 4})

	assert.Equal(t, "Summarize",
		tinystat.Summary{
			N:        4,
			Mean:     2.5,
			Variance: 1.6666666666666667,
		},
		s, epsilon)
}

func TestCompareSimilarData(t *testing.T) {
	a := tinystat.Summarize([]float64{1, 2, 3, 4})
	b := tinystat.Summarize([]float64{1, 2, 3, 4})
	d := tinystat.Compare(a, b, 80)

	assert.Equal(t, "Compare",
		tinystat.Difference{
			Delta:            0,
			CriticalValue:    1.31431116679138120,
			RelDelta:         0,
			RelCriticalValue: 0.5257244667165525,
		},
		d, epsilon)
	assert.Equal(t, "Significant", false, d.Significant())
}

func TestCompareDifferentData(t *testing.T) {
	a := tinystat.Summarize([]float64{1, 2, 3, 4})
	b := tinystat.Summarize([]float64{10, 20, 30, 40})
	d := tinystat.Compare(a, b, 80)

	assert.Equal(t, "Compare",
		tinystat.Difference{
			Delta:            22.5,
			CriticalValue:    10.624320828772332,
			RelDelta:         9.0,
			RelCriticalValue: 4.2497283315089325,
		},
		d, epsilon)
	assert.Equal(t, "Significant", true, d.Significant())
}

var epsilon = cmpopts.EquateApprox(0.001, 0.001) //nolint:gochecknoglobals // testing
