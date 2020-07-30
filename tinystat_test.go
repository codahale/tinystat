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
	assert.Equal(t, "StdDev", 1.0, s.StdDev(), epsilon)
	assert.Equal(t, "StdErr", 0.5773502691896258, s.StdErr(), epsilon)
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
	assert.Equal(t, "StdDev", 1.2909944487358056, s.StdDev(), epsilon)
	assert.Equal(t, "StdErr", 0.6454972243679028, s.StdErr(), epsilon)
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
			P:                0.20000000000000004,
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
			CriticalValue:    10.568344341563606,
			RelDelta:         9.0,
			RelCriticalValue: 4.227337736625442,
			P:                0.19999999999999996,
		},
		d, epsilon)
	assert.Equal(t, "Significant", true, d.Significant())
}

var epsilon = cmpopts.EquateApprox(0.001, 0.001) //nolint:gochecknoglobals // testing
