package tinystat_test

import (
	"testing"

	"github.com/codahale/gubbins/assert"
	"github.com/codahale/tinystat"
	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestSummarizeOdd(t *testing.T) {
	t.Parallel()

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
	t.Parallel()

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
	t.Parallel()

	a := tinystat.Summarize([]float64{1, 2, 3, 4})
	b := tinystat.Summarize([]float64{1, 2, 3, 4})
	d := tinystat.Compare(a, b, 80)

	assert.Equal(t, "Compare",
		tinystat.Difference{
			Effect:        0,
			EffectSize:    0,
			CriticalValue: 1.31431116679138120,
			PValue:        1,
			Alpha:         0.19999999999999996,
			Beta:          0,
		},
		d, epsilon)
	assert.Equal(t, "Significant", false, d.Significant())
}

func TestCompareDifferentData(t *testing.T) {
	t.Parallel()

	a := tinystat.Summarize([]float64{1, 2, 3, 4})
	b := tinystat.Summarize([]float64{10, 20, 30, 40})
	d := tinystat.Compare(a, b, 80)

	assert.Equal(t, "Compare",
		tinystat.Difference{
			Effect:        22.5,
			EffectSize:    2.452519415855564,
			CriticalValue: 10.568344341563606,
			PValue:        0.03916791618893338,
			Alpha:         0.19999999999999996,
			Beta:          0.9856216842773273,
		},
		d, epsilon)
	assert.Equal(t, "Significant", true, d.Significant())
}

var epsilon = cmpopts.EquateApprox(0.001, 0.001) //nolint:gochecknoglobals // testing
