package tinystat_test

import (
	"testing"

	"github.com/codahale/tinystat"
)

func TestSummarizeOdd(t *testing.T) {
	s := tinystat.Summarize([]float64{1, 2, 3})

	if v, want := s.N, 3.0; v != want {
		t.Errorf("N was %v, but expected %v", v, want)
	}

	if v, want := s.Min, 1.0; v != want {
		t.Errorf("Min was %v, but expected %v", v, want)
	}

	if v, want := s.Max, 3.0; v != want {
		t.Errorf("Max was %v, but expected %v", v, want)
	}

	if v, want := s.Median, 2.0; v != want {
		t.Errorf("Median was %v, but expected %v", v, want)
	}

	if v, want := s.Mean, 2.0; v != want {
		t.Errorf("Mean was %v, but expected %v", v, want)
	}

	if v, want := s.Variance, 1.0; v != want {
		t.Errorf("Variance was %v, but expected %v", v, want)
	}

	if v, want := s.StdDev, 1.0; v != want {
		t.Errorf("StdDev was %v, but expected %v", v, want)
	}
}

func TestSummarizeEven(t *testing.T) {
	s := tinystat.Summarize([]float64{1, 2, 3, 4})

	if v, want := s.N, 4.0; v != want {
		t.Errorf("N was %v, but expected %v", v, want)
	}

	if v, want := s.Min, 1.0; v != want {
		t.Errorf("Min was %v, but expected %v", v, want)
	}

	if v, want := s.Max, 4.0; v != want {
		t.Errorf("Max was %v, but expected %v", v, want)
	}

	if v, want := s.Median, 2.5; v != want {
		t.Errorf("Median was %v, but expected %v", v, want)
	}

	if v, want := s.Mean, 2.5; v != want {
		t.Errorf("Mean was %v, but expected %v", v, want)
	}

	if v, want := s.Variance, 1.6666666666666667; v != want {
		t.Errorf("Variance was %v, but expected %v", v, want)
	}

	if v, want := s.StdDev, 1.2909944487358056; v != want {
		t.Errorf("StdDev was %v, but expected %v", v, want)
	}
}

func TestCompareSimilarData(t *testing.T) {
	a := tinystat.Summarize([]float64{1, 2, 3, 4})
	b := tinystat.Summarize([]float64{1, 2, 3, 4})
	d := tinystat.Compare(a, b, 80)

	if v, want := d.Delta, 0.0; v != want {
		t.Errorf("Delta was %v, but expected %v", v, want)
	}

	if v, want := d.Error, 1.314311166777796; v != want {
		t.Errorf("Error was %v, but expected %v", v, want)
	}

	if v, want := d.StdDev, 1.2909944487358056; v != want {
		t.Errorf("StdDev was %v, but expected %v", v, want)
	}

	if v, want := d.Significant(), false; v != want {
		t.Errorf("Significance was %v, but expected %v", v, want)
	}
}

func TestCompareDifferentData(t *testing.T) {
	a := tinystat.Summarize([]float64{1, 2, 3, 4})
	b := tinystat.Summarize([]float64{10, 20, 30, 40})
	d := tinystat.Compare(a, b, 80)

	if v, want := d.Delta, 22.5; v != want {
		t.Errorf("Delta was %v, but expected %v", v, want)
	}

	if v, want := d.Error, 9.33993571056027; v != want {
		t.Errorf("Error was %v, but expected %v", v, want)
	}

	if v, want := d.StdDev, 9.17423929634859; v != want {
		t.Errorf("StdDev was %v, but expected %v", v, want)
	}

	if v, want := d.Significant(), true; v != want {
		t.Errorf("Significance was %v, but expected %v", v, want)
	}
}
