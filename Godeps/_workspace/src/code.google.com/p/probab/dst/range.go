package dst

import (
	"math/rand"
)

func RangePMF(n int64) func(i int64) float64 {
	return func(i int64) float64 {
		return fOne / float64(n)
	}
}
func LnRangePMF(n int64) func(i int64) float64 {
	return func(i int64) float64 {
		return -log(float64(n))
	}
}
func RangeNext(n int64) int64 {
	return rand.Int63n(n)
}
func Range(n int64) func() int64 {
	return func() int64 {
		return RangeNext(n)
	}
}
