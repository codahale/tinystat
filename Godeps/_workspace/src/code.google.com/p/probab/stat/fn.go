// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package stat

import (
	fn "code.google.com/p/go-fn/fn"
	"math"
)

const π = float64(math.Pi)
const ln2 = math.Ln2
const lnSqrt2π = 0.918938533204672741780329736406 // log(sqrt(2*pi))
const min64 = math.SmallestNonzeroFloat64         //   DBL_MIN
const eps64 = 1.1102230246251565e-16              // DBL_EPSILON   
const maxExp = 1024.0                             // DBL_MAX_EXP
const sqrt2 = math.Sqrt2

var nan = math.NaN()

var fZero float64 = float64(0.0)
var fOne float64 = float64(1.0)
var iZero int64 = int64(0)
var iOne int64 = int64(1)

var negInf float64 = math.Inf(-1)
var posInf float64 = math.Inf(+1)

// Functions imported from "math"
var abs func(float64) float64 = math.Abs
var floor func(float64) float64 = math.Floor
var ceil func(float64) float64 = math.Ceil
var log func(float64) float64 = math.Log
var log1p func(float64) float64 = math.Log1p
var log10 func(float64) float64 = math.Log10
var exp func(float64) float64 = math.Exp
var sqrt func(float64) float64 = math.Sqrt
var pow func(float64, float64) float64 = math.Pow
var atan func(float64) float64 = math.Atan
var tan func(float64) float64 = math.Tan
var trunc func(float64) float64 = math.Trunc
var erf func(float64) float64 = math.Erf
var erfc func(float64) float64 = math.Erfc
var isNaN func(float64) bool = math.IsNaN
var isInf func(float64, int) bool = math.IsInf

// Functions imported from "code.google.com/p/go-fn/fn"
var lnB func(float64, float64) float64 = fn.LnB
var lnΓ func(float64) float64 = fn.LnΓ

//  vAbs recalculates the data vector to absolute values.
func vAbs(x []float64) {
	for i, val := range x {
		x[i] = abs(val)
	}
}

//  vCent recalculates the data vector to centralized values.
func vCent(x []float64) {
	mu := Mean(x)
	for i, val := range x {
		x[i] = val - mu
	}
}

//  vPow recalculates the data vector to absolute values.
func vPow(x []float64, power float64) {
	for i, val := range x {
		x[i] = pow(val, power)
	}
}

//  sum returns the sum of the data vector.
func sum(x []float64) float64 {
	s := 0.0
	for _, val := range x {
		s += val
	}
	return s
}

// mean returns the mean of the data vector.
func mean(x []float64) float64 {
	μ := sum(x)
	μ /= float64(len(x))
	return μ
}

// diffMean returns the vector of centralized values.
func diffMean(x []float64) []float64 {
	d := make([]float64, len(x))
	mu := Mean(x)
	for i, val := range x {
		d[i] = val - mu
	}
	return d
}
