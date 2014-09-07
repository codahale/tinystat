// Copyright 2012 - 2013 The Fn Authors. All rights reserved. See the LICENSE file.

package fn

// Computes the log gamma correction factor for x >= 10 so that
// 
// log(x) -x + lgammacor(x)
// 
// [ lgammacor(x) is called	Del(x)	in other contexts (e.g. dcdflib)]
// 
// NOTES
// 
// This routine is a translation into C of a Fortran subroutine
// written by W. Fullerton of Los Alamos Scientific Laboratory.
// 
// SEE ALSO
// 
// similar in spirit,
// is faster and cleaner, but is only defined "fast" for half integers.

// lgammacor returns the log gamma correction factor for x >= 10
func lgammacor(x float64) float64 {

	algmcs := []float64{
		+.1666389480451863247205729650822e+0,
		-.1384948176067563840732986059135e-4,
		+.9810825646924729426157171547487e-8,
		-.1809129475572494194263306266719e-10,
		+.6221098041892605227126015543416e-13,
		-.3399615005417721944303330599666e-15,
		+.2683181998482698748957538846666e-17,
		-.2868042435334643284144622399999e-19,
		+.3962837061046434803679306666666e-21,
		-.6831888753985766870111999999999e-23,
		+.1429227355942498147573333333333e-24,
		-.3547598158101070547199999999999e-26,
		+.1025680058010470912000000000000e-27,
		-.3401102254316748799999999999999e-29,
		+.1276642195630062933333333333333e-30,
	}

	/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
	 *   xbig = 2 ^ 26.5
	 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
	const (
		nalgm = 5
		xbig  = 94906265.62425156
		xmax  = 3.745194030963158e306
	)

	if x < 10 {
		return nan
	} else if x >= xmax {
		panic("underflow")
	} else if x < xbig {
		tmp := 10 / x
		return chebyshev_eval(tmp*tmp*2-1, algmcs, nalgm) / x
	}
	return 1 / (x * 12)
}
