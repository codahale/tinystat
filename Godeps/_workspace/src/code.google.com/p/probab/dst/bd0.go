// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// bd0 the "deviance part".
func bd0(x, np float64) float64 {
	//	Evaluates the "deviance part"
	//	bd0(x,M) :=  M // D0(x/M) = M//[ x/M // log(x/M) + 1 - (x/M) ] =
	//		  =  x // log(x/M) + M - x
	//	where M = E[X] = n//p (or = lambda), for	  x, M > 0
	//
	//	in a manner that should be stable (with small relative error)
	//	for all x and M=np. In particular for x/np close to 1, direct
	//	evaluation fails, and evaluation is based on the Taylor series
	//	of log((1+v)/(1-v)) with v = (x-np)/(x+np).
	// Rewritten from the C code by Catherine Loader, catherine@research.bell-labs.com. October 23, 2000.
	if isInf(x, 0) || isInf(np, 0) || np == 0 {
		return NaN
	}
	if abs(x-np) < 0.1*(x+np) {
		v := (x - np) / (x + np)
		s := (x - np) * v
		ej := 2 * x * v
		v = v * v
		for j := 1; ; j++ { // Taylor series
			ej *= v
			s1 := s + ej/(float64(j<<1)+1)
			if s1 == s { //	last term was effectively 0
				return s1
			}
			s = s1
		}
	}
	// else:  | x - np |  is not too small
	return x*log(x/np) + np - x
}
