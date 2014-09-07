// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Planck distribution. 
// Johnson and Kotz, 1970
//
// Parameters: 
// a > 0.0	 	
// b > 0.0		
//
// Support: 
// x > 0.0

// PlanckPDF returns the PDF of the Planck distribution. 
func PlanckPDF(a, b float64) func(x float64) float64 {
	// ζ() waiting for better implementation
	ζ := ζ
	return func(x float64) float64 {
		t1 := pow(b, a+1)
		t2 := pow(x, a)
		t3 := Γ(a+1) * ζ(a+1)
		t4 := exp(b*x) - 1
		p := (t1 * t2) / (t3 * t4)
		return p
	}
}

// PlanckNext returns random number drawn from the Planck distribution. 
// Devroye 1986: 552.
// Devroye, L. 1986: Non-Uniform Random Variate Generation. Springer-Verlag, New York. ISBN 0-387-96305-7.
func PlanckNext(a, b float64) (x float64) {
	g := GammaNext(a+1, 1) // OK, consulted with Luc Devroye
	z := float64(ZetaNext(a + 1))
	return g / (b * z)
}

// Planck returns the random number generator with  Planck distribution. 
func Planck(a, b float64) func() float64 {
	return func() float64 { return PlanckNext(a, b) }
}
