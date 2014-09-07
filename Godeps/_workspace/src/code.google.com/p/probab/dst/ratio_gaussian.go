// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Gaussian ratio distribution. 

import (
	"math"
)

//GearyHinkleyTransformation transforms the ratio of two normally distributed variables to the transformed variable T would approximately have a standard Gaussian distribution. See Hinkley(1969). 
func GearyHinkleyTransformation(z, μX, σX, μY, σY, ρ float64) float64 {
	//A Geary–Hinkley transformation, under certain assumptions, returns the transformed variable T that would approximately have a standard Gaussian distribution. The approximation is good if Y is unlikely to assume negative values.
	// X = N(μX, σ2X) and Y = N(μY, σ2Y) 
	// Z = X/Y
	// http://en.wikipedia.org/wiki/Ratio_distribution#A_transformation_to_Gaussianity
	σ2X := σX * σX
	σ2Y := σY * σY

	t1 := μY*z - μX
	t2 := math.Sqrt(σ2Y*z*z - 2*ρ*σX*σY*z + σ2X)
	return t1 / t2
}

func phi(x float64) float64 {
	return ((1.0 / 2.0) * (1 + erf((x)/(sqrt2))))
}

// GaussianRatioNoCorrPDFAt returns the value of PDF of Gaussian Ratio distribution of uncorrelated variables.
func GaussianRatioNoCorrPDF(μX, σX, μY, σY float64) func(z float64) float64 {
	return func(z float64) float64 {
		σ2X := σX * σX
		σ2Y := σY * σY
		μ2X := μX * μX
		μ2Y := μY * μY
		a := math.Sqrt(z*z/σ2X + 1/σ2Y)
		b := μX*z/σ2X + μY/σ2Y
		c := math.Exp(b*b/(2*a*a) - (μ2X/σ2X+μ2Y/σ2Y)/2)
		return b*c/((a*a*a)*math.Sqrt(2*π)*σX*σY)*(2*phi(b/a)-1) + math.Exp(-(μ2X/σ2X+μ2Y/σ2Y)/2)/(a*a*π*σX*σY)
	}
}

// GaussianRatioNoCorrPDFAt returns the value of PDF of Gaussian Ratio distribution of uncorrelated variables, at x. 
func GaussianRatioNoCorrPDFAt(μX, σX, μY, σY, x float64) float64 {
	pdf := GaussianRatioNoCorrPDF(μX, σX, μY, σY)
	return pdf(x)
}

// GaussianRatioPDF returns the value of PDF of Gaussian Ratio distribution of correlated variables, at x. 
func GaussianRatioPDF(μX, σX, μY, σY, ρ float64) func(z float64) float64 {
	return func(z float64) float64 {
		α := ρ * σX / σY
		β := (σX / σY) * math.Sqrt(1-ρ*ρ)
		return β / (π*(z-α)*(z-α) + β*β)

	}
}

// GaussianRatioPDFAt returns the value of PDF of Gaussian Ratio distribution of correlated variables, at x. 
func GaussianRatioPDFAt(μX, σX, μY, σY, ρ, x float64) float64 {
	pdf := GaussianRatioPDF(μX, σX, μY, σY, ρ)
	return pdf(x)
}

// GaussianRatioApproxCDF returns the approximation  of CDF of Gaussian Ratio distribution of correlated variables. 
func GaussianRatioApproxCDF(μX, σX, μY, σY, ρ float64) func(z float64) float64 {
	// Hinkley 1969:636, Eq. 5
	return func(w float64) float64 {
		σ2X := σX * σX
		σ2Y := σY * σY

		a := math.Sqrt(w*w*σ2X - 2*ρ*w/(σX*σY) + (1 / σ2Y))
		// Hinkley 1969:636, Eq. 2

		t1 := μY*w - μX
		t2 := σX * σY * a
		return phi(t1 / t2)
	}
}
