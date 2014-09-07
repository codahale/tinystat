// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Gamma distribution, helper functions, log versions. 

func pgamma_raw_ln(x, shape float64) float64 {
	// Here, assume that  (x, shape) are not NA  &  shape > 0 .
	var res, sum float64

	if x <= 0 {
		return negInf
	}

	if x >= posInf {
		return 0
	}

	if x < 1 {
		res = pgamma_smallx_ln(x, shape)
	} else if x <= shape-1 && x < 0.8*(shape+50) {
		// incl. large shape compared to x
		sum = log(pd_upper_series(x, shape)) // = x/shape + o(x/shape)
		d := dpois_wrap_ln(shape, x)

		//	    res = log_p ? sum + d : sum * d
		res = sum + d

	} else if shape-1 < x && shape < 0.8*(x+50) {
		// incl. large x compared to shape

		d := dpois_wrap_ln(shape, x)
		if shape < 1 {
			if x*eps64 > 1-shape {
				//				sum = R_D__1
				sum = 0
			} else {
				f := pd_lower_cf(shape, x-(shape-1)) * x / shape
				// = [shape/(x - shape+1) + o(shape/(x-shape+1))] * x/shape = 1 + o(1)
				//		sum = log_p ? log (f) : f

				sum = log(f)
			}
		} else {
			sum = pd_lower_series(x, shape-1) // = (shape-1)/x + o((shape-1)/x)

			//	    sum = log_p ? log1p (sum) : 1 + sum
			sum = log1p(sum)

		}
		//	    res = log_p	? R_Log1_Exp (d + sum)	: 1 - d * sum
		res = log1Exp(d + sum)

	} else { /* x >= 1 and x fairly near shape. */
		res = ppois_asymp(shape-1, x, true)
	}
	return res
}

// Abramowitz and Stegun 6.5.29 [right]
func pgamma_smallx_ln(x, shape float64) float64 {
	var term, f2 float64
	sum := 0.0
	c := shape
	n := 0.0

	// Relative to 6.5.29 all terms have been multiplied by shape
	// and the first, thus being 1, is omitted.
	term = 1e32 // just to enter the while loop
	for abs(term) > eps64*abs(sum) {
		n++
		c *= -x / n
		term = c / (shape + n)
		sum += term
	}

	f1 := log1p(sum)

	if shape > 1 {
		f2 = dpois_raw_ln(shape, x)
		f2 = f2 + x
	} else {
		f2 = shape*log(x) - lgamma1p(shape)
	}
	return f1 + f2
}

func dpois_wrap_ln(x_plus_1, lambda float64) float64 {
	if isInf(lambda, 0) {
		return negInf
	}

	if x_plus_1 > 1 {
		return dpois_raw_ln(x_plus_1-1, lambda)

	}

	if lambda > abs(x_plus_1-1)*M_cutoff {
		return -lambda - lgammafn(x_plus_1)
	}

	d := dpois_raw_ln(x_plus_1, lambda)
	return d + log(x_plus_1/lambda)
}

func dpois_raw_ln(x, lambda float64) float64 {
	// x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
	//        lambda >= 0

	if lambda == 0 {
		if x == 0 {
			return 0
		} else {
			return negInf
		}
	}

	if isInf(lambda, 0) {
		return negInf
	}

	if x < 0 {
		return negInf
	}

	if x <= lambda*min64 {
		return -lambda
	}

	if lambda < x*min64 {
		return -lambda + x*log(lambda) - lgammafn(x+1)
	}

	return -0.5*log((π+π)*x) + (-stirlerr(x) - bd0(x, lambda))
}
