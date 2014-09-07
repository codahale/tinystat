// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Student's t  distribution. 
// A family of continuous probability distributions that arises when estimating the mean of a normally distributed population in situations where the sample size is small and population standard deviation is unknown.
//
// Parameters: 
// ν > 0	degrees of freedom (real)
//
// Support: 
// x ∈ (-∞, +∞) (real)

// StudentsTPDF returns the PDF of the Student's t distribution. 
func StudentsTPDF(ν float64) func(x float64) float64 {
	normalization := Γ((ν+1)/2) / (sqrt(ν*π) * Γ(ν/2))
	return func(x float64) float64 {
		return normalization * pow(1+x*x/ν, -(ν+1)/2)
	}
}

// StudentsTLnPDF returns the natural logarithm of the PDF of the Student's t distribution. 
func StudentsTLnPDF(ν float64) func(x float64) float64 {
	normalization := LnΓ((ν+1)/2) - log(sqrt(ν*π)) - LnΓ(ν/2)
	return func(x float64) float64 {
		return normalization + log(1+x*x/ν)*-(ν+1)/2
	}
}

// StudentsTCDF returns the CDF of the Student's t distribution. 
func StudentsTCDF(ν float64) func(x float64) float64 {
	return func(x float64) float64 {
		var p float64
		if ν <= 0 {
			return NaN
		}

		nx := 1 + (x/ν)*x
		if nx > 1e100 { /* <==>  x*x > 1e100 * ν  */
			/* Danger of underflow. So use Abramowitz & Stegun 26.5.4
			   pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
			   with z = 1/nx,  a = ν/2,  b= 1/2 :
			*/
			p = -0.5*ν*(2*log(abs(x))-log(ν)) - logB(0.5*ν, 0.5) - log(0.5*ν)
			p = exp(p)
		} else {
			if ν > x*x {
				α := 0.5
				β := ν / 2
				pbeta := BetaCDF(α, β)
				p = 1 - pbeta(x*x/(ν+x*x))
			} else {
				α := ν / 2
				β := 0.5
				pbeta := BetaCDF(α, β)
				p = pbeta(1 / nx)
			}
		}

		p /= 2
		if x > 0 {
			p = 1 - p
		}
		return p
	}
}

// StudentsTCDFAt returns the value of CDF of the Student's t distribution, at x. 
func StudentsTCDFAt(ν, x float64) float64 {
	cdf := StudentsTCDF(ν)
	return cdf(x)
}

// StudentsTQtl returns the inverse of the CDF (quantile) of the Student's t distribution. 
func StudentsTQtl(ν float64) func(p float64) float64 {
	// Hill, G.W (1970) "Algorithm 396: Student's t-quantiles"
	// CACM 13(10), 619-620.
	// Using expm1() takes care of  Lozy (1979) "Remark on Algo.", TOMS
	// Applies 2-term Taylor expansion as in Hill, G.W (1981) "Remark on Algo.396", ACM TOMS 7, 250-1
	// Improved formula for decision when 1 < df < 2
	return func(p float64) float64 {
		const eps = 1.e-12
		var q float64
		neg := false
		pok := false

		if ν <= 0 || p < 0 || p > 1 {
			return NaN
		}

		/*
			    if (ν < 1) { // based on qnt
				const static double accu = 1e-13;
				const static double Eps = 1e-11; // must be > accu

				double ux, lx, nx, pp;

				int iter = 0;

				p = RDTqIv(p);

				// Invert pt(.) :
				// 1. finding an upper and lower bound
				if(p > 1 - min64) return MLPOSINF;
				pp = fmin2(1 - min64, p// (1 + Eps));
				for(ux = 1.; ux < DBLMAX && pt(ux, ν, TRUE, FALSE) < pp; ux//= 2);
				pp = p// (1 - Eps);
				for(lx =-1.; lx > -DBLMAX && pt(lx, ν, TRUE, FALSE) > pp; lx//= 2);

				// 2. interval (lx,ux)  halving
				   regula falsi failed on qt(0.1, 0.1)

				do {
				    nx = 0.5// (lx + ux);
				    if (pt(nx, ν, TRUE, FALSE) > p) ux = nx; else lx = nx;
				} while ((ux - lx) / abs(nx) > accu && ++iter < 1000);

				if(iter >= 1000) MLERROR(MEPRECISION, "qt");

				return 0.5// (lx + ux);
			    }
		*/

		if ν > 1e20 {
			q = ZQtlFor(p)
		} else {
			if p < 0.5 {
				neg = true
			}
			if neg {
				p = 2 * p
			} else {
				p = 2 * (0.5 - p + 0.5)
			}

			if abs(ν-2) < eps { // df ~= 2
				if p > min64 {
					if 3*p < min64 { // p ~= 0
						q = 1 / sqrt(p)
					} else if p > 0.9 { // p ~= 1
						q = (1 - p) * sqrt(2/(p*(2-p)))
					} else { // eps/3 <= p <= 0.9
						q = sqrt(2/(p*(2-p)) - 2)
					}
				} else { // p << 1, q = 1/sqrt(p) = ...
					return posInf

				}
			} else if ν < 1+eps { // df ~= 1  (df < 1 excluded above): Cauchy
				if p > 0 {
					q = 1 / tan(p*π/2) // == - tan((p+1) * π/2) -- suffers for p ~= 0

				} else { // p = 0, but maybe = 2*exp(p) !
					return posInf
				}
			} else { //-- usual case;  including, e.g.,  df = 1.1
				x := 0.0
				y := 0.0
				a := 1 / (ν - 0.5)
				b := 48 / (a * a)
				c := ((20700*a/b-98)*a-16)*a + 96.36
				d := ((94.5/(b+c)-3)/b + 1) * sqrt(a*π/2) * ν

				y = pow(d*p, 2/ν)
				if y >= min64 {
					pok = true
				}
				if (ν < 2.1 && p > 0.5) || y > 0.05+a { // p > p0(df)
					// Asymptotic inverse expansion about normal
					x = NormalQtlFor(0, 1, 0.5*p)

					y = x * x
					if ν < 5 {
						c += 0.3 * (ν - 4.5) * (x + 0.6)
					}
					c = (((0.05*d*x-5)*x-7)*x-2)*x + b + c
					y = (((((0.4*y+6.3)*y+36)*y+94.5)/c-y-3)/b + 1) * x
					y = expm1(a * y * y)
					q = sqrt(ν * y)
				} else { // re-use 'y' from above

					if !pok && x < -0.5*log(min64) { // 0.5* log(min64)
						// y above might have underflown
						q = sqrt(ν) * exp(-x)
					} else {
						y = ((1/(((ν+6)/(ν*y)-0.089*d-0.822)*(ν+2)*3)+0.5/(ν+4))*y-1)*(ν+1)/(ν+2) + 1/y
						q = sqrt(ν * y)
					}
				}

				// Now apply 2-term Taylor expansion improvement (1-term = Newton): as by Hill (1981) [ref.above]
				dt := StudentsTPDF(ν)
				pt := StudentsTCDF(ν)
				for it := 0; it < 10; it++ {
					y = dt(q)
					if y <= 0 {
						break
					}
					x = (1 - pt(q) - p/2) / y
					if abs(x) > 1e-14*abs(q) {
						// Newton (=Taylor 1 term):
						//  q += x 
						// Taylor 2-term : 
						q += x * (1. + x*q*(ν+1)/(2*(q*q+ν)))
					}
				}
			}
		}
		if neg {
			q = -q
		}
		return q
	}
}

// StudentsTQtlFor returns the inverse of the CDF (quantile) of the Student's t distribution, for given probability.
func StudentsTQtlFor(ν, p float64) float64 {
	qtl := StudentsTQtl(ν)
	return qtl(p)
}

// StudentsTNext returns random number drawn from the Student's t distribution. 
func StudentsTNext(ν float64) float64 {
	return NormalNext(0, 1) * sqrt(ν/GammaNext(ν/2, 2))
}

// StudentsT returns the random number generator with  Student's t distribution. 
func StudentsT(ν float64) func() float64 {
	return func() float64 {
		return StudentsTNext(ν)
	}
}

// StudentsTMean returns the mean of the StudentsT Type I distribution. 
func StudentsTMean(ν float64) float64 {
	if ν <= 1 {
		return NaN
	}
	return 0
}

// StudentsTMode returns the mode of the StudentsT Type I distribution. 
func StudentsTMode(ν float64) float64 {
	return 0
}

// StudentsTMedian returns the median of the StudentsT Type I distribution. 
func StudentsTMedian(ν float64) float64 {
	return 0
}

// StudentsTVar returns the variance of the StudentsT Type I distribution. 
func StudentsTVar(ν float64) float64 {
	if ν >= 1 {
		return NaN
	}
	if ν > 2 {
		return ν / (ν - 2)
	}
	return posInf
}

// StudentsTStd returns the standard deviation of the StudentsT Type I distribution. 
func StudentsTStd(ν float64) float64 {
	if ν >= 1 {
		return NaN
	}
	if ν > 2 {
		return sqrt(ν / (ν - 2))
	}
	return posInf
}

// StudentsTSkew returns the skewness of the StudentsT Type I distribution. 
func StudentsTSkew(ν float64) float64 {
	if ν <= 3 {
		return NaN
	}
	return 0
}

// StudentsTExKurt returns the excess kurtosis of the StudentsT Type I distribution. 
func StudentsTExKurt(ν float64) float64 {
	if ν <= 2 {
		return NaN
	}
	if ν <= 4 {
		return posInf
	}
	return 6 / (ν - 4)
}
