// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Gamma distribution. 
// Parameters: 
// k > 0.0		shape parameter, 
// θ (Theta) > 0.0	scale parameter. 
// Alternatively, shape parameter α = k and an inverse scale parameter β = 1⁄θ, is called a rate parameter.
// If k is an integer, then the distribution represents an Erlang distribution; i.e., the sum of k  independent exponentially-distributed random variables, each of which has a mean of θ (which is equivalent to a rate parameter of 1/θ). Equivalently, if α is an integer, then the distribution again represents an Erlang distribution, i.e. the sum of α independent exponentially-distributed random variables, each of which has a mean of 1/β (which is equivalent to a rate parameter of β).
// Support: 
// x ∈ (0, ∞)

// GammaQtl returns the inverse of the CDF (quantile) of the Gamma distribution. 
func GammaQtl(alpha, scale float64) func(p float64) float64 {
	/*	This function is based on the Applied Statistics
	 *	Algorithm AS 91 ("ppchi2") and via pgamma(.) AS 239.
	 *
	 *	R core improvements:
	 *	o  lower_tail, log_p
	 *      o  non-trivial result for p outside [0.000002, 0.999998]
	 *	o  p ~ 1 no longer gives +Inf; final Newton step(s)
	 *
	 *  REFERENCES
	 *
	 *	Best, D. J. and D. E. Roberts (1975).
	 *	Percentage Points of the Chi-Squared Distribution.
	 *	Applied Statistics 24, page 385.  
	 */

	return func(p float64) float64 {

		lower_tail := true // to be removed
		log_p := false

		const (
			EPS1   = 1e-2
			EPS2   = 5e-7                /* final precision of AS 91 */
			EPS_N  = 1e-15               /* precision of Newton step / iterations */
			LN_EPS = -36.043653389117156 /* = log(.Machine$float64.eps) iff IEEE_754 */

			MAXIT = 1000

			pMIN = 1e-100      /* was 0.000002 = 2e-6 */
			pMAX = (1 - 1e-14) /* was (1-1e-12) and 0.999998 = 1 - 2e-6 */

			i420  = 1. / 420.
			i2520 = 1. / 2520.
			i5040 = 1. / 5040
		)
		var (
			p_, a, b, c, g, ch, ch0, p1         float64
			p2, q, s1, s2, s3, s4, s5, s6, t, x float64
			i, max_it_Newton                    int
		)
		/* test arguments and initialise */

		if isNaN(p) || isNaN(alpha) || isNaN(scale) {
			return p + alpha + scale
		}
		//    R_Q_P01_boundaries(p, 0., ML_POSINF)
		if p < 0 || p > 1 {
			return NaN
		}
		if p == 0 {

			return 0
		}
		if p == 1 {
			return posInf
		}

		if alpha < 0 || scale <= 0 {
			return NaN
		}

		if alpha == 0 { // all mass at 0
			return 0
		}

		max_it_Newton = 1
		if alpha < 1e-10 {
			max_it_Newton = 7 // may still be increased below
		}

		//    p_ = R_DT_qIv(p)// lower_tail prob (in any case)
		p_ = p

		g = lgammafn(alpha) // log Gamma(v/2) 

		// Phase I : Starting Approximation
		ch = qchisq_appr(p, 2*alpha, g, lower_tail, log_p, EPS1)
		if isInf(ch, 0) {
			/* forget about all iterations! */
			max_it_Newton = 0
			goto END // TO BE IMPROVED
		}
		if ch < EPS2 { /* Corrected according to AS 91; MM, May 25, 1999 */
			max_it_Newton = 20
			goto END /* and do Newton steps */ // TO BE IMPROVED
		}

		/* FIXME: This (cutoff to {0, +Inf}) is far from optimal
		 * -----  when log_p or !lower_tail, but NOT doing it can be even worse */
		if p_ > pMAX || p_ < pMIN {
			/* did return ML_POSINF or 0.;	much better: */
			max_it_Newton = 20
			goto END /* and do Newton steps */ // TO BE IMPROVED
		}

		/*----- Phase II: Iteration
		 *	Call pgamma() [AS 239]	and calculate seven term taylor series
		 */
		c = alpha - 1
		s6 = (120 + c*(346+127*c)) * i5040 /* used below, is "const" */

		ch0 = ch /* save initial approx. */
		for i = 1; i <= MAXIT; i++ {
			q = ch
			p1 = 0.5 * ch
			//	p2 = p_ - pgamma_raw(p1, alpha, /*lower_tail*/TRUE, /*log_p*/FALSE)
			p2 = p_ - pgamma_raw(p1, alpha)

			if isInf(p2, 0) || ch <= 0 {
				ch = ch0
				max_it_Newton = 27
				goto END
			}

			t = p2 * exp(alpha*Ln2+g+p1-c*log(ch))
			b = t / ch
			a = 0.5*t - b*c
			s1 = (210 + a*(140+a*(105+a*(84+a*(70+60*a))))) * i420
			s2 = (420 + a*(735+a*(966+a*(1141+1278*a)))) * i2520
			s3 = (210 + a*(462+a*(707+932*a))) * i2520
			s4 = (252 + a*(672+1182*a) + c*(294+a*(889+1740*a))) * i5040
			s5 = (84 + 2264*a + c*(1175+606*a)) * i2520

			ch += t * (1 + 0.5*t*s1 - b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))))
			if abs(q-ch) < EPS2*ch {
				goto END
			}
			if abs(q-ch) > 0.1*ch { /* diverging? -- also forces ch > 0 */
				if ch < q {
					ch = 0.9 * q
				} else {
					ch = 1.1 * q
				}
			}
		}

		/* no convergence in MAXIT iterations -- but we add Newton now... */

	END:
		/* PR# 2214 :	 From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50
		 * With a final Newton step, float64 accuracy, e.g. for (p= 7e-4; nu= 0.9)
		 *
		 * Improved (MM): - only if rel.Err > EPS_N (= 1e-15)
		 *		    - also for lower_tail = FALSE	 or log_p = TRUE
		 * 		    - optionally *iterate* Newton
		 */
		x = 0.5 * scale * ch

		if max_it_Newton != 0 {
			/* always use log scale */
			//	if (!log_p) {
			p = log(p)
			log_p = true
			//	}
			if x == 0 {
				_1_p := 1. + 1e-7
				//				_1_m := 1. - 1e-7
				x = min64
				//	    p_ = pgamma(x, alpha, scale, lower_tail, log_p)
				p_ = GammaLnCDFAt(alpha, scale, x)

				//	    if(( lower_tail && p_ > p * _1_p) || (!lower_tail && p_ < p * _1_m))
				if lower_tail && p_ > p*_1_p {
					return 0
				}
			} else { // continue, using x = min64 instead of  0
				//	    p_ = pgamma(x, alpha, scale, lower_tail, log_p)
				p_ = GammaLnCDFAt(alpha, scale, x)
			}
			if p_ == negInf {
				return 0 /* PR#14710 */
			}
			for i = 1; i <= max_it_Newton; i++ {
				p1 = p_ - p
				if abs(p1) < abs(EPS_N*p) {
					break
				}
				/* else */
				// g = dgamma(x, alpha, scale, log_p))
				g = GammaLnPDFAt(alpha, scale, x)
				//	    if g  == R_D__0 
				if g == negInf {
					break
				}
				/* else :
				 * delta x = f(x)/f'(x)
				 * if(log_p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
				 * ==> f(x)/f'(x) = f*P / P' = f*exp(p_) / P' (since p_ = log P(x))
				 */
				//	    t = log_p ? p1*exp(p_ - g) : p1/g ;

				if log_p { // ALWAYS, 
					t = p1 * exp(p_-g) // = "delta x"
				} else { // TO BE REMOVED
					t = p1 / g
				}

				//	    t = lower_tail ? x - t : x + t
				t = x - t

				//	    p_ = pgamma (t, alpha, scale, lower_tail, log_p)
				p_ = GammaLnCDFAt(alpha, scale, x)
				if abs(p_-p) > abs(p1) || (i > 1 && abs(p_-p) == abs(p1)) { // <- against flip-flop
					// no improvement
					break
				} // else : 
				//ifdef Harmful_notably_if_max_it_Newton_is_1
				// control step length: this could have started at the initial approximation 

				if t > 1.1*x {
					t = 1.1 * x
				} else if t < 0.9*x {
					t = 0.9 * x
				}
				//endif
				x = t
			}
		}
		return x
	}
}

// GammaQtlFor returns the inverse of the CDF (quantile) of the Gamma distribution, for given probability.
func GammaQtlFor(k, θ, p float64) float64 {
	cdf := GammaQtl(k, θ)
	return cdf(p)
}
