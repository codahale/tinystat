// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Gamma distribution, helper functions. 

// Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77
const scalefactor = 1.157920892373162e+77

// If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x
const M_cutoff = Ln2 * maxExp / eps64

// Continued fraction for calculation of
// 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
// auxilary in log1pmx() and lgamma1p()
func logcf(x, i, d, eps float64) float64 {
	c1 := 2 * d
	c2 := i + d
	c4 := c2 + d
	a1 := c2
	b1 := i * (c2 - i*x)
	b2 := d * d * x
	a2 := c4*c2 - b2
	b2 = c4*b1 - i*b2

	for abs(a2*b1-a1*b2) > abs(eps*b1*b2) {
		c3 := c2 * c2 * x
		c2 += d
		c4 += d
		a1 = c4*a2 - c3*a1
		b1 = c4*b2 - c3*b1

		c3 = c1 * c1 * x
		c1 += d
		c4 += d
		a2 = c4*a1 - c3*a2
		b2 = c4*b1 - c3*b2

		if abs(b2) > scalefactor {
			a1 /= scalefactor
			b1 /= scalefactor
			a2 /= scalefactor
			b2 /= scalefactor
		} else if abs(b2) < 1/scalefactor {
			a1 *= scalefactor
			b1 *= scalefactor
			a2 *= scalefactor
			b2 *= scalefactor
		}
	}
	return a2 / b2
}

// Accurate calculation of log(1+x)-x, particularly for small x.
func log1pmx(x float64) float64 {
	const minLog1Value = -0.79149064
	const two = 2.0
	const tol_logcf = 1e-14

	var r, y float64

	if x > 1 || x < minLog1Value {
		return log1p(x) - x
	} else {

		//    -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
		//   log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
		// ---------------------------------------------
		// S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
		//

		r = x / (2 + x)
		y = r * r
		if abs(x) < 1e-2 {
			return r * ((((two/9*y+two/7)*y+two/5)*y+two/3)*y - x)
		}
	}
	return r * (2*y*logcf(y, 3, 2, tol_logcf) - x)
}

// Ln(Abs(Gamma()))
func lgammafn(a float64) float64 {
	return log(abs(Γ(a)))
}

// Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5).
func lgamma1p(a float64) float64 {
	const eulers_const = 0.5772156649015328606065120900824024

	// coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 
	const N = 40
	var coeffs = [40]float64{
		0.3224670334241132182362075833230126e-0, // = (zeta(2)-1)/2
		0.6735230105319809513324605383715000e-1, // = (zeta(3)-1)/3
		0.2058080842778454787900092413529198e-1,
		0.7385551028673985266273097291406834e-2,
		0.2890510330741523285752988298486755e-2,
		0.1192753911703260977113935692828109e-2,
		0.5096695247430424223356548135815582e-3,
		0.2231547584535793797614188036013401e-3,
		0.9945751278180853371459589003190170e-4,
		0.4492623673813314170020750240635786e-4,
		0.2050721277567069155316650397830591e-4,
		0.9439488275268395903987425104415055e-5,
		0.4374866789907487804181793223952411e-5,
		0.2039215753801366236781900709670839e-5,
		0.9551412130407419832857179772951265e-6,
		0.4492469198764566043294290331193655e-6,
		0.2120718480555466586923135901077628e-6,
		0.1004322482396809960872083050053344e-6,
		0.4769810169363980565760193417246730e-7,
		0.2271109460894316491031998116062124e-7,
		0.1083865921489695409107491757968159e-7,
		0.5183475041970046655121248647057669e-8,
		0.2483674543802478317185008663991718e-8,
		0.1192140140586091207442548202774640e-8,
		0.5731367241678862013330194857961011e-9,
		0.2759522885124233145178149692816341e-9,
		0.1330476437424448948149715720858008e-9,
		0.6422964563838100022082448087644648e-10,
		0.3104424774732227276239215783404066e-10,
		0.1502138408075414217093301048780668e-10,
		0.7275974480239079662504549924814047e-11,
		0.3527742476575915083615072228655483e-11,
		0.1711991790559617908601084114443031e-11,
		0.8315385841420284819798357793954418e-12,
		0.4042200525289440065536008957032895e-12,
		0.1966475631096616490411045679010286e-12,
		0.9573630387838555763782200936508615e-13,
		0.4664076026428374224576492565974577e-13,
		0.2273736960065972320633279596737272e-13,
		0.1109139947083452201658320007192334e-13, // (zeta(40+1)-1)/(40+1)
	}

	const c = 0.2273736845824652515226821577978691e-12 // zeta(N+2)-1
	const tol_logcf = 1e-14

	if abs(a) >= 0.5 {
		return lgammafn(a + 1)
	}

	// Abramowitz & Stegun 6.1.33 : 
	// for |x| < 2, <==> log(gamma(1+x)) = 
	// -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
	// where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
	//
	// Here, another convergence acceleration trick is used to compute
	// lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n

	lgam := c * logcf(-a/2, N+2, 1, tol_logcf)
	for i := N - 1; i >= 0; i-- {
		lgam = coeffs[i] - a*lgam
	}
	return (a*lgam-eulers_const)*a - log1pmx(a)
}

//  Compute the log of a sum from logs of terms, i.e.,
//    log (exp (logx) + exp (logy))
// without causing overflows and without throwing away large handfuls
// of accuracy.
func logspace_add(logx, logy float64) float64 {
	return max(logx, logy) + log1p(exp(-abs(logx-logy)))
}

// log1Exp()
func log1Exp(x float64) float64 {
	if x > -Ln2 {
		return log(-expm1(x))
	}
	return log1p(-exp(x))
}

// Compute the log of a difference from logs of terms, i.e.,
//     log (exp (logx) - exp (logy))
// without causing overflows and without throwing away large handfuls of accuracy.
func logspace_sub(logx, logy float64) float64 {
	return logx + log1Exp(logy-logx)
}

// Abramowitz and Stegun 6.5.29 [right]
func pgamma_smallx(x, shape float64) float64 {
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

	f1 := 1 + sum

	if shape > 1 {
		f2 = dpois_raw(shape, x)
		f2 = f2 * exp(x)
	} else {
		f2 = pow(x, shape) / exp(lgamma1p(shape))
	}
	return f1 * f2
}

func pd_upper_series(x, y float64) float64 {
	term := x / y
	sum := term

	for term > sum*eps64 { // was do-while loop, OK???
		y++
		term *= x / y
		sum += term
	}

	// sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
	//	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
	//	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
	//	   ~  x/y +  o(x/y)   {which happens when shape -> Inf}
	return sum
}

// Continued fraction for calculation of
//    scaled upper-tail F_{gamma}
//  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
func pd_lower_cf(y, d float64) float64 {
	var f, of, f0, c2, c3, c4, a1, b1, a2, b2 float64
	f = 0.0
	max_it := 200000

	if y == 0 {
		return 0
	}

	f0 = y / d

	// Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE)
	if abs(y-1) < abs(d)*eps64 { // includes y < d = Inf
		return f0
	}

	if f0 > 1 {
		f0 = 1
	}

	c2 = y
	c4 = d // original (y,d), *not* potentially scaled ones!
	a1 = 0
	b1 = 1
	a2 = y
	b2 = d

	//    while NEEDED_SCALE
	for b2 > scalefactor {
		a1 /= scalefactor
		b1 /= scalefactor
		a2 /= scalefactor
		b2 /= scalefactor
	}

	i := 0
	of = -1 // far away
	for i < max_it {
		i++
		c2--
		c3 = float64(i) * c2
		c4 += 2
		// c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd
		a1 = c4*a2 + c3*a1
		b1 = c4*b2 + c3*b1

		i++
		c2--
		c3 = float64(i) * c2
		c4 += 2
		// c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even
		a2 = c4*a1 + c3*a2
		b2 = c4*b1 + c3*b2

		//	if NEEDED_SCALE
		if b2 > scalefactor {
			a1 /= scalefactor
			b1 /= scalefactor
			a2 /= scalefactor
			b2 /= scalefactor
		}

		if b2 != 0 {
			f = a2 / b2
			// convergence check: relative; "absolute" for very small f
			if abs(f-of) <= eps64*max(f0, abs(f)) {
				return f
			}
			of = f
		}
	}

	panic("NON-convergence in pgamma()'s pd_lower_cf() ")
	return f // should not happen ... 
}

func pd_lower_series(lambda, y float64) float64 {
	var f, term, sum float64
	term = 1
	sum = 0

	for y >= 1 && term > sum*eps64 {
		term *= y / lambda
		sum += term
		y--
	}

	// sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
	//	   =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
	//	   ~  y/lambda + o(y/lambda)

	if y != floor(y) {

		// The series does not converge as the terms start getting
		// bigger (besides flipping sign) for y < -lambda.
		// FIXME: in quite few cases, adding  term*f  has no effect (f too small)
		//	  and is unnecessary e.g. for pgamma(4e12, 121.1)

		f = pd_lower_cf(y, lambda+1-y)
		sum += term * f
	}
	return sum
}

// Compute the following ratio with higher accuracy that would be had
// from doing it directly.
//
//		 dnorm (x, 0, 1, FALSE)
//	   ----------------------------------
//	   pnorm (x, 0, 1, lower_tail, FALSE)
//
// Abramowitz & Stegun 26.2.12
func dpnorm(x, lp float64) float64 {
	// So as not to repeat a pnorm call, we expect
	//
	//	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
	//
	// but use it only in the non-critical case where either x is small
	// or p==exp(lp) is close to 1.

	lower_tail := true
	if x < 0 {
		x = -x
		lower_tail = !lower_tail
	}

	if x > 10 && (!lower_tail) {
		term := 1 / x
		sum := term
		x2 := x * x
		i := 1
		term = 1e32 // just to enter the while loop
		for abs(term) > eps64*sum {
			term *= float64(-i) / x2
			sum += term
			i += 2
		}

		return 1 / sum
	} else {
		//	 d := dnorm (x, 0, 1, FALSE)
		d := ZPDFAt(x)
		return d / exp(lp)
	}
	return NaN // should not happen
}

// Asymptotic expansion to calculate the probability that Poisson variate
// has value <= x.
// Various assertions about this are made (without proof) at
// http://members.aol.com/iandjmsmith/PoissonApprox.htm
func ppois_asymp(lambda, x float64, log_p bool) float64 {
	var coefs_a = [8]float64{
		-1e9, // placeholder used for 1-indexing
		2 / 3.0,
		-4 / 135.0,
		8 / 2835.0,
		16 / 8505.0,
		-8992 / 12629925.0,
		-334144 / 492567075.0,
		698752 / 1477701225.0,
	}

	var coefs_b = [8]float64{
		-1e9, // placeholder
		1 / 12.0,
		1 / 288.0,
		-139 / 51840.0,
		-571 / 2488320.0,
		163879 / 209018880.0,
		5246819 / 75246796800.0,
		-534703531 / 902961561600.0,
	}

	var (
		elfb, elfb_term                               float64
		res12, res1_term, res1_ig, res2_term, res2_ig float64
		dfm, pt_, s2pt, f, np                         float64
	)

	lower_tail := true
	dfm = lambda - x

	// If lambda is large, the distribution is highly concentrated
	//  about lambda.  So representation error in x or lambda can lead
	//   to arbitrarily large values of pt_ and hence divergence of the
	//   coefficients of this approximation.

	pt_ = -log1pmx(dfm / x)
	s2pt = sqrt(2 * x * pt_)
	if dfm < 0 {
		s2pt = -s2pt
	}
	res12 = 0
	res1_term = sqrt(x)
	res1_ig = res1_term
	res2_term = s2pt
	res2_ig = res2_term
	for i := 1; i < 8; i++ {
		res12 += res1_ig * coefs_a[i]
		res12 += res2_ig * coefs_b[i]
		res1_term *= pt_ / float64(i)
		res2_term *= 2 * pt_ / (2*float64(i) + 1)
		res1_ig = res1_ig/x + res1_term
		res2_ig = res2_ig/x + res2_term
	}

	elfb = x
	elfb_term = 1
	for i := 1; i < 8; i++ {
		elfb += elfb_term * coefs_b[i]
		elfb_term /= x
	}

	if !lower_tail { // should never happen
		elfb = -elfb
	}

	f = res12 / elfb
	np = ZCDFAt(s2pt)

	if log_p {
		//	n_d_over_p := dpnorm(s2pt, !lower_tail, np)
		n_d_over_p := dpnorm(s2pt, np)
		return np + log1p(f*n_d_over_p)
	} else {
		nd := ZPDFAt(s2pt)

		return np + f*nd
	}
	return NaN // should not happen
}

func dpois_wrap(x_plus_1, lambda float64) float64 {
	if isInf(lambda, 0) {
		return 0
	}

	if x_plus_1 > 1 {
		//		return dpois_raw(x_plus_1-1, lambda, give_log)
		return PoissonPMFAt(lambda, int64(x_plus_1-1))
	}

	if lambda > abs(x_plus_1-1)*M_cutoff {
		return exp(-lambda - lgammafn(x_plus_1))
	}
	d := PoissonPMFAt(lambda, int64(x_plus_1))
	return d * (x_plus_1 / lambda)
}

func dpois_raw(x, lambda float64) float64 {
	// x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
	//        lambda >= 0

	if lambda == 0 {
		if x == 0 {
			return 1
		} else {
			return 0
		}
	}

	if isInf(lambda, 0) {
		return 0
	}

	if x < 0 {
		return 0
	}

	if x <= lambda*min64 {
		return exp(-lambda)
	}

	if lambda < x*min64 {
		return exp(-lambda + x*log(lambda) - lgammafn(x+1))
	}

	return exp(-stirlerr(x)-bd0(x, lambda)) / sqrt((π+π)*x)
}

func pgamma_raw(x, shape float64) float64 {
	// Here, assume that  (x,shape) are not NA  &  shape > 0 . 

	var res, sum float64

	if x <= 0 {
		return 0
	}
	if x >= posInf {
		return 1
	}

	if x < 1 {
		res = pgamma_smallx(x, shape)
	} else if x <= shape-1 && x < 0.8*(shape+50) {
		// incl. large shape compared to x
		sum = pd_upper_series(x, shape) /* = x/shape + o(x/shape) */
		d := dpois_wrap(shape, x)
		res = sum * d
	} else if shape-1 < x && shape < 0.8*(x+50) {
		// incl. large x compared to shape
		d := dpois_wrap(shape, x)
		if shape < 1 {
			if x*eps64 > 1-shape {
				//				sum = R_D__1
				sum = 1
			} else {
				f := pd_lower_cf(shape, x-(shape-1)) * x / shape
				// = [shape/(x - shape+1) + o(shape/(x-shape+1))] * x/shape = 1 + o(1)
				sum = f
			}
		} else {
			sum = pd_lower_series(x, shape-1) // = (shape-1)/x + o((shape-1)/x)
			sum = 1 + sum
		}
		res = 1 - d*sum
	} else { // x >= 1 and x fairly near shape.
		res = ppois_asymp(shape-1, x, false)
	}

	// We lose a fair amount of accuracy to underflow in the cases
	// where the final result is very close to min64.	
	//  In those cases, simply redo via logarithm.
	if res < min64/eps64 {
		return exp(pgamma_raw_ln(shape, x))
	}
	return res
}

func qchisq_appr(p, nu, g float64, lower_tail, log_p bool, tol float64) float64 {
	// g  = log Gamma(nu/2)

	const (
		C7  = 4.67
		C8  = 6.66
		C9  = 6.73
		C10 = 13.32
	)

	var (
		alpha, a, c, ch, lgam1pa, p1 float64
		p2, q, t, x                  float64
	)

	if isNaN(p) || isNaN(nu) {
		return p + nu
	}

	if (log_p && p > 0) || (!log_p && (p < 0 || p > 1)) {
		return NaN
	}

	if nu <= 0 {
		return NaN
	}

	alpha = 0.5 * nu // = [pq]gamma() shape
	c = alpha - 1
	p1 = log(p)

	if nu < (-1.24)*(p1) { // for small chi-squared

		// log(alpha) + g = log(alpha) + log(gamma(alpha)) =
		//        = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
		//  catastrophic cancellation when alpha << 1

		//	 lgam1pa = (alpha < 0.5) ? lgamma1p(alpha) : (log(alpha) + g)
		if alpha < 0.5 {
			lgam1pa = lgamma1p(alpha)
		} else {
			lgam1pa = log(alpha) + g
		}
		ch = exp((lgam1pa+p1)/alpha + Ln2)
	} else if nu > 0.32 {

		//  using Wilson and Hilferty estimate
		//	x = qnorm(p, 0, 1, lower_tail, log_p)
		x = ZQtlFor(p)
		p1 = 2. / (9 * nu)
		ch = nu * pow(x*sqrt(p1)+1-p1, 3)

		// approximation for p tending to 1
		if ch > 2.2*nu+6 {
			//	    ch = -2*(R_DT_Clog(p) - c*log(0.5*ch) + g)
			ch = -2 * (log1p(-p) - c*log(0.5*ch) + g)
		}
	} else { // "small nu" : 1.24*(-log(p)) <= nu <= 0.32

		ch = 0.4
		//	a = R_DT_Clog(p) + g + c*M_LN2
		a = log1p(-p) + g + c*Ln2

		q = 1 // to enter the while loop
		for abs(q-ch) > tol*abs(ch) {
			q = ch
			p1 = 1. / (1 + ch*(C7+ch))
			p2 = ch * (C9 + ch*(C8+ch))
			t = -0.5 + (C7+2*ch)*p1 - (C9+ch*(C10+3*ch))/p2
			ch -= (1 - exp(a+0.5*ch)*p2*p1) / t
		} //while
	}
	return ch
}
