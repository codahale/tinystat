// Copyright 2012 The Probab Authors. All rights reserved. See the LICENSE file.

package dst

// Random variates from the Poisson distribution.
// rewritten from  Mathlib : A C Library of Special Functions
// Ahrens and Dieter (1982).

import (
	"math/rand"
)

// PoissonNext returns random number drawn from the Poisson distribution. 
func PoissonNext(λ float64) int64 {
	const (
		a0     = -0.5
		a1     = 0.3333333
		a2     = -0.2500068
		a3     = 0.2000118
		a4     = -0.1661269
		a5     = 0.1421878
		a6     = -0.1384794
		a7     = 0.1250060
		one_7  = 0.1428571428571428571
		one_12 = 0.0833333333333333333
		one_24 = 0.0416666666666666667
	)

	var (
		b1, b2, c, c0, c1, c2, c3       float64
		p0, p, q, s, d, omega           float64
		big_l                           float64 // integer "w/o overflow"
		del, fx, fy, g, px, py, t, v, x float64
		m, l, k                         int
	)

	// Factorial Table (0:9)! 
	fact := []float64{1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.}
	pp := make([]float64, 36)
	difmuk := 0.0
	E := 0.0
	fk := 0.0
	u := 0.0
	pois := -1.0

	if isInf(λ, 1) || λ < 0.0 {
		//		return NaN
		panic("bad lambda")
	}

	if λ == 0.0 {
		return 0
	}

	k = 0
	kflag := false
	big_mu := false
	new_big_mu := false
	stepF := false

	muprev := 0.0
	muprev2 := 0.0

	if λ >= 10 {
		big_mu = true
		new_big_mu = false
	} else {
		big_mu = false
	}

	//	if !(big_mu && mu == muprev) { // maybe compute new pars

	if big_mu {
		new_big_mu = true

		// Case A. (recalculation of s,d,l	because λ has changed):
		// The poisson probabilities pk exceed the discrete normal
		// probabilities fk whenever k >= m(λ).

		muprev = λ
		s = sqrt(λ)
		d = 6. * λ * λ
		big_l = floor(λ - 1.1484)
		// = an upper bound to m(λ) for all λ >= 10.
	} else { // Small λ ( < 10) -- not using normal approx.

		// Case B. (start new table and calculate p0 if necessary)

		if λ != muprev {
			muprev = λ
			m = imax2(1, int(λ))
			l = 0 // pp[] is already ok up to pp[l]
			p = exp(-λ)
			p0 = p
			q = p
		}

		for {
			// Step U. uniform sample for inversion method
			u := rand.Float64()
			if u <= p0 {
				return 0
			}

			// Step T. table comparison until the end pp[l] of the
			//   pp-table of cumulative poisson probabilities
			//   (0.458 > ~= pp[9](= 0.45792971447) for λ=10 )
			if l != 0 {
				kk := 1
				if u > 0.458 {
					kk = imin2(l, m)
				}
				for k = kk; k <= l; k++ {
					if u <= pp[k] {
						return int64(k)
					}
				}
				if l == 35 { // u > pp[35]
					continue
				}
			}
			// Step C. creation of new poisson
			//   probabilities p[l..] and their cumulatives q =: pp[k] 
			l++
			for k = l; k <= 35; k++ {
				p *= λ / float64(k)
				q += p
				pp[k] = q
				if u <= q {
					l = k
					return int64(k)
				}
			}
			l = 35
		} // end(for)
	} // λ < 10
	//	} // end initialize persistent vars

	// Only if λ >= 10

	// Step N. normal sample
	g = λ + s*rand.NormFloat64() // norm_rand() ~ N(0,1), standard normal

	if g >= 0. {
		pois = floor(g)
		// Step I. immediate acceptance if pois is large enough
		if pois >= big_l {
			return int64(pois)
		}
		// Step S. squeeze acceptance
		fk = pois
		difmuk = λ - fk
		u = rand.Float64() // ~ U(0,1) - sample
		if d*u >= difmuk*difmuk*difmuk {
			return int64(pois)
		}
	}

	// Step P. preparations for steps Q and H.
	//   (recalculations of parameters if necessary)

	if new_big_mu || λ != muprev2 {
		// Careful! muprev2 is not always == muprev
		//  because one might have exited in step I or S
		muprev2 = λ
		omega = M_1_SQRT_2PI / s

		// The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
		// approximations to the discrete normal probabilities fk.

		b1 = one_24 / λ
		b2 = 0.3 * b1 * b1
		c3 = one_7 * b1 * b2
		c2 = b2 - 15.*c3
		c1 = b1 - 6.*b2 + 45.*c3
		c0 = 1. - b1 + 3.*b2 - 15.*c3
		c = 0.1069 / λ // guarantees majorization by the 'hat'-function.
	}
	if g >= 0. {
		// 'Subroutine' F is called (kflag=0 for correct return)
		kflag = false
		//		goto Step_F
		stepF = true
	}

	for {
		if !stepF {
			// Step E. Exponential Sample

			E = rand.ExpFloat64() // ~ Exp(1) (standard exponential)

			//  sample t from the laplace 'hat'
			//    (if t <= -0.6744 then pk < fk for all λ >= 10.)
			u = 2*rand.Float64() - 1.
			t = 1.8 + fsign(E, u)
		}
		if t > -0.6744 || stepF {
			if !stepF {
				pois = floor(λ + s*t)
				fk = pois
				difmuk = λ - fk

				// 'subroutine' F is called (kflag=1 for correct return)
				kflag = true
			}
			// Step_F:  calculation of px,py,fx,fy.

			if pois < 10 { // use factorials from table fact[] 
				px = -λ
				py = pow(λ, pois) / fact[int(pois)]
			} else {
				// Case pois >= 10 uses polynomial approximation
				//  a0-a7 for accuracy when advisable
				del = one_12 / fk
				del = del * (1. - 4.8*del*del)
				v = difmuk / fk
				if abs(v) <= 0.25 {
					px = fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
				} else { // |v| > 1/4 
					px = fk*log(1.+v) - difmuk - del
				}
				py = M_1_SQRT_2PI / sqrt(fk)
			}
			x = (0.5 - difmuk) / s
			x *= x
			fx = -0.5 * x
			fy = omega * (((c3*x+c2)*x+c1)*x + c0)
			if kflag {
				// Step H. Hat acceptance (E is fored on rejection)
				if c*abs(u) <= py*exp(px+E)-fy*exp(fx+E) {
					//					break L1
					return int64(pois)
				}
			} else {
				// Step Q. Quotient acceptance (rare case)
				if fy-u*fy <= py*exp(px-fx) {
					//					break L1
					return int64(pois)
				}
			}
			stepF = false
		} // t > -.67.. 
	}

	return int64(pois)
}
