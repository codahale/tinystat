package stat

//Frequentist Confidence Interval for binomial parameter

import "code.google.com/p/probab/dst"

// BinomPConfInt returns a one-sided frequentist Confidence Interval for binomial parameter estimated from a random sample.
// Ref.: Hahn & Meeker (1991).
func BinomPConfInt(n, k int64, alpha float64) (lo, hi float64) {
	// Arguments: 
	// 	n - sample size	
	//	k - observed number of successes (p=n/k)
	//	alpha - lervel of confidence
	//
	// Returns: 
	//	lo	lower confidence limit
	//	hi	upper confidence limit

	if k <= 0 {
		lo = 0.0
	} else {
		lo = 1.0 / (1.0 + float64((n-k+1))*dst.FQtlFor(2*n-2*k+2, 2*k, alpha)/float64(k))
	}

	if k >= n {
		hi = 1.0
	} else {
		hi = 1.0 / (1.0 + float64((n-k))/(float64((k+1))*dst.FQtlFor(2*k+2, 2*n-2*k, alpha)))
	}
	return
}
