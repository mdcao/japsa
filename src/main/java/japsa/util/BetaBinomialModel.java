/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * 05/08/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.util;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathUtils;

/**
 * Implementation of a beta-binomial model
 * 
 * @see <a href="http://en.wikipedia.org/wiki/Beta-binomial_distribution">Beta-binomial distribution (Wikipedia)</a>
 * @author minhduc
 *
 */
public class BetaBinomialModel {


	static int psuedoCount = 1; 

	/*OK, let do the maths
	 * count_flank = #read in the flanks
	 * count_repeat = #reads in repeats
	 * count_all = count_flank + count_repeat
	 * 
	 * assume count_flank follows Bin(count_all,p), so p would follow the 
	 * beta distribution Beta(alpha,beta) which alpha= count_flank, beta=count_repeat
	 * 
	 * 
	 * 
	 */
	/**
	 * Compute the ratio distribution of ps/pr where ps and pr are the params of
	 * the binomial distributionsm from countR vs totR/countS vs totS 
	 * 
	 * Use a normal distribution for now. Use sampling to compute a nomal distribution
	 * @param countS
	 * @param totS
	 * @param countR
	 * @param totR
	 * @param numSamples
	 * @return
	 */
	public static NormalDistribution ratioDistribution(double countR, double totR, double countS, double totS, int numSamples){		

		BetaDistribution betaR = 
			new BetaDistribution(countR + psuedoCount, totR - countR + psuedoCount);

		BetaDistribution betaS = 
			new BetaDistribution(countS + psuedoCount, totS - countS + psuedoCount);
		
		//betaR models the  

		double sum = 0, sq = 0;
		//int count = 0;
		for (int i = 0; i < numSamples; i++){
			double r = betaR.sample();

			if (r == 0 ){//wont happen{
				i --;continue;				
			}

			double ratio = betaS.sample()/ r;
			sum += ratio;
			sq += ratio*ratio;			
		}
		double mean = sum/numSamples;

		return new NormalDistribution(mean, Math.sqrt(sq/numSamples - mean*mean));
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {


	}	

	/**
	 * This class is copied from org.apache.commons.math3.distribution.SaddlePointExpansion
	 * so that some methods (such as logBinomialProbability) are accessible from here.
	 * <p>
	 * Utility class used by various distributions to accurately compute their
	 * respective probability mass functions. The implementation for this class is
	 * based on the Catherine Loader's <a target="_blank"
	 * href="http://www.herine.net/stat/software/dbinom.html">dbinom</a> routines.
	 * </p>
	 * <p>
	 * This class is not intended to be called directly.
	 * </p>
	 * <p>
	 * References:
	 * <ol>
	 * <li>Catherine Loader (2000). "Fast and Accurate Computation of Binomial
	 * Probabilities.". <a target="_blank"
	 * href="http://www.herine.net/stat/papers/dbinom.pdf">
	 * http://www.herine.net/stat/papers/dbinom.pdf</a></li>
	 * </ol>
	 * </p>
	 *
	 * @since 2.1
	 * @version $Id: SaddlePointExpansion.java 1244107 2012-02-14 16:17:55Z erans $
	 */
	static public final class SaddlePointExpansion {

		/** 1/2 * log(2 &#960;). */
		private static final double HALF_LOG_2_PI = 0.5 * FastMath.log(MathUtils.TWO_PI);

		/** exact Stirling expansion error for certain values. */
		private static final double[] EXACT_STIRLING_ERRORS = { 0.0, /* 0.0 */
			0.1534264097200273452913848, /* 0.5 */
			0.0810614667953272582196702, /* 1.0 */
			0.0548141210519176538961390, /* 1.5 */
			0.0413406959554092940938221, /* 2.0 */
			0.03316287351993628748511048, /* 2.5 */
			0.02767792568499833914878929, /* 3.0 */
			0.02374616365629749597132920, /* 3.5 */
			0.02079067210376509311152277, /* 4.0 */
			0.01848845053267318523077934, /* 4.5 */
			0.01664469118982119216319487, /* 5.0 */
			0.01513497322191737887351255, /* 5.5 */
			0.01387612882307074799874573, /* 6.0 */
			0.01281046524292022692424986, /* 6.5 */
			0.01189670994589177009505572, /* 7.0 */
			0.01110455975820691732662991, /* 7.5 */
			0.010411265261972096497478567, /* 8.0 */
			0.009799416126158803298389475, /* 8.5 */
			0.009255462182712732917728637, /* 9.0 */
			0.008768700134139385462952823, /* 9.5 */
			0.008330563433362871256469318, /* 10.0 */
			0.007934114564314020547248100, /* 10.5 */
			0.007573675487951840794972024, /* 11.0 */
			0.007244554301320383179543912, /* 11.5 */
			0.006942840107209529865664152, /* 12.0 */
			0.006665247032707682442354394, /* 12.5 */
			0.006408994188004207068439631, /* 13.0 */
			0.006171712263039457647532867, /* 13.5 */
			0.005951370112758847735624416, /* 14.0 */
			0.005746216513010115682023589, /* 14.5 */
			0.005554733551962801371038690 /* 15.0 */
		};

		/**
		 * Default constructor.
		 */
		private SaddlePointExpansion() {
			super();
		}

		/**
		 * Compute the error of Stirling's series at the given value.
		 * <p>
		 * References:
		 * <ol>
		 * <li>Eric W. Weisstein. "Stirling's Series." From MathWorld--A Wolfram Web
		 * Resource. <a target="_blank"
		 * href="http://mathworld.wolfram.com/StirlingsSeries.html">
		 * http://mathworld.wolfram.com/StirlingsSeries.html</a></li>
		 * </ol>
		 * </p>
		 *
		 * @param z the value.
		 * @return the Striling's series error.
		 */
		static double getStirlingError(double z) {
			double ret;
			if (z < 15.0) {
				double z2 = 2.0 * z;
				if (FastMath.floor(z2) == z2) {
					ret = EXACT_STIRLING_ERRORS[(int) z2];
				} else {
					ret = Gamma.logGamma(z + 1.0) - (z + 0.5) * FastMath.log(z) +
						z - HALF_LOG_2_PI;
				}
			} else {
				double z2 = z * z;
				ret = (0.083333333333333333333 -
					(0.00277777777777777777778 -
						(0.00079365079365079365079365 -
							(0.000595238095238095238095238 -
								0.0008417508417508417508417508 /
								z2) / z2) / z2) / z2) / z;
			}
			return ret;
		}

		/**
		 * A part of the deviance portion of the saddle point approximation.
		 * <p>
		 * References:
		 * <ol>
		 * <li>Catherine Loader (2000). "Fast and Accurate Computation of Binomial
		 * Probabilities.". <a target="_blank"
		 * href="http://www.herine.net/stat/papers/dbinom.pdf">
		 * http://www.herine.net/stat/papers/dbinom.pdf</a></li>
		 * </ol>
		 * </p>
		 *
		 * @param x the x value.
		 * @param mu the average.
		 * @return a part of the deviance.
		 */
		static double getDeviancePart(double x, double mu) {
			double ret;
			if (FastMath.abs(x - mu) < 0.1 * (x + mu)) {
				double d = x - mu;
				double v = d / (x + mu);
				double s1 = v * d;
				double s = Double.NaN;
				double ej = 2.0 * x * v;
				v = v * v;
				int j = 1;
				while (s1 != s) {
					s = s1;
					ej *= v;
					s1 = s + ej / ((j * 2) + 1);
					++j;
				}
				ret = s1;
			} else {
				ret = x * FastMath.log(x / mu) + mu - x;
			}
			return ret;
		}

		/**
		 * Compute the logarithm of the PMF for a binomial distribution
		 * using the saddle point expansion.
		 *
		 * @param x the value at which the probability is evaluated.
		 * @param n the number of trials.
		 * @param p the probability of success.
		 * @param q the probability of failure (1 - p).
		 * @return log(p(x)).
		 */
		static public  double logBinomialProbability(int x, int n, double p, double q) {
			double ret;
			if (x == 0) {
				if (p < 0.1) {
					ret = -getDeviancePart(n, n * q) - n * p;
				} else {
					ret = n * FastMath.log(q);
				}
			} else if (x == n) {
				if (q < 0.1) {
					ret = -getDeviancePart(n, n * p) - n * q;
				} else {
					ret = n * FastMath.log(p);
				}
			} else {
				ret = getStirlingError(n) - getStirlingError(x) -
					getStirlingError(n - x) - getDeviancePart(x, n * p) -
					getDeviancePart(n - x, n * q);
				double f = (MathUtils.TWO_PI * x * (n - x)) / n;
				ret = -0.5 * FastMath.log(f) + ret;
			}
			return ret;
		}
	}

}
