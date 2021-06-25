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

/*                           Revision History                                
 * 28/05/2014 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsadev.tools.misc;

import java.io.BufferedReader;

import japsa.seq.SequenceReader;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathUtils;

/**
 * @author minhduc
 *
 */
public class VNTRReadDepth {	
	public static double logBin(int k, int n, double p){
		return 
		logBinomialProbability(k,n, p, 1.0 - p);		 
	}
	
	public static void main(String[] args) throws Exception{		
		String totFile = "/home/minhduc/Projects/data/Genomes/TB/collections/TB_samples.fdepth";
		String depthFile = "/home/minhduc/Projects/data/Genomes/TB/collections/TB_samples.rdepth";
		String lengthFile = "/home/minhduc/Projects/data/Genomes/TB/collections/TB_samples.length";		
		
		BufferedReader inAns = SequenceReader.openFile(totFile);
		int [][] tot    = new int[24][32];
		int [][] depth  = new int[24][32];
		int [][] length = new int[24][32];
		
		
		BufferedReader in = SequenceReader.openFile(totFile);
		in.readLine();//header
		for (int r = 0; r < 24; r++){
			String [] toks = in.readLine().trim().split("\t");
			for (int l = 0; l < 32;l++)
				tot[r][l] = Integer.parseInt(toks[l]);
		}
		in.close();
		
		in = SequenceReader.openFile(depthFile);
		in.readLine();//header
		for (int r = 0; r < 24; r++){
			String [] toks = in.readLine().trim().split("\t");
			for (int l = 0; l < 32;l++)
				depth[r][l] = Integer.parseInt(toks[l]);
		}
		in.close();
		
		in = SequenceReader.openFile(lengthFile);
		in.readLine();//header
		for (int r = 0; r < 24; r++){
			String [] toks = in.readLine().trim().split("\t");
			for (int l = 0; l < 32;l++)
				length[r][l] = Integer.parseInt(toks[l]);
		}
		in.close();
		
		int index = 5;
		
		int wrong = 0, right = 0;
		for (int r = 0; r < 24; r++){
			for (int c = 6; c < 32; c++){
				double p =  (depth[r][c] * 1.0 / tot[r][c]) / (depth[r][index] * 1.0 / tot[r][index]);
				double t = length[r][c] * 1.0 / length[r][index];
				System.out.printf("%4.2f vs %4.2f\t",p,t);
				if (0.9 < t && t < 1.1){//Nova
					if (0.8 < p && p < 1.2) right ++;
					else wrong ++;
				}
				if (t>=1.1){
					if (p >= 1.1) right ++;
					else wrong ++;
				}
				if (t<=0.9){
					if (p <0.9) right ++;
					else wrong ++;
				}
				
			}//for c
			System.out.println();
		}//for r		
		System.out.println(right + "  " + wrong);
		inAns.close();
	}
	
	
	
	public static void main1(String[] args) throws Exception{
		int readLength = 100;
		
		String ansFile = "/home/minhduc/Projects/data/Genomes/TB/collections/Miru24.dat";
		String datFile = "/home/minhduc/Projects/data/Genomes/TB/collections/TB_samples.depths";
				
		BufferedReader inAns = SequenceReader.openFile(ansFile);
		
		String aLine = inAns.readLine().trim();
		
		String [] strains = {"H37Rv", "Erdman", "Haarlem", "KZN_1435", "W_148"};
		double [] strainType = new double[strains.length];
		inAns.readLine();		
		
		BufferedReader inDat = SequenceReader.openFile(datFile);
		String line = inDat.readLine().trim().substring(3);
		
		String [] samples = line.split("\t");
		
		double [] tot = new double[samples.length];
		double [] depth = new double[samples.length];		
		
		int numRep = 24;
		
		//Read in total depth
		line = inDat.readLine().trim();
		String[] toks = line.split("\t");
		for (int i = 0; i < tot.length; i++)
			tot[i] = Double.parseDouble(toks[i]) / readLength;
				
		int rightP = 0, wrongP = 0; 
		int noChange = 0;
		for (int r = 0 ; r < numRep; r++){
			//get rep information
			aLine = inAns.readLine().trim();
			toks = aLine.split("\t");
			
			double repeatRef  = Double.parseDouble(toks[2]) - Double.parseDouble(toks[1]);
			double repeatUnit = Double.parseDouble(toks[3]);
			
			for (int i = 0; i< strains.length; i++)
				strainType[i] = Double.parseDouble(toks[6+i]);		
			
			//Read read depth			
			line = inDat.readLine().trim();
			toks = line.split("\t");
			for (int i = 0; i< depth.length; i++){
				depth[i] = Double.parseDouble(toks[i]) / readLength;				
			}
			
			//Analysis
			int ref = 4;
			for (int i = 5; i< depth.length; i++){
				//System.out.print(depthBin2(depth[i], tot[i], depth[ref], tot[ref],1)+"," + depthBin2(depth[i], tot[i], depth[ref], tot[ref],100)+"\t");
				
				int index = 0;
				if (samples[i].startsWith("E"))
					index = 1;
				else if (samples[i].startsWith("H"))
					index = 2;
				else if (samples[i].startsWith("K"))
					index = 3;
				else if (samples[i].startsWith("W"))
					index = 4;
				
				double actualLength = repeatRef + repeatUnit * strainType[index];
				//expct actualLength / repeatRef ~ (Ds/Ts)/(Dr/Tr)
				//System.out.print(actualLength/repeatRef + "," + actualLength + "," +repeatRef + "," + ((depth[i]/tot[i])/(depth[ref]/tot[ref]))+"," + strainType[index] +  
				//		"#" + (depth[i]+1) + "," + (tot[i]+1) + "," + (depth[ref] + 1) + "," + (tot[ref] + 1)+"\t");
					
				double ratio = depthBeta(depth[i], tot[i], depth[ref], tot[ref]);
				if (actualLength / repeatRef > 1.1) {
					if (ratio > 1){
						rightP ++;
					}else
						wrongP ++;					
				}else if (repeatRef / actualLength   > 1.1) {
					if (ratio > 1){
						wrongP ++;
					}else
						rightP ++;					
				} else{
					//if (ratio > 1.1 || ratio < 0.9)
					//	wrongP ++;
					//else
					//	rightP++;
						
				}
				
				
				System.out.print(actualLength/repeatRef + " vs " + depthBeta(depth[i], tot[i], depth[ref], tot[ref]) + "\t");
							
				
				//double odd = - depthBin(depth[i], tot[i], depth[ref], tot[ref],1) + depthBin3(depth[i], tot[i], depth[ref], tot[ref],actualLength, repeatRef,1);
				//if (odd > 0.001) {
				//	rightP ++;
				//}
				//if (odd < -0.001) {
				//	wrongP ++;
				//}				 
			}			
			System.out.println();			
		}					
		System.out.println(rightP + " vs " + wrongP + " vs " + noChange);
	}
	
	public static double depthPoisson(int depthS, int totDepthS, double depthR, double totDepthR, int numIt){
		double lamda = totDepthS * depthR / totDepthR;		
		//if (numIt <=1){
			PoissonDistribution pd = new PoissonDistribution(lamda);
			return pd.probability(depthS);			
		//}		
	}
	
	public static double depthBeta(double depthS, double totDepthS, double depthR, double totDepthR){
		double ps = (depthS + 1) / (totDepthS + 1);
		double pr = (depthR + 1) / (totDepthR + 1);
		
		return  ps/ pr;		
	}
	
	//public static double depthBin4(int depthS, double totDepthS, double depthR, double totDepthR, double s, double r, int numIt){
	//	
	//	
	//}
	
	/**
	 * log likelihood of P(dS | totS, dR, totR) ~ Bin(dS|totS, p) where p ~ Beta(dR + 1, totR - dR + 1)
	 * @param depthS
	 * @param totDepthS
	 * @param depthR
	 * @param totDepthR
	 * @return
	 */
	public static double depthBin(double depthS, double totDepthS, double depthR, double totDepthR, int numIt){
		double p = depthR / totDepthR;
		if (numIt <=1)
			return - logBin((int)depthS, (int)totDepthS, p);
		
		BetaDistribution beta = new BetaDistribution(depthR + 1, totDepthR - depthR + 1);
		int sum = 0;
		int ids = (int) depthS;
		int its = (int) depthS;
		for (int i = 0; i < numIt; i++){
			p = beta.sample();
			sum -=  logBin(ids, its, p);
		}
		//Note: this is the average of logProb, LC's version is average of prob
		return sum / numIt;		
	}
	

	/**
	 * P(dS | totS, dR, totR) ~ Bin (dS|dS + dR, p) where p ~ Beta(totS+1, totR+1)
	 * @param depthS
	 * @param totDepthS
	 * @param depthR
	 * @param totDepthR
	 * @param numIt
	 * @return
	 */
	public static double depthBin2(double depthS, double totDepthS, double depthR, double totDepthR, int numIt){		
		int ids = (int) depthS;		
		int idr = (int) depthR;		
		
		double p = totDepthS / (totDepthS + totDepthR);
		
		if (numIt <=1)
			return -logBin(ids, (ids + idr), p);
		
		BetaDistribution beta = new BetaDistribution(totDepthS + 1, totDepthR + 1);
		int sum = 0;
		for (int i = 0; i < numIt; i++){
			p = beta.sample();
			sum -=  logBin(ids, (ids + idr), p);
		}
		//Note: this is the average of logProb, LC's version is average of prob
		return sum / numIt;		
	}
	

	public static double depthBin3(double depthS, double totDepthS, double depthR, double totDepthR, double s, double r, int numIt){
		double p = (depthR / totDepthR) * (s/r);
		
		if (numIt <=1)
			return - logBin((int)depthS, (int) totDepthS, p);
		
		BetaDistribution beta = new BetaDistribution(depthR + 1, totDepthR - depthR + 1);
		int sum = 0;
		int ids = (int) depthS;
		int its = (int) totDepthS;
		for (int i = 0; i < numIt; i++){
			p = beta.sample() * s/r;
			sum -=  logBin(ids, its, p);
		}
		//Note: this is the average of logProb, LC's version is average of prob
		return sum / numIt;		
	}
	
	/************************************************************************
	double countnormal = 1, totnormal = 2;

	public double probability(double lrr, double baf) {
		double cellularity = vals[0]; 
		double ratio = vals[1];
		if(ratioAsLevels) ratio = ratios[backCN]*ratio;
		else if(cellAsLevels){
				cellularity = Math.min(1.0, ratios[backCN]*cellularity);
				//if(cellularity>1) cellularity = 1.0/cellularity;
		}
	    double 	mult = ((rcn/ratio)*cellularity) + (1-cellularity);
		
		double k = baf;//counttumour;
		double n = lrr;// tottumour;
		
	    if(numIt==1){
	    	  double p =  (countnormal/totnormal)* mult;
	    	  double v = k *Math.log(p)+(n-k)*Math.log(1-p);
	    	  return  v; 
	    	
	    }else{
	    	 double maxv = Double.NEGATIVE_INFINITY;
	    	 b.setState((countnormal+1)/betaDownWeight, ((totnormal-countnormal)+1)/betaDownWeight);
		    for(int j = 0; j<numIt; j++){
				  double p =  b.nextDouble()* mult;
				  double v = k *Math.log(p)+(n-k)*Math.log(1-p);
				
				  if(v>maxv) maxv = v;
				  logprobs[j] = v;
		    }
			double sum=0;
			for(int j =0; j<numIt; j++){
				sum+=Math.exp(logprobs[j]-maxv);
			}
			return (Math.log(sum/numIt)+maxv);
	    }
	    }
	/************************************************************************/	    


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
    static double logBinomialProbability(int x, int n, double p, double q) {
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
