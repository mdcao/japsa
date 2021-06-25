package japsadev.tools;

import org.apache.commons.math3.special.Gamma;

public class BetaBinomial {
	 double alpha, beta, trials;
	 public void set(double alpha, double beta, double trials){
		 this.alpha = alpha;
		 this.beta  = beta;
		 this.trials = trials;
	 }
	 public double logdensity(double k) {
         //int k = (int) Math.rint(x);
         
//    if (k < 0 | k > trials) return 0;
		 if( alpha < 1e-10){
	    	 return Double.NEGATIVE_INFINITY;
	    	/*  double res1 =Gamma.logGamma(k+alpha);
	    	  double res2 = Gamma.logGamma(trials-k+beta);
	    	  double res3 = Gamma.logGamma(alpha+beta);
	    	  double res4 = Gamma.logGamma(trials+2);
	    	  double res5 = Math.log(trials+1);
	    	  double res6 = Gamma.logGamma(alpha+beta+trials);
	    	  double res7 = Gamma.logGamma(alpha);
	    	  double res8 = Gamma.logGamma(beta);
	    	  double res9 = Gamma.logGamma(k+1);
	    	  double res10 = Gamma.logGamma(trials-k+1);
	    	 throw new RuntimeException("os na");*/
	     }
 double res = (Gamma.logGamma(k+alpha)+Gamma.logGamma(trials-k+beta)+Gamma.logGamma(alpha+beta)+Gamma.logGamma(trials+1)) - 
         (
        		 //	Math.log(trials+1)+
        		 		Gamma.logGamma(alpha+beta+trials)+Gamma.logGamma(alpha)+Gamma.logGamma(beta)+Gamma.logGamma(k+1)+Gamma.logGamma(trials-k+1));
 
 //double res = (Gamma.logGamma(k+alpha)+Gamma.logGamma(trials-k+beta)+Gamma.logGamma(alpha+beta)+Gamma.logGamma(trials+2)) 
/* if(Double.isNaN(res)){
   throw new RuntimeException("!!");
}*/
 	return res ;//+ norm;
 
}
}
