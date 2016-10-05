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
 * 26/08/2016 - Minh Duc Cao: Start                                        
 *  
 ****************************************************************************/


package japsa.util;

import java.util.Random;

/**
 * A library for simulation
 * @author minhduc
 *
 */
public class Simulation {
	
	/**
	 * Generate a random seed if the input <=0
	 * @param seed
	 * @return
	 */
	public static int seed(int seed){
		if (seed <= 0)
			seed = new Random().nextInt();
		
		//make sure seed is not negative
		if (seed <0 )
			seed = - seed;
		
		return seed;
	}
	
	/**
	 * This function return a sample from log-logistic distribution (aka Fisk 
	 * distribution) with a scale parameter alpha and shape parameter beta.
	 * 
	 * See wikipedia on Fisk Distribution
	 * 
	 * @param alpha: scale parameter
	 * @param beta: shape parameter
	 * @param rnd: random generator
	 * @return
	 */
	public static double logLogisticSample(double alpha, double beta, Random rnd){
		double u = rnd.nextDouble();

		if (u <= 0.5)
			return alpha * Math.pow (u / (1.0 - u), 1.0 / beta);
		else
			return alpha / Math.pow ((1.0 - u)/ u, 1.0 / beta);	
	}
	
	public static double logLogisticPDF(double x, double alpha, double beta){
		x = x / alpha;
		return Math.pow(beta * x, -beta -1) *  Math.pow(( 1 + Math.pow(x, -beta)), -2);
	}	
}
