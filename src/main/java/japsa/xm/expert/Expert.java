/*****************************************************************************
 * Copyright (c) 2010 Minh Duc Cao, Monash University.  All rights reserved. *
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
 * 3. Neither the name of Monash University nor the names of its contributors*
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

package japsa.xm.expert;

import japsa.seq.Alphabet;


public abstract class Expert {

	public static int CONTEXT_LENGTH = 20;// Just a defaut value
	public static int HASH_SIZE = 11;
	
	// Evaluation of experts
	protected double posteriorProb = 1;
	// History of prediction
	private double[] history;
	private int ind = 0;
	private static Alphabet alphabet = null;

	/**
	 * Set the dna only it had not been set
	 * 
	 * @param anAlphabet
	 */
	public static void setAlphabet(Alphabet anAlphabet) {
		if (alphabet != null && alphabet != anAlphabet) {
			throw new RuntimeException("Alphabet aready set to " + alphabet);
		}
		alphabet = anAlphabet;

	}

	public static Alphabet alphabet() {
		return alphabet;
	}


	public Expert() {
		this(1.0);
		
	}
			
	public Expert(double prior) {
		history = new double[CONTEXT_LENGTH];
		posteriorProb = prior;

		//pre-compute this so no need to recompute 1.0/seq.alphabet().size() again
		history[0] =  1.0/ alphabet.size();
		posteriorProb *= history[0]; 

		for (int i = 1; i < CONTEXT_LENGTH; i++) {
			history[i] = history[0]; 
			posteriorProb *= history[i];
		}
	}	
	
	/**
	 * 
	 * @param prob
	 */
	protected void updateCost(double prob) {
		posteriorProb = posteriorProb * prob / history[ind];
		history[ind] = prob;

		ind = (ind + 1) % history.length;
	}

	
	// More is actually more
	public double posteriorProb() {
		return posteriorProb;
	}

	
	/**
	 * The probability of current charactor
	 * 
	 * @param character
	 * @return
	 */
	public abstract double probability(int character);

	/**
	 * When update a position, the expert need to get the cost of current
	 * prediction, this cost needed for rating expert return neg if out of bound
	 */

	public abstract double update(int actual);	

	
	
	public static void setContext(int context) {
		CONTEXT_LENGTH = context;
	}

	public static void setParams(int alSize, int context) {		
		setContext(context);
	}
}
