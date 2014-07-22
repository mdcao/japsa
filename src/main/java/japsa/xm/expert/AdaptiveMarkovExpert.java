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


public class AdaptiveMarkovExpert extends Expert {
	private AdaptiveMarkov markov;

	public AdaptiveMarkovExpert(int order, int historySize) {
		this(order, historySize, 1.0);
	}
	
	public AdaptiveMarkovExpert(int order, int historySize,  double pseudocount) {
		super();
		markov = new AdaptiveMarkov(order, historySize, pseudocount);
	}
	

	public double probability(int character) {
		return markov.probability(character);
	}

	public double update(int actual) {
		double cost = (markov.probability(actual));
		updateCost(cost);
		markov.update(actual);
		return cost;
	}
	
	
	public String toString() {
		return "AME";
	}
}



class AdaptiveMarkov {
	private double[] charCounts;
	private double[] countTotal;
	
	
	private int currentInd = 0;// index of current context	

	int[] history;
	int backInd = 0;
	int ind = 0;
	
	
	public AdaptiveMarkov(int order, int historySize) {
		this(order, historySize, 1.0);
	}
	/**
	 * Create a Markov model with a order and a psuedocount
	 * @param order
	 * @param psudecount
	 */
	public AdaptiveMarkov(int order, int historySize, double psuedecount) {
		//this.order = order;
		int MASK = (int) Math.pow(Expert.alphabet().size(), order);
		if (order >= 0) {
			charCounts = new double[MASK * Expert.alphabet().size()];
			countTotal = new double[MASK];
			
			history = new int[historySize];			
			
			//Initlinise
			for (int i = 0; i < charCounts.length; i++)
				charCounts[i] = psuedecount;
			
			for (int i = 0; i < countTotal.length; i++)
				countTotal[i] = psuedecount * Expert.alphabet().size();
		}
	}

	private boolean yes = false;
	public void update(int a) {
		countTotal[currentInd]++;
		currentInd = currentInd * Expert.alphabet().size() + a;
		charCounts[currentInd]++;
		currentInd = currentInd % countTotal.length;
		
		//Remove those outside the history window
		if (yes) {
			countTotal[backInd]--;
			backInd = backInd * Expert.alphabet().size() + history[ind];
			charCounts[backInd]--;
			backInd = backInd % countTotal.length;
			history[ind] = a;
			ind = (ind + 1) % history.length;
		} else {
			history[ind] = a;
			ind = (ind + 1) % history.length;
			if (ind == 0) {
				yes = true;
			}
		}
	}

	public double probability(int a) {
		return  charCounts[currentInd * Expert.alphabet().size() + a]
				/ countTotal[currentInd];
	}
}

