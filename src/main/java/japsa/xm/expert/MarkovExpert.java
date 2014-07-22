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

public class MarkovExpert extends Expert {
	private Markov markov;

	public MarkovExpert(int order) {
		this(order, 1);
	}
	
	public MarkovExpert(int order,double pseudocount) {
		super();
		markov = new Markov(order, pseudocount);
	}	

	public void learn(byte[] aSeq) {
		for (int i = 0; i < aSeq.length; i++) {
			markov.update(aSeq[i]);
		}

	}
	
	public double probability(int character) {
		return markov.probability(character);
	}


	public double update(int actual){
		double costActual = markov.probability(actual);
		updateCost(costActual);
		markov.update(actual);

		return costActual;

	}
	


	public String toString() {
		return "ME";
	}

	public void learn() {
	}

}

class Markov {
	private double[] charCounts;
	private double[] countTotal;	
	
	private int currentInd = 0;// index of current context	

	public Markov(int order) {
		this(order,1);
	}
	/**
	 * Create a Markov model with a order and a psuedocount
	 * @param order
	 * @param psudecount
	 */
	public Markov(int order, double pseudocount) {
		//this.order = order;
		int MASK = (int) Math.pow(Expert.alphabet().size(), order);
		if (order >= 0) {
			charCounts = new double[MASK * Expert.alphabet().size()];
			countTotal = new double[MASK];
			
			//Initlinise
			for (int i = 0; i < charCounts.length; i++)
				charCounts[i] = pseudocount;
			
			for (int i = 0; i < countTotal.length; i++)
				countTotal[i] = pseudocount * Expert.alphabet().size();
		}
	}

	public void update(int a) {
		countTotal[currentInd]++;
		currentInd = currentInd * Expert.alphabet().size() + a;
		charCounts[currentInd]++;
		currentInd = currentInd % countTotal.length;
	}

	public double probability(int a) {
		return  charCounts[currentInd * Expert.alphabet().size() + a]
				/ countTotal[currentInd];
	}
}


