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

package japsadev.xm.genome;

public class AdaptiveMarkovExpertLong extends ExpertLong {
	private AdaptiveMarkovLong markov;

	public AdaptiveMarkovExpertLong(GenomeSequence seq, int order) {
		super(seq);
		markov = new AdaptiveMarkovLong(order);
		// posProb = 1.0 /3;
	}

	public void resurrect(GenomeSequence work, long i, int past) {
	}

	public void resign() {
	}

	public double probability(int character) {
		return markov.probability(character);
	}

	public int copyFrom(int i) {
		return 0;
	}

	public double update(int actual) {
		double cost = (markov.probability(actual));
		updateCost(cost);
		markov.update(actual);
		return cost;
	}

	public String toString() {
		return "ME";
	}

	public static void main(String[] args) {

	}

	public int copyFrom() {
		return 0;
	}

	public void learn() {
	}
}

class AdaptiveMarkovLong {
	int[] charCounts;
	int[] countTotal;
	int order;
	int currentInd = 0;// index of current context
	int MASK;// matrix size

	int[] history;
	int HIS_SIZE = 256;
	int backInd = 0;
	int ind = 0;

	public AdaptiveMarkovLong(int order) {
		this.order = order;
		MASK = (int) Math.pow(ExpertLong.ALPHABET_SIZE, order);

		history = new int[HIS_SIZE];
		if (order >= 0) {
			charCounts = new int[MASK * ExpertLong.ALPHABET_SIZE];
			countTotal = new int[MASK];
			for (int i = 0; i < charCounts.length; i++)
				charCounts[i] = 1;
			for (int i = 0; i < countTotal.length; i++)
				countTotal[i] = ExpertLong.ALPHABET_SIZE;
		}
	}

	boolean yes = false;

	public void update(int a) {
		// double res = encodeLen(a);
		countTotal[currentInd]++;// = DistributionExpert.ADD;
		currentInd = currentInd * ExpertLong.ALPHABET_SIZE + a;
		charCounts[currentInd]++;
		currentInd = currentInd % MASK;

		if (yes) {
			countTotal[backInd]--;
			backInd = backInd * ExpertLong.ALPHABET_SIZE + history[ind];
			charCounts[backInd]--;
			backInd = backInd % MASK;
			history[ind] = a;
			ind = (ind + 1) % HIS_SIZE;
		} else {
			history[ind] = a;
			ind = (ind + 1) % HIS_SIZE;
			if (ind == 0) {
				yes = true;
			}
		}
	}

	public double probability(int a) {
		return ((double) charCounts[currentInd * ExpertLong.ALPHABET_SIZE + a])
				/ countTotal[currentInd];
	}
}
