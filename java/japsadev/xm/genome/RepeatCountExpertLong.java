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

/**
 * Written by Chris Mears, modified and maintained by Minh Duc Cao
 */
package japsadev.xm.genome;

public class RepeatCountExpertLong extends RepeatExpertLong {

	int countRight;
	int count;

	double[] probs = new double[ExpertLong.ALPHABET_SIZE];

	// int currentBase;

	public RepeatCountExpertLong(GenomeSequence seq, long start,
			MyBitSetLong b, int type) {
		super(seq, start, b, type);

		expertType = type;
		count = 1;
		countRight = 0;

	}

	/**
	 * A method resemble constructor, called to reuse experts
	 */
	public void reuseExpert(GenomeSequence seq, long startPos, MyBitSetLong b,
			int type) {
		reset(seq, startPos, b, type);//

		count = 1;
		countRight = 0;
	}

	protected void computeProbs() {
		double prob = (countRight + .01) / (count);
		if (prob <= 0 || prob >= 1) {
			System.err.println("Error " + prob + "  " + count + "   "
					+ countRight);
			(new Exception()).printStackTrace();
			System.exit(1);
		}

		double other = (1.0 - prob) / (ExpertLong.ALPHABET_SIZE - 1);

		for (int i = 0; i < ExpertLong.ALPHABET_SIZE; i++) {
			if (expertType == COPY_TYPE) {
				probs[i] = (i == this.currentBase) ? prob : other;
			} else {
				probs[i] = (i == (3 - this.currentBase)) ? prob : other;
			}
		}
	}

	public double probability(int character) {
		return probs[character];
		/**************************************************************
		 * double prob = (countRight + .01 ) / (count); if (prob <=0 || prob >=
		 * 1){ System.err.println("Error " + prob + "  " + count + "   " +
		 * countRight); (new Exception()).printStackTrace(); System.exit(1); }
		 * 
		 * int match = character;
		 * 
		 * if (expertType == PALIN_TYPE) match = 3 - character;
		 * 
		 * 
		 * if (currentBase == match) return prob; else return (1 - prob) /
		 * (ExpertLong.ALPHABET_SIZE - 1); /
		 **************************************************************/
	}

	public double update(int actual) {

		if ((currentPointer - start) * expertType >= length) {
			// System.out.println("hmmmm");
			return -1;
		}

		int match = actual;

		if (expertType == PALIN_TYPE)
			match = 3 - actual;

		// Recent prediction
		double prob = probability(actual);

		// Remove previous history element
		// countRight += ( genSeq.getBase(currentPointer) == match)? 1:0;
		countRight += (currentBase == match) ? 1 : 0;
		count++;

		currentPointer += expertType;
		currentBase = genSeq.getBase(currentPointer);
		computeProbs();

		updateCost(prob);

		return prob;
	}

	public RepeatExpertLong duplicate(GenomeSequence seq, long startPos,
			MyBitSetLong b) {
		return new RepeatCountExpertLong(seq, startPos, b, expertType);
	}

}
