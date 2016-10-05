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

import japsa.util.Distribution;
import japsa.util.JapsaMath;

public abstract class ExpertLong {

	// global variables, params for the algorithm
	protected static int CONTEXT_LENGTH = 20;// Just a defaut value
	public static int ALPHABET_SIZE = 4;// should be 4 for dna
	public static double INV_ALPHABET_SIZE = .25;
	public static int HASH_SIZE = 11;
	private static double[] DEF_HISTORY;
	private static double DEF_MSGLEN;

	//
	ExpertLong next = null; // to implement a linked list of experts
	protected long id;// ID to identify it self
	private int counter = 0;

	GenomeSequence genSeq;
	// private double prior = 1.0;

	// Evaluation of experts
	protected double msgLenProb = 1;
	// private double prior = 1;

	// History of prediction
	private double[] history;
	private int ind = 0;

	/**
	 * For summary of experts
	 */
	// #DEBUG_BEGIN
	// public static double baseCost = 0;

	// debug
	// protected double encodeCost = 0;
	// protected int encodeCount = 0;
	// protected double encodeCostLuck = 0;
	// #DEBUG_END

	// public double infoGain = 0.0;
	// A dummy one
	public ExpertLong() {
	}

	// public double getPrior() {
	// return prior;
	// }

	// public void setPrior(double prior) {
	// this.prior = prior;
	// }

	public ExpertLong(GenomeSequence gen) {
		this.genSeq = gen;
		msgLenProb = DEF_MSGLEN;//

		history = new double[CONTEXT_LENGTH];

		for (int i = 0; i < CONTEXT_LENGTH; i++) {
			history[i] = INV_ALPHABET_SIZE;// DEFAULT_LEN;//JapsaMath.log2(ALPHABET_SIZE);
			// posteriorProb *= history[i];
		}
	}

	public void reset(GenomeSequence gen) {
		this.genSeq = gen;
		msgLenProb = DEF_MSGLEN;//

		for (int i = 0; i < CONTEXT_LENGTH; i++) {
			history[i] = INV_ALPHABET_SIZE;// DEFAULT_LEN;//JapsaMath.log2(ALPHABET_SIZE);
		}
		ind = counter = 0;
	}

	public void updateCost(double prob) {
		msgLenProb = msgLenProb * prob / history[ind];
		history[ind] = prob;

		ind = (ind + 1) % history.length;
	}

	// Rating based on msg Length
	// More is less, less is more
	public double rate() {
		// return msgLen;
		return -JapsaMath.log2(msgLenProb);

	}

	// More is actually more
	public double rateProb() {
		return msgLenProb;
	}

	public long getID() {
		return id;
	}

	public void setID(long id) {
		this.id = id;
	}

	// Operate on linked list
	public ExpertLong getNext() {
		return next;
	}

	public void setNext(ExpertLong next) {
		this.next = next;
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

	// Should return the current rate after resurrect
	// Resurrect when working on workSeq, at position posSrc
	public abstract void resurrect(GenomeSequence workSeq, long pos, int past);

	public abstract void resign();

	/**
	 * Counter
	 * 
	 * @return
	 */
	public int getCounter() {
		return counter;
	}

	public void setCounter(int counter) {
		this.counter = counter;
	}

	public void resetCounter() {
		counter = 0;
	}

	public void incrementCounter() {
		counter++;
	}

	// These methods are for setting global parameters
	public static void setAphabetSize(int size) {
		ALPHABET_SIZE = size;
		INV_ALPHABET_SIZE = 1.0 / ALPHABET_SIZE;
	}

	public static void setContext(int context) {
		CONTEXT_LENGTH = context;
	}

	public static void setParams(int alSize, int context) {
		setAphabetSize(alSize);
		setContext(context);

		DEF_HISTORY = new double[CONTEXT_LENGTH];
		DEF_MSGLEN = 1.0;

		for (int i = 0; i < CONTEXT_LENGTH; i++) {
			DEF_HISTORY[i] = INV_ALPHABET_SIZE;// DEFAULT_LEN;//JapsaMath.log2(ALPHABET_SIZE);
			DEF_MSGLEN *= DEF_HISTORY[i];
		}
	}

	public Distribution getDistribution() {
		Distribution dist = new Distribution(ExpertLong.ALPHABET_SIZE);
		for (byte i = 0; i < ExpertLong.ALPHABET_SIZE; i++)
			dist.setWeight(i, this.probability(i));

		return dist;
	}

}
