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

import japsa.seq.AbstractSequence;
import japsa.util.MyBitSet;

/**
 * @author Minh Duc Cao
 * 
 */
public abstract class RepeatExpert extends Expert {

	public static int COPY_TYPE = 1;
	public static int PALIN_TYPE = -1;

	
	protected int start;// Starting point of expert on the sequence

	/**
	 * The current pointer of expert, japsa.seq[currentPointer] is expert prediction
	 * mobile elements
	 */
	protected int currentPointer;// the pointer to current position of expert

	/**
	 * The length of the expert, ideally this is the lenght of the mobile
	 * elements
	 */
	protected int length;// maximum possible length
	
	MyBitSet bitSet;
	protected int expertType = 1;// Copy or palindrome
	
	protected AbstractSequence seq;
	
	protected int id;// ID to identify it self
	public int getID() {
		return id;
	}

	public void setID(int id) {
		this.id = id;
	}
	private int counter = 0;
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


	protected RepeatExpert(AbstractSequence seq, int start, MyBitSet b, int type) {
		super();
		this.seq = seq;
		bitSet = b;
		this.start = start;
		currentPointer = start;// Currently point to the very first char

		expertType = type;

		// Get the length of the expert
		if (type == PALIN_TYPE)
			length = start - 1;// at most as long as this position
		else
			length = seq.length() - start - 1;
	}

	public void resign() {
		bitSet.clear(id);
	}

	public abstract RepeatExpert duplicate(AbstractSequence seq1, int start1, MyBitSet b);

/***	
	public void resurrect(Sequence workSeq, int pos, int past, double[] markovCost) {

		bitSet.set(id);
		// make sure pos >=context_length
		currentPointer = start - past * expertType;

		for (int i = pos - past; i <= pos; i++) {
			double prob = update(workSeq.symbolAt(i));
			this.infoGain += markovCost[i] + JapsaMath.log2(prob);
		}

	}
	****************/
	

	public void resurrect(AbstractSequence workSeq, int pos, int past) {
		bitSet.set(id);
		// make sure pos >=context_length
		currentPointer = start - past * expertType;

		for (int i = pos - past; i <= pos; i++) {
			update(workSeq.symbolAt(i));
		}
	}

	public int getStart() {
		return start;
	}

	public int getCurrentPointer() {
		return currentPointer;
	}

	public int getLength() {
		return length;
	}

	public int getExpertType() {
		return expertType;
	}
	
	public AbstractSequence getSeq(){
		return seq;
	}
}
