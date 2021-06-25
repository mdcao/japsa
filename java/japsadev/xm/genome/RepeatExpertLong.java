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

/**
 * @author Minh Duc Cao
 * 
 */
public abstract class RepeatExpertLong extends ExpertLong {

	public static int COPY_TYPE = 1;
	public static int PALIN_TYPE = -1;

	protected long start;// Starting point of expert on the sequence
	/**
	 * The current pointer of expert, japsa.seq[currentPointer] is expert
	 * prediction mobile elements
	 */
	protected long currentPointer;// the pointer to current position of expert
	protected int currentBase;// the base that the pointer points to

	/**
	 * The length of the expert, ideally this is the lenght of the mobile
	 * elements
	 */
	protected long length;// maximum possible length

	MyBitSetLong bitSet;
	protected int expertType = 1;// Copy or palindrome

	RepeatExpertLong(GenomeSequence seq, long start, MyBitSetLong b, int type) {
		super(seq);
		bitSet = b;
		this.start = start;
		currentPointer = start;// Currently point to the very first char
		currentBase = this.genSeq.getBase(currentPointer);

		expertType = type;

		// Get the length of the expert
		if (type == PALIN_TYPE)
			length = start - 1;// at most as long as this position
		else
			length = seq.getLength() - start - 1;
	}

	void reset(GenomeSequence seq, long startPos, MyBitSetLong b, int type) {
		reset(seq);
		bitSet = b;
		this.start = startPos;
		currentPointer = start;// Currently point to the very first char

		// currentBase = this.genSeq.getBase(currentPointer);

		expertType = type;

		// Get the length of the expert
		if (type == PALIN_TYPE)
			length = start - 1;// at most as long as this position
		else
			length = seq.getLength() - start - 1;
	}

	public void resign() {
		bitSet.clear(id);
	}

	public void setID(long lid) {
		id = lid;
	}

	public long getID() {
		return id;
	}

	public abstract RepeatExpertLong duplicate(GenomeSequence seq1,
			long start1, MyBitSetLong b);

	public abstract void reuseExpert(GenomeSequence seq1, long start1,
			MyBitSetLong b, int eType);

	protected abstract void computeProbs();

	public void resurrect(GenomeSequence workSeq, long pos, int past) {
		// resurrect(byte[] workSeq, int posSrc, int past) {
		bitSet.set(id);
		// make sure posSrc >=context_length
		currentPointer = start - past * expertType;
		currentBase = this.genSeq.getBase(currentPointer);
		computeProbs();

		for (long i = pos - past; i <= pos; i++) {
			update(workSeq.getBase(i));
		}
	}

	public long getStart() {
		return start;
	}

	public long getCurrentPointer() {
		return currentPointer;
	}

	public long getLength() {
		return length;
	}

	public int getExpertType() {
		return expertType;
	}
}
