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

package japsadev.xm.hash;

import japsa.util.IntIterator;
import japsa.xm.hash.PatternStore;

/**
 * @author Minh Duc Cao
 * 
 */
public abstract class PrefixArrayAbstract implements PatternStore {
	private int hashSize = 10;
	// Suffix Array
	int[] pa;
	// and its invert
	int[] invPA;// = mySeq[sa[i] - 1];The charracter at the position previous of
				// suffix i (i)
	// A pointer to the sequence being indexed
	byte[] mySeq;
	// LCP
	int[] lcp;

	int startRange, endRange;
	int nextStartRange, nextEndRange;
	int matchLength;

	int[] posMatches = new int[1024];// Store all matches for the next proposing
	int posIndex = 0;

	// The sequence to match
	byte seqToMatch[];
	int matchIndex;// The index of the current sequence

	public PrefixArrayAbstract() {
		// Nothing
	}

	public void setHashSize(int hs) {
		this.hashSize = hs;
	}

	public void setSeqToMatch(byte[] anSeq) {
		seqToMatch = anSeq;
		matchIndex = -1;
	}

	void computeLcp() {
		lcp = new int[mySeq.length];
		int i, h = 0;
		h = 0; /* visit in string order */
		for (i = mySeq.length - 1; i >= 0; i--) { /* omit last, least suff */
			int x = invPA[i];
			if (x > 0) {
				int j = pa[x - 1];

				int p1 = i - h;
				int p0 = j - h;

				while (p1 >= 0 && p0 >= 0
						&& isMatched(mySeq[p1--], mySeq[p0--]))
					h++;

				// if (lcp[x] != h)
				// System.out.println("NO at i = " + i + " x = " + x + " h = " +
				// h + " lcp = " + lcp[x]);

				lcp[x] = h;
				if (h > 0)
					h--;
			}
		}
		lcp[0] = 0; /* least suffix has no predecessor */
	}

	public void printSummary() {
	};

	abstract protected boolean isMatched(byte b1, byte b2);

	// {
	// return b1 == b2;
	// }

	protected void shrink(int nextByte) {
		if (nextStartRange > nextEndRange) {
			this.nextStartRange = 0;
			this.nextEndRange = mySeq.length - 1;

			while (!isMatched(mySeq[pa[this.nextStartRange]],
					seqToMatch[matchIndex]))
				this.nextStartRange++;

			while (!isMatched(mySeq[pa[this.nextEndRange]],
					seqToMatch[matchIndex]))
				this.nextEndRange--;

			matchLength = 0;
		}

		matchLength++;// Add another char to the match
		startRange = nextStartRange;// startNext[nextByte];
		endRange = nextEndRange;// endNext[nextByte];

		nextStartRange = mySeq.length;
		nextEndRange = 0;
	}

	byte prevSymbol(int aPos) {
		if (pa[aPos] == mySeq.length - 1)
			return -1;
		else
			return mySeq[pa[aPos] + 1];
	}

	int prePosition(int aPos) {
		if (pa[aPos] == mySeq.length - 1)
			return 0;
		else
			return invPA[pa[aPos] + 1];
	}

	/**
	 * The next byte from the current encoding symbol
	 * 
	 * @param nextByte
	 */
	public void nextKey(int nextByte) {
		// Try to make the next shrink valid. you cant shrink while the next
		// range is not valid
		while (nextEndRange < nextStartRange) {
			if (!expand())
				break;
		}

		matchIndex++;// Another one
		// shrinking
		shrink(nextByte);

		// expand expand
		browse();
	}

	// Shrink the match length to get more candidates
	private boolean expand() {
		matchLength--;
		if (matchLength == 0)
			return false;// The character has not seen before

		int tmp;
		// expand upward
		while (lcp[startRange] >= matchLength) {
			startRange--;

			// Put the new value into the nextRange
			if (posIndex < posMatches.length)
				posMatches[posIndex++] = pa[startRange];

			// To see if the nextRange valide
			if (isMatched(seqToMatch[matchIndex + 1], prevSymbol(startRange))
					&& (tmp = prePosition(startRange)) < nextStartRange)
				nextStartRange = tmp;
		}
		// expand download

		while (endRange + 1 < mySeq.length && lcp[endRange + 1] >= matchLength) {
			endRange++;

			// Put the new value into the nextRange
			if (posIndex < posMatches.length)
				posMatches[posIndex++] = pa[endRange];
			// To see if the nextRange valide
			if (isMatched(seqToMatch[matchIndex + 1], prevSymbol(endRange))
					&& (tmp = prePosition(endRange)) > nextEndRange)
				nextEndRange = tmp;
		}

		return true;

	}

	// Browse throw the list of all matches

	public void browse() {
		posIndex = 0;
		int tmp;
		for (int i = startRange; i <= endRange; i++) {
			// Put into matches table
			if (posIndex < posMatches.length)
				posMatches[posIndex++] = pa[i];

			if (isMatched(seqToMatch[matchIndex + 1], prevSymbol(i))) {
				if ((tmp = prePosition(i)) > nextEndRange)
					nextEndRange = tmp;
				if ((tmp = prePosition(i)) < nextStartRange)
					nextStartRange = tmp;
			}
			// Update nextStart/nextEnd range
		}
	}

	/**
	 * Construct the SA object: compute the pa, lcp etc
	 * 
	 * @param aSeq
	 */

	abstract public void print();

	// Methods required by PatternStore interface
	public void clear() {
	};

	public void putCurrentValue(int val) {
	};

	/**
	 * Return the iterator of all previous candidate copy experts
	 * 
	 * @return
	 */
	public IntIterator copyIterator() {
		return null;
	};

	/**
	 * Return the iterator of all previous candidate palin experts
	 * 
	 * @return
	 */
	public IntIterator palinIterator() {
		return null;
	};

	/**
	 * Return the iterator of all previous candidate copy and palin experts They
	 * are in fact the combination of both copy and palin experts
	 * 
	 * @return
	 */
	public IntIterator iterator() {
		return new CopyIterator();
	}

	class CopyIterator implements IntIterator {
		int myIndex;

		public CopyIterator() {
			myIndex = 0;
		}

		public boolean hasNext() {

			if (hashSize > matchLength)
				return false;

			if (myIndex < posIndex)
				return true;

			/****************************************************************/
			while (hashSize < matchLength) {
				expand();
				if (myIndex < posIndex)
					return true;
			}
			/****************************************************************/
			return false;
		}

		public int sizeAvailable() {
			return posIndex;
		}

		public int next() {
			return posMatches[myIndex++];
		}
	}
}
