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

//import java.util.Random;

/**
 * @author Minh Duc Cao
 * 
 */
public class PrefixArray extends PrefixArrayAbstract {

	public PrefixArray(byte[] aSeq) {
		this(aSeq, false);

	}

	public PrefixArray(byte[] aSeq, boolean isReverse) {
		mySeq = aSeq;
		initilise(isReverse);
	}

	public void initilise(boolean isReverse) {
		int lastIndex = mySeq.length - 1;
		if (!isReverse) {
			// Reverse the sequence
			for (int i = 0; i < mySeq.length / 2; i++) {
				byte tmpB = mySeq[i];
				mySeq[i] = mySeq[lastIndex - i];
				mySeq[lastIndex - i] = tmpB;

			}
		}
		pa = SuffixArrayConstruction.construct(mySeq);// Suffix Array

		// Convert the suffix array to inverse prefix array
		// pa = new int[mySeq.length];
		invPA = new int[mySeq.length];

		for (int i = 0; i < mySeq.length; i++) {
			// Reverse s
			pa[i] = lastIndex - pa[i + 1];
			invPA[pa[i]] = i;
		}

		// reverse mySequence back
		for (int i = 0; i < mySeq.length / 2; i++) {
			byte tmpB = mySeq[i];
			mySeq[i] = mySeq[lastIndex - i];
			mySeq[lastIndex - i] = tmpB;

		}
		/***/
		startRange = 0;
		endRange = mySeq.length;
		matchLength = 1;// No match at the moment
	}

	// Impplement parent abstract
	protected boolean isMatched(byte b1, byte b2) {
		// boolean x = (b1 - b2) % 2 == 0;
		return b1 == b2;
	}

	/**
	 * Construct the SA object: compute the pa, lcp etc
	 * 
	 * @param aSeq
	 */

	public void print() {
		System.out.printf("\n   i  lcp pa in  \n");

		for (int i = 0; i < mySeq.length; i++) {
			System.out.printf("%4d%4d%4d%4d ", i, lcp[i], pa[i], invPA[i]);
			for (int j = pa[i]; j >= 0; j--) {
				System.out.printf("%d", mySeq[j]);
			}
			System.out.println();
		}
	}

}
