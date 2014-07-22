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

package japsa.xm.hash;

import java.util.StringTokenizer;

public class GappedHashtable extends MyHashtable {

	// /////////////////////////////////////////////////////////////
	private long currentLongKey = 0;
	private long currentPalinLongKey = 0;

	long[] islands;
	int[] islandSize;

	// /////////////////////////////////////////////////////////////

	String mask;

	public GappedHashtable(String aMask) {
		this(aMask, 2);
	}

	public GappedHashtable(String aMask, int bpp) {
		super();
		this.initilise(aMask, bpp);

	}

	public void initilise() {
		super.initilise();
	}

	public void initilise(String aMask, int bpp) {
		super.initilise(hashLength(aMask), bpp);
		mask = aMask;

		preCompute = (_bitPerSymbol * (mask.length() - 1));

		mask = aMask;

		// mask should be in form [01]+
		StringTokenizer st = new StringTokenizer(mask, "0");
		islands = new long[st.countTokens()];
		islandSize = new int[st.countTokens()];

		int j = -1;
		int last = 0;
		int count = 0;// Count of Zeros

		// Break the mask into islands
		for (int i = 0; i < mask.length(); i++) {
			if (mask.charAt(mask.length() - i - 1) == '1') {
				if (last == 0) {
					j++;
					islandSize[j] = count * 2;
					islands[j] = (1l << (i * 2)) | (1l << (i * 2 + 1));
					last = 1;
				} else
					islands[j] = islands[j] | (1l << (2 * i))
							| (1l << (i * 2 + 1));

			} else {
				last = 0;
				count++;
			}
		}
		// Construct maskP
		StringBuffer sb = new StringBuffer();
		for (int i = mask.length() - 1; i >= 0; i--) {
			sb.append(mask.charAt(i));
		}

		System.out.println(mask);
		for (j = 0; j < islands.length; j++) {
			System.out.printf("%s %d \n", Long.toBinaryString(islands[j]),
					islandSize[j]);
		}

	}

	public static int hashLength(String aMask) {
		int _count = 0;
		for (int i = 0; i < aMask.length(); i++) {
			if (aMask.charAt(i) == '1')
				_count++;
		}

		return _count;
	}

	/**
	 * Compute the new key as a new base is read in
	 * 
	 * @param baseInd
	 */

	// TODO
	public void nextKey(byte baseInd) {
		// System.out.printf("%2d ",baseInd);
		currentLongKey <<= _bitPerSymbol;
		currentLongKey |= baseInd;// same as + baseInd

		currentPalinLongKey >>= _bitPerSymbol;
		currentPalinLongKey |= (complement - baseInd) << preCompute;// (_bitPerSymbol
																	// *
																	// (hashSize
																	// -1));

		long tmpKey = 0;
		for (int j = 0; j < islands.length; j++) {
			tmpKey |= ((currentLongKey & islands[j]) >> islandSize[j]);
		}
		super.setCurrentKey((int) tmpKey);
		// System.out.printf("%10d => %6d  | ",currentLongKey, tmpKey);

		tmpKey = 0;
		for (int j = 0; j < islands.length; j++) {
			tmpKey |= ((currentPalinLongKey & islands[j]) >> islandSize[j]);
		}
		super.setPsuedoPalinKey((int) tmpKey);
		// System.out.printf("%10d => %6d  \n ",currentPalinLongKey, tmpKey);
	}

	public static void main(String[] args) {
		GappedHashtable myhash = new GappedHashtable(args[0], 2);
		byte[] seq = { 0, 1, 0, 1, 0, 2, 3, 1, 0, 2, 3, 0, 1, 3, 2, 0, 1, 3, 2,
				3, 2, 3 };

		for (int i = 0; i < seq.length; i++) {
			// myhash.nextKey((byte)(i % 4));
			// myhash.putCurrentValue(i);
			myhash.nextKey(seq[i]);
			// System.out.println();
		}
	}
}
