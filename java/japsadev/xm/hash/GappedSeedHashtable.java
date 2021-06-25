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

/**
 * On 10 Jul MD added gapped seed as key
 * Note: Maximum 1 << 30 cells as a limit of java array or 15 
 * Maximum =
 */
import japsa.util.JapsaMath;
import japsa.xm.expert.Expert;

import java.util.StringTokenizer;

public class GappedSeedHashtable {
	/**
	 * Maximum capacity as of in HashMap class
	 */

	static final int INITIAL_CELL_CAPACITY = 16;
	private long currentLongKey = 0;
	private long currentPalinLongKey = 0;

	long currentKey;
	long palinKey;

	// private int _usableBitMask = 0;
	private int _bitPerSymbol = 2;// For DNA, for protein mush be 5

	long[] islands;
	int[] islandSize;
	static final int MAX_STORE = 7000003;// A big prime number
	private HashCell[] cells;
	private int _bitLength;

	public GappedSeedHashtable(String mask) {
		_bitPerSymbol = (int) Math.ceil(JapsaMath
				.log2(Expert.alphabet().size()));
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

		System.out.println(mask);
		for (j = 0; j < islands.length; j++) {
			System.out.printf("%s %d \n", Long.toBinaryString(islands[j]),
					islandSize[j]);
		}

		// Change to count of 1s
		count = mask.length() - count;

		// Inililise usable bit mask
		// for (int i = 0; i < count * _bitPerSymbol; i++) {
		// _usableBitMask <<= 1;
		// _usableBitMask |= 1; // same as ++
		// }// assert _usableBitMask = 11..11: hashSize * _bitPerSymbol bit

		cells = new HashCell[MAX_STORE];
		_bitLength = (mask.length() - 1) * _bitPerSymbol;
	}

	/**
	 * Compute the new key as a new base is read in
	 * 
	 * @param baseInd
	 */

	public void nextKey(int baseInd) {
		// Get the sunsequence under mask
		currentLongKey <<= _bitPerSymbol;
		currentLongKey |= baseInd;// same as + baseInd

		currentKey = 0;
		for (int j = 0; j < islands.length; j++) {
			currentKey |= ((currentLongKey & islands[j]) >> islandSize[j]);
		}

	}

	public void nextPalinKey(long baseInd) {
		baseInd <<= (_bitLength);
		currentPalinLongKey >>= _bitPerSymbol;
		currentPalinLongKey |= baseInd;// same as + baseInd

		palinKey = 0;
		for (int j = 0; j < islands.length; j++) {
			palinKey |= ((currentPalinLongKey & islands[j]) >> islandSize[j]);
		}
	}

	// Get and create
	public HashCell getCell(long key) {
		int index = (int) (key % MAX_STORE);

		HashCell cell = cells[index];
		if (cell == null) {
			cells[index] = new HashCell(index);
			return cells[index];
		}
		while (cell.id != key && cell.next != null) {
			cell = cell.next;
		}
		// Iether cell.next == null or cell.id = key
		if (cell.id == key)
			return cell;
		// else
		cell.next = new HashCell(key);
		return cell.next;
	}

	public HashCell getCurrentCell() {
		return getCell(currentKey);
	}

	public HashCell getPalinCell() {
		return getCell(palinKey);
	}

	public void putCurrentCell(int val) {
		HashCell cell = getCell(currentKey);
		cell.putValue(val);
	}

	public class HashCell {
		long id;
		HashCell next = null;
		int count = -1;
		int[] values = null;

		HashCell(long id) {
			values = new int[INITIAL_CELL_CAPACITY];
			this.id = id;
		}

		/**
		 * Put a number into the hash table at the current key
		 * 
		 * @param val
		 */
		public void putValue(int val) {
			if (count < 0) {
				values = new int[INITIAL_CELL_CAPACITY];
				count = 0;
			} else if (count == values.length) {// Full
				int newArray[] = new int[values.length << 1];
				System.arraycopy(values, 0, newArray, 0, values.length);
				values = newArray;
			}
			values[count] = val;
			count++;
		}

		public int getCount() {
			return count;
		}

		public int[] getValues() {
			return values;
		}
	}
}
