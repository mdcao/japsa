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
 *
 */

import japsa.util.IntIterator;

import java.util.Random;

public class MyBigHashtableLong {// implements PatternStore{

	public Random rnd = new Random(1);
	int hashSize;
	int numSeqs;

	/**
	 * Maximum capacity as of in HashMap class
	 */
	static final int MAXIMUM_CAPACITY = 1 << 30;
	static final int INITIAL_CELL_CAPACITY = 16;

	// private int hashSize;
	private int currentKey = 0;
	private int psuedoPalinKey = 0;

	// values[key][seqID][posSrc]
	private int[][][] values;

	// valueCount[key][seqID]
	private int[][] valueCount;
	private int _usableBitMask = 0;

	protected int _bitPerSymbol = 2;// For DNA, for protein must be 5
	protected int complement; // == largest base complement(x) = complement - x.
								// In case of dna, complement = 3

	private long alloc = 0, used = 0, reallo = 0, init = 0;
	protected int preCompute;
	private MyBigHashLongIterator copyIter = null;

	// A dummy one, should not call this
	protected MyBigHashtableLong() {

	}

	public MyBigHashtableLong(int nSeqs, int hashSize) {
		this(nSeqs, hashSize, 2);
	}

	public MyBigHashtableLong(int nSeqs, int hashSize, int bpp) {
		this.numSeqs = nSeqs;
		copyIter = new MyBigHashLongIterator();
		initilise(hashSize, bpp);
	}

	public void initilise() {
		initilise(hashSize, _bitPerSymbol);
	}

	public void initilise(int hSize, int bpp) {
		this.hashSize = hSize;
		// Inililise usable bit mask
		_bitPerSymbol = bpp;

		complement = (byte) ((1 << bpp) - 1);// 1 << 2 - 1 = 3

		for (int i = 0; i < hashSize * _bitPerSymbol; i++) {
			_usableBitMask <<= 1;
			_usableBitMask |= 1; // same as ++
		}// assert _usableBitMask = 11..11: hashSize * _bitPerSymbol bit

		preCompute = (_bitPerSymbol * (hashSize - 1));

		clear();
	}

	public void clear() {
		// Create an array of entries
		long s1 = System.currentTimeMillis();
		valueCount = new int[1 << (_bitPerSymbol * hashSize)][numSeqs];
		values = new int[valueCount.length][numSeqs][];
		long s2 = System.currentTimeMillis();

		System.out.println((s2 - s1) + "=============");

		// Set all count to -1
		for (int x = 0; x < valueCount.length; x++) {
			for (int y = 0; y < numSeqs; y++)
				valueCount[x][y] = -1;
		}
	}

	public void reinitialise_optimise() {
		// Create an array of entries
		for (int x = 0; x < valueCount.length; x++) {
			for (int y = 0; y < numSeqs; y++) {
				valueCount[x][y]++;
				if (valueCount[x][y] > 0) {
					values[x][y] = new int[valueCount[x][y]];
					alloc += values[x][y].length;
					init++;
				}
				valueCount[x][y] = 0;
			}
		}
		this.currentKey = this.psuedoPalinKey = 0;
	}

	public void initiliseCurrentValue() {
		for (int y = 0; y < numSeqs; y++)
			initiliseValue(currentKey, y);
	}

	/**
	 * Only call when the key has not been initilised
	 * 
	 * @param key
	 * @param sid
	 */
	public void initiliseValue(int key, int sid) {
		values[key][sid] = new int[INITIAL_CELL_CAPACITY];
		valueCount[key][sid] = 0;

		alloc += INITIAL_CELL_CAPACITY;
		init++;
	}

	/**
	 * Compute the new key as a new base is read in
	 * 
	 * @param baseInd
	 */

	public void nextKey(int baseInd) {
		// currentKey = currentKey x
		currentKey <<= _bitPerSymbol;
		currentKey &= _usableBitMask;
		currentKey |= baseInd;// same as + baseInd

		psuedoPalinKey >>= _bitPerSymbol;
		// baseInd = ;
		psuedoPalinKey |= (complement - baseInd) << preCompute;// (_bitPerSymbol
																// * (hashSize
																// -1));

	}

	/**
	 * Put a number into the hash table at the current key
	 * 
	 * @param val
	 */
	public void putCurrentValue_psuedo(int sid, int val) {
		valueCount[currentKey][sid]++;
	}

	public void printMemoryNeeded() {
		long sum = 0;
		for (int i = 0; i < valueCount.length; i++) {
			for (int y = 0; y < numSeqs; y++)
				if (valueCount[i][y] >= 0) {
					sum += (valueCount[i][y] + 1);
				}
		}
		System.out.println(valueCount.length + " arrays of total " + sum);
	}

	/**
	 * Put a number into the hash table at the current key
	 * 
	 * @param val
	 */
	public void putCurrentValue(int sid, int val) {
		if (valueCount[currentKey][sid] < 0) {
			initiliseValue(currentKey, sid);
			// initiliseCurrentValue();
		}

		putValue(val, currentKey, sid);
	}

	public void putValue(int val, int key, int sid) {
		// invariat valueCount[key] = values[key].length
		if (valueCount[key][sid] == values[key][sid].length) {
			// If no the array is full, reallocate array, doule size
			// int newArray[] = new int[values[key].length << 1];
			int newArray[] = new int[(int) (values[key][sid].length * 1.5)];
			System.arraycopy(values[key][sid], 0, newArray, 0,
					values[key][sid].length);

			alloc += newArray.length;
			alloc -= values[key].length;

			reallo++;

			values[key][sid] = newArray;
		}
		values[key][sid][valueCount[key][sid]] = val;
		valueCount[key][sid]++;

		used++;
	}

	// ///////////////////////////////////////////
	/**
	 * Return the size of the current hash value
	 * 
	 * @return
	 */
	// public int getCurrentCount(){
	// return valueCount[currentKey];
	// }

	// public int getCount(int key){
	// return valueCount[key];
	// }
	/**
	 * Return the current key
	 * 
	 * @return
	 */
	public int getCurrentKey() {
		return currentKey;
	}

	public void setCurrentKey(int key) {
		this.currentKey = key;
	}

	public int getPsuedoPalinKey() {
		return psuedoPalinKey;
	}

	public void setPsuedoPalinKey(int psuedoPalinKey) {
		this.psuedoPalinKey = psuedoPalinKey;
	}

	public void printSummary() {
		System.out
				.printf(" Init %d(%d), reallocation %d, Allocatate %d, used %d, (%.2f %% )\n",
						init, valueCount.length, reallo, alloc, used,
						(used * 100.0 / alloc));
	}

	public MyBigHashLongIterator iterator() {
		copyIter.reset();
		return copyIter;
	}

	public static void main(String[] args) {
		// MyBigHashtableLong myhash = new MyBigHashtableLong(5,2);

		// for ( int i = 0; i < 100; i++){
		// myhash.nextKey((byte)(i % 4));
		// myhash.putCurrentValue(i);
		// System.out.println();
		// }
	}

	public int hashSize() {
		return hashSize;
	}

	class MyBigHashLongIterator implements IntIterator {

		// v[sid][posSrc]
		int[][] copies;// = new int[2 * numSeqs][];
		int[] totals;// totals for all

		int total;// total number of posSrc

		boolean isPalin = false;
		int sid;

		public MyBigHashLongIterator() {
			copies = new int[2 * numSeqs][];
			totals = new int[2 * numSeqs];
		}

		/**
		 * Set an array of arrays
		 */
		public void reset() {
			total = 0;
			for (int y = 0; y < numSeqs; y++) {
				copies[y] = values[currentKey][y];
				totals[y] = valueCount[currentKey][y];

				if (totals[y] < 0)
					totals[y] = 0;

				total += totals[y];

				copies[y + numSeqs] = values[psuedoPalinKey][y];
				totals[y + numSeqs] = valueCount[psuedoPalinKey][y];
				if (totals[y + numSeqs] < 0)
					totals[y + numSeqs] = 0;
				total += totals[y + numSeqs];
			}
		}

		public boolean isPalin() {
			return isPalin;
		}

		public int sequenceID() {
			return 0;
		}

		public boolean hasNext() {
			return total > 0;
			// return cIter.hasNext() || pIter.hasNext();
		}

		public int sizeAvailable() {
			return total;
			// return cIter.total + pIter.total;
		}

		public int next() {
			int random_index = rnd.nextInt(total);

			sid = 0;
			while (random_index >= totals[sid]) {
				random_index -= totals[sid];
				sid++;
			}

			// assert: accTotal < randomIndex

			int[] myArray = copies[sid];

			totals[sid]--;
			int myTotal = totals[sid];// point to the last availble

			if (random_index < myTotal) {// Move the expert selected to end of
											// the list
				// so wont be selected again
				int tmp = myArray[random_index];
				myArray[random_index] = myArray[myTotal];
				myArray[myTotal] = tmp;
			}// else random_index +1 == total

			if (sid < numSeqs) {
				isPalin = false;
			} else {
				isPalin = true;
				sid -= numSeqs;
			}

			total--;
			return myArray[myTotal];

		}
	}

}
