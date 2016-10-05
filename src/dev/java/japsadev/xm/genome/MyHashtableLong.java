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
 * On 10 Jul MD added gapped seed as key
 * Note: Maximum 1 << 30 cells as a limit of java array or 15 
 * Maximum =
 */

import japsa.util.IntIterator;
import japsa.xm.hash.PatternStore;

import java.util.Arrays;
import java.util.Random;

public class MyHashtableLong implements PatternStore {
	public Random rnd = new Random(1);
	int hashSize;

	/**
	 * Maximum capacity as of in HashMap class
	 */
	static final int MAXIMUM_CAPACITY = 1 << 30;
	static final int INITIAL_CELL_CAPACITY = 16;

	// private int hashSize;
	private int currentKey = 0;
	private int psuedoPalinKey = 0;

	private int[][] values;
	private int[] valueCount;
	private int _usableBitMask = 0;

	protected int _bitPerSymbol = 2;// For DNA, for protein must be 5
	protected byte complement; // == largest base complement(x) = compliment - x

	private int alloc = 0, used = 0, reallo = 0, init = 0;
	protected int preCompute;

	// A dummy one, should not call this
	protected MyHashtableLong() {

	}

	public MyHashtableLong(int hashSize) {
		this(hashSize, 2);
	}

	public MyHashtableLong(int hashSize, int bpp) {
		initilise(hashSize, bpp);
	}

	public void initilise() {
		initilise(hashSize, _bitPerSymbol);
	}

	public void initilise(int hSize, int bpp) {
		this.hashSize = hSize;
		// Inililise usable bit mask
		_bitPerSymbol = bpp;

		complement = (byte) ((1 << bpp) - 1);

		for (int i = 0; i < hSize * _bitPerSymbol; i++) {
			_usableBitMask <<= 1;
			_usableBitMask |= 1; // same as ++
		}// assert _usableBitMask = 11..11: hashSize * _bitPerSymbol bit

		preCompute = (_bitPerSymbol * (hSize - 1));

		clear();
	}

	public void clear() {
		// Create an array of entries
		valueCount = new int[1 << (_bitPerSymbol * hashSize)];
		values = new int[valueCount.length][];

		// Set all count to -1
		Arrays.fill(valueCount, -1);
	}

	public void reinitialise_optimise() {
		// Create an array of entries
		for (int x = 0; x < valueCount.length; x++) {
			valueCount[x]++;
			if (valueCount[x] > 0) {
				values[x] = new int[valueCount[x]];
				alloc += values[x].length;
				init++;
			}

			valueCount[x] = 0;
		}
		this.currentKey = this.psuedoPalinKey = 0;
	}

	public void initiliseCurrentValue() {
		initiliseValue(currentKey);
	}

	public void initiliseValue(int key) {
		values[key] = new int[INITIAL_CELL_CAPACITY];
		valueCount[key] = 0;

		alloc += INITIAL_CELL_CAPACITY;
		init++;
	}

	/**
	 * Compute the new key as a new base is read in
	 * 
	 * @param baseInd
	 */

	public void nextKey(byte baseInd) {
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
	public void putCurrentValue_psuedo(int val) {
		valueCount[currentKey]++;
	}

	public void printMemoryNeeded() {
		long sum = 0;
		for (int i = 0; i < valueCount.length; i++) {
			if (valueCount[i] >= 0) {
				sum += (valueCount[i] + 1);
			}
		}

		System.out.println(valueCount.length + " arrays of total " + sum);
	}

	/**
	 * Put a number into the hash table at the current key
	 * 
	 * @param val
	 */
	public void putCurrentValue(int val) {
		if (valueCount[currentKey] < 0) {
			initiliseCurrentValue();
		}

		putValue(val, currentKey);
	}

	public void putValue(int val, int key) {
		// invariat valueCount[key] = values[key].length
		if (valueCount[key] == values[key].length) {
			// If no the array is full, reallocate array, doule size
			// int newArray[] = new int[values[key].length << 1];
			int newArray[] = new int[(int) (values[key].length * 1.5)];
			System.arraycopy(values[key], 0, newArray, 0, values[key].length);

			alloc += newArray.length;
			alloc -= values[key].length;

			reallo++;

			values[key] = newArray;
		}
		values[key][valueCount[key]] = val;
		valueCount[key]++;

		used++;
	}

	/**
	 * Return the size of the current hash value
	 * 
	 * @return
	 */
	public int getCurrentCount() {
		return valueCount[currentKey];
	}

	public int getCount(int key) {
		return valueCount[key];
	}

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

	public IntIterator iterator() {
		return new MyHashLongIterator();
	}

	public MyHashLongIterator getLongIterator() {
		return new MyHashLongIterator();
	}

	public IntIterator copyIterator() {
		return new CopyIterator();
	}

	public IntIterator palinIterator() {
		return new PalinIterator();
	}

	public static void main(String[] args) {
		MyHashtableLong myhash = new MyHashtableLong(5, 2);

		for (int i = 0; i < 100; i++) {
			myhash.nextKey((byte) (i % 4));
			myhash.putCurrentValue(i);
			System.out.println();
		}
	}

	public int hashSize() {
		return hashSize;
	}

	class CopyIterator implements IntIterator {
		private int total;
		private int[] myArray = null;

		public CopyIterator() {
			total = valueCount[currentKey];
			if (total < 0)
				total = 0;
			myArray = values[currentKey];
		}

		public boolean hasNext() {
			return total > 0;
		}

		public int sizeAvailable() {
			return total;
		}

		public void remove() {
		}

		public int next() {
			int random_index = rnd.nextInt(total);
			// Shuffle the array

			if (random_index + 1 < total) {// Move the expert selected to end of
											// the list
				// so wont be selected again
				int tmp = myArray[random_index];
				myArray[random_index] = myArray[total - 1];
				myArray[total - 1] = tmp;
			}// else random_index +1 == total

			total--;
			return myArray[total];
		}

		public int next(int random_index) {
			if (random_index + 1 < total) {// Move the expert selected to end of
											// the list
				// so wont be selected again
				int tmp = myArray[random_index];
				myArray[random_index] = myArray[total - 1];
				myArray[total - 1] = tmp;
			}// else random_index +1 == total

			total--;
			return myArray[total];
		}
	}

	class PalinIterator implements IntIterator {
		private int total;
		private int[] myArray = null;

		public PalinIterator() {
			total = valueCount[psuedoPalinKey];
			if (total < 0)
				total = 0;
			myArray = values[psuedoPalinKey];
		}

		public boolean hasNext() {
			return total > 0;
		}

		public int sizeAvailable() {
			return total;
		}

		public int next() {
			// TODO: can reuse the random number generated before
			int random_index = rnd.nextInt(total);
			// Shuffle the array
			if (random_index + 1 < total) {// Move the expert selected to end of
											// the list
				// so wont be selected again
				int tmp = myArray[random_index];
				myArray[random_index] = myArray[total - 1];
				myArray[total - 1] = tmp;
			}// else random_index +1 == total
			total--;

			return (myArray[total]);
		}

		public int next(int random_index) {
			if (random_index + 1 < total) {// Move the expert selected to end of
											// the list
				// so wont be selected again
				int tmp = myArray[random_index];
				myArray[random_index] = myArray[total - 1];
				myArray[total - 1] = tmp;
			}// else random_index +1 == total

			total--;
			return myArray[total];
		}

	}

	class MyHashLongIterator implements IntIterator {
		PalinIterator pIter;
		CopyIterator cIter;
		boolean isPalin = false;

		public MyHashLongIterator() {
			pIter = new PalinIterator();
			cIter = new CopyIterator();
		}

		public boolean isPalin() {
			return isPalin;
		}

		public boolean hasNext() {
			return cIter.hasNext() || pIter.hasNext();
		}

		public int sizeAvailable() {
			return cIter.total + pIter.total;
		}

		public int next() {
			int random_index = rnd.nextInt(sizeAvailable());
			if (random_index < cIter.sizeAvailable()) {
				isPalin = false;
				return cIter.next();
				// return cIter.next(random_index);
			} else {
				isPalin = true;
				return pIter.next();
				// return pIter.next(random_index - cIter.sizeAvailable());
			}
		}
	}

}
