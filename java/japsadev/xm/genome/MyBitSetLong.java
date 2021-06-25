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

import java.io.Serializable;
import java.util.Arrays;

/**
 * @author minhduc
 * 
 */
public class MyBitSetLong implements Serializable {
	private static final long serialVersionUID = 1L;

	public static final int SIZE_INT = 32;
	public static final int LOG_INT_SIZE = 5;
	public static final long FIVE_1s = 31;// = 0x11111

	public static int[] POS = new int[SIZE_INT];
	static {
		POS[0] = 1;
		for (int i = 1; i < SIZE_INT; i++) {
			POS[i] = POS[i - 1] << 1;
		}
	}

	int array[];

	public MyBitSetLong(long length) {
		// System.out.println("Alocating " + (length / SIZE_INT + 1) +
		// " ints ");
		int array_length = (int) (length / SIZE_INT + 1);
		array = new int[array_length];
		Arrays.fill(array, 0);
	}

	public void set(long pos) {
		// int ind = (int) (posSrc / SIZE_INT), place = (int) (posSrc % SIZE_INT) ;
		int ind = (int) (pos >> LOG_INT_SIZE), place = (int) (pos & FIVE_1s);
		array[ind] |= POS[place];
	}

	public void clear(long pos) {
		// int ind = (int) (posSrc / SIZE_INT), place = (int) (posSrc % SIZE_INT) ;
		int ind = (int) (pos >> LOG_INT_SIZE), place = (int) (pos & FIVE_1s);
		array[ind] &= ~POS[place];
	}

	public boolean get(long pos) {
		// int ind = (int) (posSrc / SIZE_INT), place = (int) (posSrc % SIZE_INT);
		int ind = (int) (pos >> LOG_INT_SIZE), place = (int) (pos & FIVE_1s);
		return ((array[ind] & POS[place]) != 0);
	}

	public boolean equal(MyBitSetLong other) {
		if (array.length != other.array.length)
			return false;
		for (int i = 0; i < array.length; i++) {
			if (array[i] != other.array[i])
				return false;
		}

		return true;
	}

	// Other must not shorter than me
	public void copy(MyBitSetLong other) {
		for (int i = 0; i < array.length; i++)
			other.array[i] = array[i];
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		MyBitSetLong aBitSet = new MyBitSetLong(4000000000l);
		for (long i = 4000000000l - 100; i < 4000000000l; i++) {
			aBitSet.set(i);
		}

	}

}
