/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
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
 * 3. Neither the names of the institutions nor the names of the contributors*
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

/*                           Revision History                                
 * 11/01/2012 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsa.util;

import java.util.Arrays;

/**
 * A re-implementation of BitSet class with cheaper operations by bypassing
 * safety checks
 * 
 * @author Minh Duc Cao
 * 
 */
public class MyBitSet {

	public static final int SIZE_INT = 32;
	public static int[] POS = new int[SIZE_INT];
	static {
		POS[0] = 1;
		for (int i = 1; i < SIZE_INT; i++) {
			POS[i] = POS[i - 1] << 1;
		}
	}

	private int array[];
	
	/**
	 * Create a bitset containing at least length bits
	 * @param length
	 */

	public MyBitSet(int length) {
		// System.out.println("Alocating " + (length / SIZE_INT + 1) +
		// " ints ");
		array = new int[length / SIZE_INT + 1];
		//Arrays.fill(array, 0);
	}

	/**
	 * Set bit at position pos
	 * @param pos
	 */
	public void set(int pos) {
		int ind = pos / SIZE_INT, place = pos % SIZE_INT;
		array[ind] |= POS[place];
	}
	
	/**
	 * Clear bit at position pos
	 * @param pos
	 */

	public void clear(int pos) {
		int ind = pos / SIZE_INT, place = pos % SIZE_INT;
		array[ind] &= ~POS[place];
	}

	/**
	 * Get bit at position pos
	 * @param pos
	 * @return
	 */
	public boolean get(int pos) {
		int ind = pos / SIZE_INT, place = pos % SIZE_INT;
		return ((array[ind] & POS[place]) != 0);
	}

	public boolean equal(MyBitSet other) {
		return Arrays.equals(array, other.array);
	}

	// Other must not shorter than me
	public void copy(MyBitSet other) {
		
		for (int i = 0; i < array.length; i++)
			other.array[i] = array[i];
	}	
}
