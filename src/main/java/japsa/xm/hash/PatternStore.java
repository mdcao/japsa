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

import japsa.util.IntIterator;

/**
 * @author Minh Duc Cao
 * 
 */
public interface PatternStore {

	/**
	 * Clear every thing in the hash table
	 */
	public void clear();

	/**
	 * Next index key designed from a base
	 * 
	 * @param baseInd
	 *            : a new base
	 */
	public void nextKey(int baseInd);

	/**
	 * Put a value (position) in to current key
	 * 
	 * @param val
	 */
	public void putCurrentValue(int val);

	/**
	 * Return the iterator of all previous candidate copy experts
	 * 
	 * @return
	 */
	public IntIterator copyIterator();

	/**
	 * Return the iterator of all previous candidate palin experts
	 * 
	 * @return
	 */
	public IntIterator palinIterator();

	/**
	 * Return the iterator of all previous candidate copy and palin experts They
	 * are in fact the combination of both copy and palin experts
	 * 
	 * @return
	 */
	public IntIterator iterator();

	// Summary of what is stored
	public void printSummary();

}
