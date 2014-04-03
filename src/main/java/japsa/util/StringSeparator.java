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

/**************************     REVISION HISTORY    **************************
 * File: StringSeparator.java
 * 14/11/2013 - Minh Duc Cao: Created
 *
 ****************************************************************************/

package japsa.util;

import java.util.Iterator;


/**
 * Separate a string to substrings with charater separator (such as space, 
 * tab, etc). This class is designed to overcome the inefficiency of 
 * String.split() and StringTokenizer.
 * 
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 */
public class StringSeparator implements Iterator<String>{
	/**
	 * The string to separate
	 */
	private String input;
	
	/**
	 * The delimiter
	 */
	private char delimiter;
	/**
	 * maximum number of tokens
	 */
	private int limit;
	
	// The next are for internal processing
	private int count = 1;
	private int currentIndex = 0, //it should point to the next character after the last delimiter 
			     nextIndex; //point to the next delimiter
	/**
	 * Create a new separator with at most limit number of tokens
	 * 
	 * @param input
	 * @param delimiter
	 * @param limit
	 */
	public StringSeparator(String input, char delimiter, int limit){
		this.input = input;
		this.delimiter = delimiter;	
		this.limit = limit;
		
		seekNext();//Prepare for the first token		
	}
	
	public StringSeparator(String input, char delimiter){
		this(input, delimiter, 0);
	}
	
	private void seekNext(){
		if (count == limit){
			
		}
		nextIndex = currentIndex;
		while (nextIndex < input.length()){
			if (input.charAt(nextIndex) == delimiter)
				break;
			
			nextIndex ++;
		}		
	}
	
	/**
	 * Reset the Separator with a new delimiter
	 * @param newDelim
	 */
	public void reset(char newDelim){
		delimiter = newDelim;
		reset();
	}
	/**
	 * Reset the separator
	 */
	public void reset(){
		currentIndex = 0;
		count = 1;
		seekNext();		
	}	

	/* (non-Javadoc)
	 * @see java.util.Iterator#hasNext()
	 */
	@Override
	public boolean hasNext() {		
		return nextIndex <= input.length();
	}

	/* (non-Javadoc)
	 * @see java.util.Iterator#next()
	 */
	@Override
	public String next() {
		//construct the current token
		String ret = input.substring(currentIndex, nextIndex);		
		
		//Move to the next token
		currentIndex = nextIndex + 1;
		seekNext();		
		
		return ret;
	}

	/* (non-Javadoc)
	 * @see java.util.Iterator#remove()
	 */
	@Override
	public void remove() {
		//Doing nothing		
	}
}
