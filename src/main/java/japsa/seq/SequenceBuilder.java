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
 * 04/01/2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.seq;

import japsa.seq.Alphabet;
import japsa.util.ByteArray;

import java.util.Arrays;




/**
 * Implement sequence class that size can be changed by adding/removing bases
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com) * 
 */
public class SequenceBuilder extends AbstractSequence {	
	/**
	 * The array to hold the sequence
	 */
	private byte[] byteSeq;
	int length = 0;

	/**
	 * Create an empty sequence with a maximum length
	 * 
	 * @param dna
	 * 
	 */
	public SequenceBuilder(Alphabet alphabet, int maxLength) {
		super(alphabet);
		byteSeq = new byte[maxLength];	
	}

	public SequenceBuilder(Alphabet alphabet, int maxLength, String name) {
		super(alphabet, name);
		byteSeq = new byte[maxLength];
	}

	/**
	 * Construct a sequence from a ByteArray object.
	 * @param alphabet
	 * @param bArray
	 * @param name
	 */
	public SequenceBuilder(Alphabet alphabet, ByteArray bArray, String name) {
		super(alphabet, name);		
		byteSeq = bArray.toArray();
		length = byteSeq.length;
	}
	
	public SequenceBuilder(Sequence seq, int length){
		this (seq.alphabet(), length);
		append(seq);
	}


	/**
	 * Copy the byte array up to the length
	 * 
	 * @param dna
	 * @param byteArray
	 * @param length
	 */

	public SequenceBuilder(Alphabet alphabet, byte[] byteArray, int length) {
		super(alphabet);
		byteSeq = Arrays.copyOf(byteArray, length);
		this.length = length;
	}

	public SequenceBuilder(Alphabet alphabet, byte[] byteArray, int length, String name) {
		super(alphabet, name);
		byteSeq = Arrays.copyOf(byteArray, length);
		this.length = length;
	}

	public SequenceBuilder(Alphabet alphabet, byte[] byteArray) {
		this(alphabet, byteArray, byteArray.length);
	}

	public SequenceBuilder(Alphabet alphabet, byte[] byteArray, String name) {
		this(alphabet, byteArray, byteArray.length, name);
	}

	/**
	 * Construct a sequence with an dna from the string represent the 
	 * sequence and the name
	 * @param dna
	 * @param seqStr
	 * @param name
	 */

	public SequenceBuilder(Alphabet alphabet, String seqString, String name) {
		this(alphabet,seqString.length(),name);	

		for (int i = 0; i < byteSeq.length; i++) {
			append((byte) alphabet.char2int(seqString.charAt(i)));
		}
	}

	/**
	 * Append a base to the end of the sequence
	 * @param base
	 */
	public void append(byte base){
		if (length >= byteSeq.length){
			//extend
			int newLength = byteSeq.length * 2;
			if (newLength < 0) {
				newLength = Integer.MAX_VALUE;
			}
			// if the array is extended
			if (newLength <= length) {
				// in.close();
				throw new RuntimeException(
					"Sequence is too long to handle");
			}
			// create new array of byte
			byteSeq = Arrays.copyOf(byteSeq, newLength);
		}
		byteSeq[length] = base;
		length ++;
	}

	/**
	 * append another sequence from start (inclusive) to end (exclusive)
	 * @param seq
	 * @param start
	 * @param end
	 */
	public void append(AbstractSequence seq, int start, int end){
		if (length + end - start > byteSeq.length){
			//extend			
			int newLength = byteSeq.length * 2;
			while (true){
				newLength = newLength * 2;
				if (newLength < 0) {
					newLength = Integer.MAX_VALUE;
					byteSeq = Arrays.copyOf(byteSeq, newLength);
					break;
				}
				
				//TODO: if the array is extended
				//if (newLength < length + end - start) {
				//	// in.close();
				//	throw new RuntimeException(
				//		"Sequence is too long to handle");
				//}
				// create new array of byte
				if (length + end - start <= newLength){
					byteSeq = Arrays.copyOf(byteSeq, newLength);
					break;
				}
			}
		}
		for (int i = start; i < end; i++)
			byteSeq[length++] = (byte) seq.symbolAt(i);		
	}

	public void append(AbstractSequence seq){
		append(seq,0,seq.length());		
	}

	/**
	 * Return the length of the sequence
	 * 
	 * @see japsa.seq.AbstractSequence#length()
	 */
	@Override
	public int length() {
		return length;
	}

	public Sequence toSequence(){
		return new Sequence(alphabet(), byteSeq, length, getName());
	}

	public void reset(){
		length = 0;
	}

	/*
	 * @see bio.seq.AbstractSequence#symbolAt(int)
	 */
	@Override
	public int symbolAt(int loc) {
		if (loc < 0 || loc >= length())
			return -1;

		return byteSeq[loc];// This should convert into an int
	}

	public byte getBase(int loc){
		if (loc < 0 || loc >= length())
			return -1;

		return byteSeq[loc];
	}

	/**
	 * @param loc
	 * @param base
	 */
	public byte setBase(int loc, byte base) {
		if (loc < 0 || loc >= length()){
			throw new RuntimeException("Wrong location (max " + length + "):" + loc);
		}
		byteSeq[loc] = base;
		return base;
	}

	public void setSymbol(int loc, int symbol){
		if (loc < 0 || loc >= length())
			throw new RuntimeException("Wrong location (max " + length + "):" + loc);

		byteSeq[loc] = (byte) symbol;
	}

	/* (non-Javadoc)
	 * @see java.lang.CharSequence#subSequence(int, int)
	 */
	@Override
	public SequenceBuilder subSequence(int start, int end) {
		byte [] newSeq = Arrays.copyOfRange(byteSeq, start, end); 
		return new SequenceBuilder(alphabet(), newSeq);
	}
}
