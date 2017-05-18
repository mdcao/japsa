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
 * 04/01/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.seq;

import japsa.seq.Alphabet;
import japsa.util.ByteArray;

import java.util.Arrays;
import java.util.Random;




/**
 * Implement a sequence in which the data are stored in a byte array, a byte for
 * a element of the sequence. The byte array object is not changleble, though the
 * value in each element (base) can be changed.
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com) * 
 */
public class Sequence extends AbstractSequence implements Cloneable {	
	/**
	 * The array to hold the sequence
	 */
	private final byte[] byteSeq;

	/**
	 * Create an empty sequence with a specified length
	 * 
	 * @param alphabet
	 * 
	 */
	public Sequence(Alphabet alphabet, int length) {
		super(alphabet);
		byteSeq = new byte[length];
	}

	public Sequence(Alphabet alphabet, int length, String name) {
		super(alphabet, name);
		byteSeq = new byte[length];
	}

	/**
	 * Construct a sequence from a sequence of characters.
	 * @param alphabet
	 * @param charSeq
	 * @param name
	 */	
	public Sequence(Alphabet alphabet, char[] charSeq, String name) {
		super(alphabet, name);
		byteSeq = new byte[charSeq.length];
		for (int i = 0; i < byteSeq.length; i++) {
			byteSeq[i] = (byte) alphabet.char2int(charSeq[i]);
		}
	}
	
	/**
	 * Construct a sequence from a ByteArray object.
	 * @param alphabet
	 * @param bArray
	 * @param name
	 */
	public Sequence(Alphabet alphabet, ByteArray bArray, String name) {
		super(alphabet, name);		
		byteSeq = bArray.toArray();
	}


	/**
	 * Copy the byte array up to the length
	 * 
	 * @param alphabet
	 * @param byteArray
	 * @param length
	 */

	public Sequence(Alphabet alphabet, byte[] byteArray, int length) {
		super(alphabet);
		byteSeq = Arrays.copyOf(byteArray, length);
	}

	public Sequence(Alphabet alphabet, byte[] byteArray, int length, String name) {
		super(alphabet, name);
		byteSeq = Arrays.copyOf(byteArray, length);
	}

	public Sequence(Alphabet alphabet, byte[] byteArray) {
		this(alphabet, byteArray, byteArray.length);
	}

	public Sequence(Alphabet alphabet, byte[] byteArray, String name) {
		this(alphabet, byteArray, byteArray.length, name);
	}
	
	//public Sequence(Alphabet dna, byte[] byteArray, String name, String desc) {
	//	this(dna, byteArray, byteArray.length, name);
	//	setDesc(desc);
	//}
	
	/**
	 * Construct a sequence with an dna from the string represent the 
	 * sequence and the name
	 * @param alphabet
	 * @param seqString
	 * @param name
	 */
	
	public Sequence(Alphabet alphabet, String seqString, String name) {
		//this(dna, seqStr.toCharArray(),name);		
		this(alphabet,seqString.length(),name);	
		
		for (int i = 0; i < byteSeq.length; i++) {
			byteSeq[i] = (byte) alphabet.char2int(seqString.charAt(i));
		}
	}
	

	/**
	 * Return the length of the sequence
	 * 
	 * @see japsa.seq.AbstractSequence#length()
	 */
	@Override
	public int length() {
		return byteSeq.length;
	}

	/*
	 * @see bio.seq.AbstractSequence#symbolAt(int)
	 */
	@Override
	public int symbolAt(int loc) {
		return byteSeq[loc];// This should convert into an int
	}

	public byte getBase(int loc) {
		return byteSeq[loc];
	}

	/**
	 * Concatenate a sequence with another
	 * @param anotherSeq
	 */

	public Sequence concatenate(Sequence anotherSeq) {				
		if (this.alphabet() != anotherSeq.alphabet())
			throw new RuntimeException("The alphabets do not match");
		
		Sequence newSeq = new Sequence(alphabet(), length() + anotherSeq.length());
		
		for (int j = 0; j < anotherSeq.length(); j++) {
			newSeq.byteSeq[byteSeq.length + j] = anotherSeq.byteSeq[j];
		}

		return newSeq;
	}
	
	public byte[] toBytes() {
		return byteSeq;
	}

	/**
	 * @param loc
	 * @param base
	 */
	public byte setBase(int loc, byte base) {
		byteSeq[loc] = base;
		return base;
	}
	
	public void setSymbol(int loc, int symbol){
		byteSeq[loc] = (byte) symbol;
	}
		
	/**
	 * Clone the sequence
	 */
	public Sequence clone(){
		Sequence seq = new Sequence(alphabet(), byteSeq, getName());
		seq.setDesc(getDesc());
		return seq;
	}
	
	/**
	 * Create a random sequence with some length and some frequency distribution
	 * @param alphabet
	 * @param length
	 * @param freqs
	 * @param rand a random generator
	 * @return
	 */
	public static Sequence random(Alphabet alphabet, int length, double [] freqs, Random rand){
		if (freqs.length != alphabet.size()){
			throw new RuntimeException("Frequencies array should have the same size as dna");
		}
		Sequence seq = new Sequence(alphabet, length);
		
		//setting accumulation distribution
		double [] accum = new double[freqs.length];		
		accum[0] = freqs[0];
		for (int i = 1; i < freqs.length;i++){
			accum[i] = accum[i-1] + freqs[i];
		}
		
		//normalise, just in case
		double sum = accum[accum.length-1];
		for (int i = 0; i < freqs.length;i++){
			accum[i] /= sum;
		}		
		//java.util.Random rand = new java.util.Random();		
		
		for (int x = 0; x < seq.length();x++){
			double r = rand.nextDouble();
			for (byte i = 0; i < accum.length; i++){
				if (accum[i] >= r){
					seq.byteSeq[x] = i;
					break;
				}
			}
		}
		
		return seq;
	}
	/**
	 * Create a random sequence
	 * @param alphabet
	 * @param length
	 * @param freqs
	 * @return
	 */
	public static Sequence random(Alphabet alphabet, int length, double [] freqs){
		return random(alphabet, length, freqs, new Random());
	}

	/* (non-Javadoc)
	 * @see java.lang.CharSequence#subSequence(int, int)
	 */
	@Override
	public Sequence subSequence(int start, int end) {
		byte [] newSeq = Arrays.copyOfRange(byteSeq, start, end); 
		return new Sequence(alphabet(), newSeq);
	}
	/**
	 * Start: 0-index inclusive,
	 * end:0-index, exclusive
	 * @param start
	 * @param end
	 * @param flank
	 * @return
	 */
	public Sequence subsequenceWithFlank(int start, int end, int flank){
		Sequence featureSeq = new Sequence(alphabet(), end - start + 2 * flank);

		for (int i = 0; i < featureSeq.length();i++){
			int j = start - flank + i; 
			if (j < 0 || j>= length())
				featureSeq.setSymbol(i, Alphabet.DNA.N);
			else 
				featureSeq.setSymbol(i, getBase(j));
		}//for		
		return featureSeq;
	}

}
