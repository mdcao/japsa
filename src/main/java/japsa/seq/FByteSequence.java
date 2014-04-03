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

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;




/**
 * Implement a sequence in which the data are stored in a byte array, a byte for
 * a element of the sequence. The byte array object is not changleble, though the
 * value in each element (base) can be changed.
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com) * 
 */
public class FByteSequence extends AbstractSequence {	
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
	public FByteSequence(Alphabet alphabet, int length) {
		super(alphabet);
		byteSeq = new byte[length];
	}

	public FByteSequence(Alphabet alphabet, int length, String name) {
		super(alphabet, name);
		byteSeq = new byte[length];
	}

	/**
	 * Construct a sequence from a sequence of characters.
	 * @param alphabet
	 * @param charSeq
	 * @param name
	 */
	@Deprecated
	public FByteSequence(Alphabet alphabet, char[] charSeq, String name) {
		super(alphabet, name);
		byteSeq = new byte[charSeq.length];
		for (int i = 0; i < byteSeq.length; i++) {
			byteSeq[i] = (byte) alphabet.char2int(charSeq[i]);
		}
	}

	@Deprecated
	public FByteSequence(Alphabet alphabet, char[] charSeq) {
		this(alphabet, charSeq, "");
	}

	/**
	 * Copy the byte array up to the length
	 * 
	 * @param alphabet
	 * @param byteArray
	 * @param length
	 */

	public FByteSequence(Alphabet alphabet, byte[] byteArray, int length) {
		super(alphabet);
		byteSeq = Arrays.copyOf(byteArray, length);
	}

	public FByteSequence(Alphabet alphabet, byte[] byteArray, int length, String name) {
		super(alphabet, name);
		byteSeq = Arrays.copyOf(byteArray, length);
	}

	public FByteSequence(Alphabet alphabet, byte[] byteArray) {
		this(alphabet, byteArray, byteArray.length);
	}

	public FByteSequence(Alphabet alphabet, byte[] byteArray, String name) {
		this(alphabet, byteArray, byteArray.length, name);
	}
	
	public FByteSequence(Alphabet alphabet, byte[] byteArray, String name, String desc) {
		this(alphabet, byteArray, byteArray.length, name);
		setDesc(desc);
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
		if ( byteSeq[loc] >=0)
			return byteSeq[loc];
		else
			return 256 + byteSeq[loc];
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see bio.seq.AbstractSequence#toBytes()
	 */
	public byte[] toBytes() {
		return byteSeq;
	}

	/**
	 * @param loc
	 * @param base
	 */
	public void setBase(int loc, int base) {
		byteSeq[loc] = (byte) base;
	}
	
	public void setSymbol(int loc, int symbol){
		byteSeq[loc] = (byte) symbol;
	}

		
	/**
	 * Clone the sequence
	 */
	public FByteSequence clone(){
		return new FByteSequence(alphabet(), byteSeq, getName(), getDesc());
	}
	
	/**
	 * Create a random sequence with some length and some frequency distribution
	 * @param alphabet
	 * @param length
	 * @param freqs
	 * @param rand a random generator
	 * @return
	 */
	public static FByteSequence random(Alphabet alphabet, int length, double [] freqs, Random rand){
		if (freqs.length != alphabet.size()){
			throw new RuntimeException("Frequencies array should have the same size as alphabet");
		}
		FByteSequence seq = new FByteSequence(alphabet, length);
		
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
	
	public static FByteSequence read(String fileName) throws IOException{
		BufferedInputStream  in = new BufferedInputStream (new  FileInputStream(fileName));
		byte [] seqByte = new byte[8192];
		int	seqIndex = 0;
		
		while(true){
			int n = in.read();
			if (n < 0)
				break;
			
			if (seqIndex >= seqByte.length) {// Full
				int newLength = seqByte.length * 2;
				if (newLength < 0) {
					newLength = Integer.MAX_VALUE;
				}
				// if the array is extended
				if (newLength <= seqIndex) {
					in.close();
					throw new RuntimeException(
							"Sequence is too long to handle");				}
				seqByte = Arrays.copyOf(seqByte, newLength);
			}
			//next symbol
			seqByte[seqIndex++] = (byte)n;			
		}
		in.close();
		return new FByteSequence(new Alphabet.AllByte(), seqByte, seqIndex, "ad");		
	}
	
	public static void main(String[] args) throws IOException{
/*********************************************		
		FByteSequence seq = FByteSequence.read("../data/enwik8");
		System.out.println(seq.length());
		
		for (int i = 0; i < 1000; i ++){
			System.out.print(seq.charAt(i));;
		}
		
		ExpertModel eModel = new ExpertModel(2, seq.alphabet() , 20, 0, 0.15, 3, false);
		
		eModel.setHashType("hash");
		eModel.setSelfRep(true);
		
		eModel.offSetSeed = new RepeatCountExpert(new Sequence(null,0), 0, null,
				RepeatExpert.COPY_TYPE);
		eModel.palinSeed = new RepeatCountExpert(new Sequence(null,0), 0, null,
				RepeatExpert.PALIN_TYPE);
		
		FByteSequence [] seqs = new FByteSequence[1];
		seqs[0] = seq;
		System.out.println(eModel.encode1(seqs));
/*********************************************/		
		
	}
	
}
