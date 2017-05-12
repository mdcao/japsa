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
 * 25/03/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.seq;

import java.io.IOException;
import java.util.Arrays;




/**
 * A fastq record is basically a sequence with a quality string, represented as
 * an array of bytes.
 * 
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
public class FastqSequence extends Sequence {
	
	/**
	 * Array of base quality
	 */
	private final byte [] quality;	
	
	/**
	 * Create a fastq record
	 * @param dna
	 * @param byteArray
	 * @param quality
	 * @param name
	 */
	public FastqSequence(Alphabet alphabet, byte[] byteArray, byte[] quality, String name) {
		super(alphabet, byteArray, name);
		this.quality = Arrays.copyOf(quality, length());
	}
		
	/**
	 * Create a fastq record with a fix length
	 * @param dna
	 * @param byteArray
	 * @param quality
	 * @param length
	 * @param name
	 */
	
	public FastqSequence(Alphabet alphabet, byte[] byteArray, byte[] quality, int length,  String name) {
		super(alphabet, byteArray, length, name);
		this.quality = Arrays.copyOf(quality, length);		
	}
	/**
	 * Create a fastq from 4 strings: name, seq, +, and quality
	 * @param dna
	 * @param toks
	 */
	public FastqSequence(Alphabet alphabet, String [] toks) {	
		super(alphabet, toks[1], toks[0]);
		this.quality = toks[3].getBytes();		
	}

	
	
	/* (non-Javadoc)
	 * @see japsa.seq.AbstractSequence#write(japsa.seq.SequenceOutputStream)
	 * 
	 * Write in fastq format to output stream
	 */
	public void print(SequenceOutputStream out) throws IOException{		
		out.print('@');
		out.print(getName());
		out.print('\n');
		
		for (int i = 0; i < length();i++)
			out.print(charAt(i));
		out.print('\n');
		
		out.print('+');		
		out.print('\n');				
		
		out.write(quality);
		out.print('\n');		
	}
	
	public byte getQualByte(int loc){
		return quality[loc];
	}
	
	/**
	 * Convert to a string. This method is only for convenient and is very slow
	 */
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append('@');
		sb.append(getName());
		sb.append('\n');
		
		for (int i = 0; i < length();i++)
			sb.append(charAt(i));
		sb.append("\n+\n");
						
		
		for (int i = 0; i < quality.length; i++)
			sb.append(quality[i]);
		//out.print('\n');		
		
		return sb.toString();
	}
}
