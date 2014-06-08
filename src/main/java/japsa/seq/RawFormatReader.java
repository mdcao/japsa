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
 * 21/10/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/


package japsa.seq;


import japsa.seq.Alphabet;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import java.util.Arrays;




/**
 * A class to support reading sequences from raw, i.e.: every char is a symbol
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 *
 */
public class RawFormatReader extends SequenceReader{
	static int INI_SEQ_SIZE = 8192;//1 << 20;
	//Temporary byte array for the sequence
	private byte [] seq = new byte[INI_SEQ_SIZE];	
	private int seqIndex = 0;

	/**
	 * 
	 * Construct the reader from a file
	 */
	public RawFormatReader(InputStream ins) throws IOException{
		super(ins);
	}
	/**
	 * Construct the reader from a file
	 * @param fileName
	 * @throws IOException
	 */
	public RawFormatReader(String fileName) throws IOException{
		this(new FileInputStream(fileName));
	}


	public boolean hasNext() throws IOException{		
		return !eof;
	}

	public static Sequence read (InputStream ins, Alphabet alphabet) throws IOException{
		RawFormatReader reader = new RawFormatReader(ins);
		Sequence seq = reader.nextSequence(alphabet);
		reader.close();
		return seq;
	}

	/**
	 * Fast method to read a fasta file 
	 * Precondition:it is a fasta file. If some invalid character exists, an 
	 * IOException will be thrown.
	 * 
	 * @param in
	 * @param alphabet: the alphabet, if not specified, the default is DNA16
	 * @return
	 * @throws IOException
	 */
	public Sequence nextSequence(Alphabet alphabet) throws IOException{
		if (eof)
			return null;

		if (alphabet == null){
			alphabet = Alphabet.DNA16();//The most conservative
		}

		while (true){
			int nucleotide = alphabet.byte2index(currentByte);
			if (nucleotide >= 0){//valid nucleotide				
				//ensure the sequence array is big enough
				if (seqIndex >= seq.length) {// Full
					int newLength = seq.length * 2;
					if (newLength < 0) {
						newLength = Integer.MAX_VALUE;
					}
					// if the array is extended
					if (newLength <= seqIndex) {
							// in.close();
						throw new RuntimeException("Sequence is too long to handle");
					}
					// create new array of byte
					seq = Arrays.copyOf(seq, newLength);
				}
				//next symbol
				seq[seqIndex++] = (byte) nucleotide;
				//seqIndex++;					
			}else if(nucleotide == -1){
				throw new RuntimeException("Unexecpected character '" + (char) currentByte + "' for alphabet {" + alphabet + "} at the line  " + lineNo);
			}			
			if (!nextByte())
				break;
		}
		return new Sequence(alphabet, seq, seqIndex, "READ");
	}
}

