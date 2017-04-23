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


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;



/**
 * A class to support fast reading sequences from fasta files
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 *
 */
public class FastqReader extends SequenceReader{
    private static final Logger LOG = LoggerFactory.getLogger(FastqReader.class);
	
	private byte [] seq = new byte[1024];//estimated max read length
	private byte [] qual = new byte[1024];
	private byte [] name = new byte[512];		
	
	/**
	 * Construct the reader from a file
	 * @param fileName
	 * @throws IOException
	 */
	public FastqReader(String fileName) throws IOException{
		super(fileName);
	}
	
	public FastqReader(InputStream ins) throws IOException{
		super(ins);
	}

	
	/* (non-Javadoc)
	 * @see japsa.seq.SequenceFileReader#nextSequence()
	 */
	@Override
	public FastqSequence nextSequence(Alphabet alphabet) throws IOException {
		//no more to read
		if (eof)
			return null;
		
		if (currentByte != 64)			
			throw new RuntimeException("@ is expected at the start line  " + lineNo);
			
		int nameIndex = 0, seqIndex = 0, qualIndex = 0;
		
		while (nextByte() && !eol){//read header
			//ensure nameBuf is big enough
			if (nameIndex >= name.length){
				int newLength = name.length * 2;
				if (newLength < 0) {
					newLength = Integer.MAX_VALUE;
				}
				// if the array is extended
				if (newLength <= seqIndex) {
					throw new RuntimeException(
							"Name is too long to handle at line " + lineNo);
				}//if
			}//if
			//assert nameIndex < name.length			
			name[nameIndex ++ ] = currentByte;
		}
		
		//read sequence		
		while (nextByte() && !eol){//read header			
			byte nucleotide = alphabet.byte2index(currentByte);
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
						throw new RuntimeException(
								"Sequence is too long to handle");
					}
					// create new array of byte
					seq = Arrays.copyOf(seq, newLength);
					qual = new byte[newLength];
				}
				//next symbol
				seq[seqIndex++] = nucleotide;
				//seqIndex++;					
			}else{
				throw new RuntimeException("Unexecpected character '" + (char) currentByte + "' for dna {" + alphabet + "} at the line  " + lineNo); 
			}			
		}//while
		
		//skip the next line
		while (nextByte() && !eol){//read header			
		}	
		
		//quality line
		while (nextByte() && !eol){//read header
			qual[qualIndex ++] = currentByte;
		}
		
		if (seqIndex != qualIndex){
			//throw new RuntimeException("Lengths of sequence and quality strings do not match at line " + lineNo + " : " + seqIndex + " vs " + qualIndex);
			LOG.warn("Lengths of sequence and quality strings do not match at line " + lineNo + " : " + seqIndex + " vs " + qualIndex);
		}
		
		//Read the next byte from the stream (expecting a @ or eof
		nextByte();
		
		return new FastqSequence(alphabet, seq, qual, seqIndex, new String(name,0, nameIndex));
	}
	
}
