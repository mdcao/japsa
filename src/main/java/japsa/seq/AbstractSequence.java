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
 * 01/01/2012 - Minh Duc Cao: Created
 * 30/12/2012 - Minh Duc Cao: revised                                        
 *  
 ****************************************************************************/

package japsa.seq;



import japsa.seq.Alphabet;
import japsa.util.JapsaMath;

import java.io.IOException;





/**
 * Abstract class for a sequence
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 * 
 */

public abstract class AbstractSequence implements
Comparable<japsa.seq.AbstractSequence>, CharSequence {

	/**
	 * The dna of the sequence, this should not be changed during the time
	 * life of the sequence
	 */
	private final Alphabet alphabet;

	// Some default ID and description of the sequence
	private String name = "", desc = "";

	public AbstractSequence(Alphabet alphabet) {
		this.alphabet = alphabet;
	}

	public AbstractSequence(Alphabet alphabet, String name) {
		this.alphabet = alphabet;
		this.name = name;
	}

	public AbstractSequence(Alphabet alphabet, String name, String desc) {
		this.alphabet = alphabet;
		this.name = name;
		this.desc = desc;
	}


	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @param name 
	 *    the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * @return the desc
	 */
	public String getDesc() {
		return desc;
	}

	/**
	 * @param desc
	 *            the desc to set
	 */
	public void setDesc(String desc) {
		this.desc = desc;
	}


	/**
	 * Return the dna of the sequence
	 * 
	 * @return
	 */
	public Alphabet alphabet() {
		return alphabet;
	}
	/**
	 * 
	 * @return the length of the sequence
	 */
	public abstract int length();

	/**
	 * 
	 * @param loc the location
	 * @return the index of the symbol at location loc
	 */
	public abstract int symbolAt(int loc);
	public abstract void setSymbol(int loc, int symbol);
	
	public abstract byte getBase(int loc);
	public abstract byte setBase(int loc, byte base);
	

	
	/**
	 * Convert the sequence to a byte array
	 * @return
	 */
	//public abstract byte[] toBytes();

	/**
	 * Convert the sequence to an array of chars
	 * @return
	 */
	public char[] charSequence() {
		char[] charSeq = new char[length()];
		for (int i = 0; i < charSeq.length; i++) {
			charSeq[i] = charAt(i);
		}
		return charSeq;
	}

	/**
	 * 
	 * @param loc
	 * @return character at the location (0-base)
	 */
	public char charAt(int loc) {
		return alphabet.int2char(symbolAt(loc));
	}


	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(japsa.seq.AbstractSequence o) {
		// TODO Auto-generated method stub
		return 0;
	}

	/**
	 * Check if this sequence matches another sequence. Return 0 if they match, 
	 * return a negative value if the alphabets or lengths differ, a positive 
	 * number of the first mismatched base. The definition of match is 
	 * characterised by the dna. For example, if DNA5, A matches A and
	 * does not match C or N, while N matches every other symbol and itself.
	 * 
	 * Note that if japsa.seq A matches japsa.seq B, the reverse does not neccessarily hold.
	 * @param another
	 * @return 0 if the two sequences identical, -2 if different alphabet, 
	 * -1 if different length, and a positive number as the first position
	 * the two sequences differ 
	 */
	public int match(Sequence another){
		if (this.alphabet() != another.alphabet())
			return -2;
		if (length() != another.length())
			return -1;

		for (int i = 0; i < length() ; i++)
			if (!alphabet().match(symbolAt(i), another.symbolAt(i)))
				return i + 1;

		return 0;
	}
	
	/**
	 * Write the sequence to a stream on a line, no formating
	 * @param out
	 * @throws IOException
	 */
	public void writeStream(SequenceOutputStream out) throws IOException {
		for (int i = 0; i < length(); i++) {
			out.print(charAt(i));			
		}
	}


	/**
	 * Write the sequence to an output stream in fasta format. Make sure that 
	 * this is the start of the file, or a newline was just written. For big 
	 * sequence, a buffered stream should be used.
	 * 
	 * @param out: The stream to write the sequence to
	 * @throws IOException
	 */
	public void writeFasta(SequenceOutputStream out) throws IOException {
		out.print('>');
		out.print(name);
		if (desc.length() > 0){
			out.print(' ');
			out.print(desc.replaceAll("\n+", ";"));//make sure no extra new line is written 
		}

		for (int i = 0; i < length(); i++) {
			if (i % 100 == 0)
				out.print('\n');
			out.print(charAt(i));			
		}
		out.print('\n');
	}


	/**
	 * Write the sequence to a file with filename
	 * @param fileName: name of the file to write to
	 * @throws IOException
	 */
	public void writeFasta(String fileName) throws IOException {
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(fileName);
		writeFasta(out);
		out.close();
	}

	/**
	 * Write any description of the sequence to the stream
	 * @param out
	 * @throws IOException
	 */
	public void writeDescription(SequenceOutputStream out) throws IOException{
		if (desc.length() > 0){
			String [] toks = desc.trim().split("\n");
			for (int i = 0; i < toks.length; i++){
				out.print(JapsaFileFormat.SEQUENCE_COMMENT);
				out.print(toks[i]);
				out.print('\n');
			}		
		}
	}
	
	/**
	 * Write in JSA format to an output stream. For speed, a buffered stream should
	 * be used. The pointer should be at the beginning of a line
	 * TODO: test this
	 * @param out
	 * @throws IOException
	 */
	public void writeJSA(SequenceOutputStream out) throws IOException{
		out.print(JapsaFileFormat.HEADER);
		out.print(JapsaFileFormat.DELIMITER);
		out.print(name);
		out.print(JapsaFileFormat.DELIMITER);
		out.print(length());
		out.print(JapsaFileFormat.DELIMITER);
		out.print(alphabet.getName());
		out.print('\n');
		
		writeDescription(out);
		
		int radixSize = JapsaMath.radixSize(length());

		for (int x = 0; x < length(); x++) {
			char c = charAt(x);
			if (x % JapsaFileFormat.CHAR_PER_LINE == 0) {
				out.print('\n');
				out.print(x+1, radixSize);
				out.print(' ');
				out.print(' ');				
			} else if (x % JapsaFileFormat.CHAR_PER_BLOCK == 0) {
				out.print(' ');
			}
			out.print(c);
		}
		out.print('\n');
		out.print('\n');
	}


	/**
	 * Write the sequence to a file in JSA format
	 * @param fileName
	 * @throws IOException
	 */
	public void writeJSA(String fileName) throws IOException {
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(fileName);
		writeJSA(out);
		out.close();
	}


	/**
	 * Write to an output stream. For speed, a buffered stream should
	 * be used. The pointer should be at the beginning of a line
	 * @param out
	 * @throws IOException
	 */
	public void print(SequenceOutputStream out) throws IOException{
		//writeJSA(out);
		writeFasta(out);
	}

	
	/**
	 * Write the sequence to a file in JSA format
	 * @param fileName
	 * @throws IOException
	 */
	public void write(String fileName) throws IOException {
		writeJSA(fileName);		
	}	

	/**
	 * Return the string represent the sequence
	 */		

	public String toString(){
		StringBuilder sb = new StringBuilder(length());
		for (int i = 0; i < length(); i++)
			sb.append(this.charAt(i));
		return sb.toString();
	}

}

