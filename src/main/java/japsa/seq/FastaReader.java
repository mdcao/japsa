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

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import java.util.Arrays;
import java.util.zip.GZIPInputStream;




/**
 * A class to support reading sequences from fasta files
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 *
 */
public class FastaReader extends SequenceReader{
    static int INI_SEQ_SIZE = 8191;//1 << 20-1;
	//private int seqNo = 1;//keep tract of sequence number	

	//Temporary byte array for the sequence
	private byte [] seq = new byte[INI_SEQ_SIZE];
	private int seqIndex = 0;

	/**
	 * 
	 * Construct the reader from a file
	 */
	public FastaReader(InputStream ins) throws IOException{
		super(ins);
	}
	/**
	 * Construct the reader from a file
	 * @param fileName
	 * @throws IOException
	 */
	public FastaReader(String fileName) throws IOException{
		this(new FileInputStream(fileName));
	}

	/* (non-Javadoc)
	 * @see japsa.seq.SequenceFileReader#hasNext()
	 */
	//@Override

	public boolean hasNext() throws IOException{		
		while (true){
			if (currentByte == 62 && eol)
				return true;
			if (!nextByte()) return false;
		}
	}
	/**
	 * Read the first sequence from an input stream
	 * @param ins
	 * @param alphabet
	 * @return
	 * @throws IOException
	 */

	public static Sequence read (InputStream ins, Alphabet alphabet) throws IOException{
		FastaReader reader = new FastaReader(ins);
		Sequence seq = (reader).nextSequence(alphabet);
		reader.close();
		return seq;
	}

	/**
	 * Fast method to read a fasta file 
	 * Precondition:it is a fasta file. If some invalid character exists, an 
	 * IOException will be thrown.
	 * 
	 * @param alphabet: the alphabet
	 * @return
	 * @throws IOException
	 */
	public Sequence nextSequence(Alphabet alphabet) throws IOException{
		if (eof)
			return null;


		if (alphabet == null){
			alphabet = Alphabet.DNA16();//The most conservative
		}

		StringBuilder header = new StringBuilder();
		seqIndex = 0;//the number of nucletiodes read

		if (currentByte != 62){// 62 = '>'
			throw new RuntimeException("> is expected at the start line " + lineNo + ", found " + ((char)currentByte));
		}

		//Read the header
		while (!eol){
			nextByte();
			header.append((char) currentByte);
		}

		while (true){
			if (nextByte()){
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
					}
					//next symbol
					seq[seqIndex++] = nucleotide;
					//seqIndex++;					
				}else  if(currentByte == 62){
					//start of a new sequence
					//seqNo ++;
					return makeSequence(alphabet, seq, seqIndex, header.toString().trim());
				}else{
					if (nucleotide == -1){
						throw new RuntimeException("Unexecpected character '" + (char) currentByte + "' for dna {" + alphabet + "} at the line  " + lineNo);						
					}
				}
			}//if
			else{
				//seqNo ++;
				return makeSequence(alphabet, seq, seqIndex, header.toString().trim());
			}			

		}
	}
	

	static private Sequence makeSequence
	  (Alphabet alphabet, byte[] byteArray, int length, String name){
		String[] toks = name.split("\\s",2);
		if (toks.length >= 2){
			Sequence seq = new Sequence(alphabet, byteArray, length, toks[0]);
			seq.setDesc(toks[1]);
			return seq;
		}else{
			return new Sequence(alphabet, byteArray, length, name);
		}
	}
	/**
	 * A class to support fast reading sequences from fasta files by passing some
	 * of the checking. The class assumes a fasta sequence starts from the 
	 * character '>' regardless of where it is in the line.
	 * 
	 * This class should be used only for very large files, and with special care. 
	 * For most cases, FastaReader class should be used.
	 *	   
	 * @author Minh Duc Cao (http://www.caominhduc.org/)
	 * FIXME: For reading from stdin, this reader does not terminate upon end of
	 * stream signal.
	 */

	public static class Faster{	

		/**
		 * The internal buffer array where the data is stored. When necessary,
		 * it may be replaced by another array of
		 * a different size.
		 */
		private byte [] buff = new byte[BUFF_SIZE];
		private int count = -1; //the index of valid buffer
		private int pos = 0;//place to read next
		private int seqNo = 1;//keep tract of sequence number	

		//Temporary byte array for the sequence
		private byte [] seq = new byte[INI_SEQ_SIZE];
		private int seqIndex = 0;

		/**
		 * The underlying input stream 
		 */
		private InputStream in;	
		/**
		 * 
		 * Construct the reader from a file
		 */
		public Faster(InputStream ins) throws IOException{
			/**************************************************************/
			if (!ins.markSupported()) {
				ins = new BufferedInputStream(ins);
			}
			ins.mark(2);
			int magic = ins.read() & 0xff | ((ins.read() << 8) & 0xff00);
			ins.reset();


			if (magic == GZIPInputStream.GZIP_MAGIC)
				this.in = new GZIPInputStream(ins);		
			else
				/**************************************************************/
				this.in = ins;
		}
		/**
		 * Construct the reader from a file
		 * @param fileName
		 * @throws IOException
		 */
		public Faster(String fileName) throws IOException{
			this(new FileInputStream(fileName));
		}

		public void close() {		
			//in.close();//DO I need this?
		}


		/**
		 * This method moves the pointer to the start of the next sequence in the 
		 * stream if it does not currently point to the beginning of a sequence.
		 */
		public boolean hasNext() throws IOException{
			//if has processed all bytes read
			if (pos >= count){
				count = in.read(buff);
				if (count <= 0) return false;
				pos = 0;
			}
			//assert pos < count

			//if the reader points to the start of the next sequence
			if (buff[pos] == 62){// 62 = '>'
				return true;
			}

			//if not, looking for the next '>'
			while (true){
				pos ++;
				if (pos >= count){
					count = in.read(buff);
					if (count <= 0) return false;
					pos = 0;
				}
				if (buff[pos] == 62){// 62 = '>'
					return true;
				}				
			}


		}

		public static Sequence read (InputStream ins, Alphabet alphabet) throws IOException{
			return (new Faster(ins)).nextSequence(alphabet);
		}
		/**
		 * Fast method to read a fasta file 
		 * Precondition:it is a fasta file. If some invalid character exists, an 
		 * IOException will be thrown.
		 * 
		 * @return
		 * @throws IOException
		 */
		public Sequence nextSequence(Alphabet alphabet) throws IOException{			 				
			//array to hold indexes of nucleotides
			if (alphabet == null){
				alphabet = Alphabet.DNA16();//The most conservative
			}

			StringBuilder header = new StringBuilder();

			seqIndex = 0;//the number of nucletiodes read		
			if (pos >= count){
				count = in.read(buff);
				pos = 0;
			}
			//No more in the stream
			if (count <=0 ) return null;

			if (buff[pos++] != 62){// 62 = '>'
				pos --;
				throw new RuntimeException("> is expected at the start of Fasta No " + seqNo);
			}

			boolean seqMode = false; 

			for (;;){
				//make sure there is something in the buffer
				//this could be replaced with nextByte()
				if (pos >= count){
					count = in.read(buff);
					if (count <= 0) break;//no more to read
					pos = 0;
				}
				//process buffer
				while (pos < count){
					byte currentByte = buff[pos++];				

					if (seqMode){//reading header
						//reading sequence					
						byte nucleotide = alphabet.byte2index(currentByte);
						//assert nucleotide < dna.size()
						if (nucleotide >= 0){
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
							}
							seq[seqIndex++] = nucleotide;
							//seqIndex++;					
						}else  if(currentByte == 62){
							//start of a new sequence
							pos --;
							seqNo ++;
							return new Sequence(alphabet, seq, seqIndex, header.toString());
						}else{
							if (nucleotide == -1){
								throw new RuntimeException("Unexecpected character '" + (char) currentByte + "' for dna {" + alphabet + "} for sequence  " + seqNo);							
							}
						}
					}//mode == 0
					else{
						if (currentByte == 10 ){//CR 
							seqMode = true;
							//						lineNo ++;
							continue;
						}else if (currentByte == 13){//LF
							seqMode = true;
							continue;
						}
						header.append((char) currentByte);
					}
				}//while					
			}//for
			seqNo ++;
			return new Sequence(alphabet, seq, seqIndex, header.toString());
		}	
	}
}

