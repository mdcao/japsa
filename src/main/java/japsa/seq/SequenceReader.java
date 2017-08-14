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
 * 22/03/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.seq;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;


/**
 * This class provides a scheleton for reading sequence from various file
 * formats.
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
public abstract class SequenceReader implements Closeable{

	/**
	 * Provide a skeleton of a sequence file reader. Assume that the file
	 * may contains more than one sequences with the same format.
	 * Known format: 
	 * Fasta: start with '>"
	 * Fastq: start with @
	 * @param args
	 */		
	
	//implement new line character: CR (Mac), LF (Unix) of CR follows by LF (Wins)
	protected static final byte CR = 10, LF = 13;	
	static  final int BUFF_SIZE = 8192;

	/**
	 * The internal buffer array where the data is stored. When necessary,
	 * it may be replaced by another array of
	 * a different size.
	 */
	private byte [] buff = new byte[BUFF_SIZE];	
	private int count = -1; //the index of valid buffer
	private int pos = 0;//place to read next
	
	protected int lineNo = 0;	
	protected byte currentByte = -1;
	protected boolean eol = true;	
	protected boolean eof = false;	

	/**
	 * The underlying input stream 
	 */
	private InputStream in;
	
	/**
	 * 
	 * Construct the reader from a file
	 */
	public SequenceReader(InputStream ins) throws IOException{
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
		
		//Read the first byte
		nextByte();
	}
	/**
	 * Construct the reader from a file
	 * @param fileName
	 * @throws IOException
	 */
	public SequenceReader(String fileName) throws IOException{
		this(new FileInputStream(fileName));
	}
	
	/**
	 * Close the stream
	 * @throws IOException
	 */

	public void close() throws IOException{		
		in.close();
	}	
	

	protected byte[] nextLine = new byte[1024];
	protected int nextLineLength = 0;
	/**
	 * Read the next line (from current pointer to the next eol)
	 * @return
	 */
	
	protected int nextLine() throws IOException{	
		if (eof) return 0;
		nextLineLength = 0;
		if (currentByte != LF && currentByte != CR)
			nextLine[nextLineLength ++] = currentByte;
		
		while (nextByte()){	
			if (nextLineLength >= nextLine.length){
				//double the buffer
				int newLength = nextLine.length * 2;
				if (newLength < 0)
					newLength = Integer.MAX_VALUE;
				nextLine = Arrays.copyOf(nextLine, newLength );
				
			}
			nextLine [nextLineLength ++] = currentByte;
			
			if (eol) 
				return nextLineLength;
		}
		return nextLineLength;
	}
	
	public String getCurrentLine(){
		return new String(nextLine, 0, nextLineLength);
	}
	
	/**
	 * Read the next byte in the stream into currentByte
	 * @return false if end of stream
	 * @throws IOException
	 */
	protected final boolean nextByte() throws IOException{
		
		//advance the pointer 
		pos ++;
		
		//ensure more data from the stream
		if (pos >= count){
			count = in.read(buff);
			if (count <= 0) {
				eof = true;
				eol = true;//eof implies eol
				return false;
			}
			pos = 0;
		}		
		//assert pos < count
		
		//Check if an end of line is reached
		if (buff[pos] == CR){
			eol = true;lineNo ++;			
		}else if (buff[pos] == LF){
			if (currentByte != CR){
				eol = true;
				lineNo ++;
			}else
				//effectively ignore the LF that follows CR
				return nextByte();
		}else{
			eol = false;
		}		
		currentByte = buff[pos];
		
		return true;
	}
	
	/**
	 * 
	 * Read the next n bytes to store in the buff. Return the number of bytes
	 * read.
	 * @param buf
	 * @throws IOException
	 */
	public int read(byte[] buf) throws IOException{
		int n = 0;
		
		for (; n < buf.length && !eof; n++){
			buf[n] = currentByte;
			nextByte();
		}			
		
		return n;
		
	}
		
	/**
	 * Return the next sequence from the stream
	 * @return
	 */
	public abstract Sequence nextSequence(Alphabet alphabet)  throws IOException;
		
	
	
	/**
	 * Open a file for other read operations. If the file is zipped, this method
	 * open an unzipped reader stream for normal reading. If the file name 
	 * supplied as a "-", it will return the stream from standard input.	
	 * @param fileName
	 *            Name of the file to open
	 * @return the buffered reader associate with the file
	 * @throws IOException
	 */
	public static BufferedReader openFile(String fileName) throws IOException {
		
		InputStream is;
		if ("-".equals(fileName))
			is = new BufferedInputStream(System.in);
		else
			is = new BufferedInputStream(new FileInputStream(fileName));
		
		return openFile(is);
	}
	

	/**
	 * Open an input stream for other read operations. If stream is zipped, this method
	 * open an unzipped reader stream for normal reading.
	 * 
	 * TODO: it seems like this method returns too many level or wrap
	 * @param is the input stream
	 * @return the buffered reader associate with the file
	 * @throws IOException
	 */
	
	public static BufferedReader openFile(InputStream is) throws IOException {
		int magic;		
		
		is.mark(4);
		magic = is.read();
		magic += is.read() * 256;
		is.reset();
		
		if (magic == GZIPInputStream.GZIP_MAGIC)
			is = new GZIPInputStream(is);
		
		return new BufferedReader(new InputStreamReader(is));
	}
	
	/**
	 * Predict the format of the file and invoke the appropriate Reader. Return
	 * null if the format is not recognised.
	 * This class currently only recognises JAPSA, Fasta and Fastq format. 
	 * 
	 * @param filename
	 * @return the file reader, or null if the format is not recognised
	 * @throws IOException
	 */
	
	public static SequenceReader getReader (String filename) throws IOException{
		if ("-".equals(filename))
			return getReader(System.in);
		return getReader(new FileInputStream(filename));
	}
	
	/**
	 * Predict the format of the file and invoke the appropriate Reader. Return
	 * null if the format is not recognised.
	 * This class currently only recognises JAPSA, Fasta and Fastq format. 
	 * 
	 * @param filename
	 * @return the file reader, or null if the format is not recognised
	 * @throws IOException
	 */	
	public static SequenceReader getReader (InputStream ins) throws IOException{
		//Make sure the stream supports mark
		if (!ins.markSupported()) {
			ins = new BufferedInputStream(ins);
		}
		
		//check if a zip
		ins.mark(10);
		byte [] myBuf = new byte[10];
		ins.read(myBuf);
		ins.reset();		
		int magic = myBuf[0] & 0xff | ((myBuf[1] << 8) & 0xff00);
		
		if (magic == GZIPInputStream.GZIP_MAGIC){
			ins = new BufferedInputStream(new GZIPInputStream(ins));
			ins.mark(10);
			ins.read(myBuf);
			ins.reset();
		}
		
		String format = new String(myBuf);		
		
		if (format.startsWith(">"))
			return new FastaReader(ins);
		else if (format.startsWith("@"))
			return new FastqReader(ins);
		else if (format.startsWith(JapsaFileFormat.HEADER))
			return new JapsaFileFormat(ins);
			
		return null;
	}	
	
	/**
	 * Read all sequences and return an array list
	 * 
	 * @param fileName
	 * @param dna
	 * @return
	 * @throws IOException
	 */
	public static ArrayList<Sequence> readAll(String fileName, Alphabet alphabet) throws IOException{
		ArrayList<Sequence> seqs = new ArrayList<Sequence>();
		SequenceReader reader = getReader(fileName);
		File f = new File(fileName);
		if(!f.exists()) throw new IOException ("File does not exist "+fileName);
		if(reader ==null) throw new IOException("Reader is null");
		if(alphabet ==null) throw new IOException("alphabet is null");
		while (true){
			Sequence seq = reader.nextSequence(alphabet);
			if (seq == null)
				break;
			else
				seqs.add(seq);
		}
		return seqs;
	}
	
}
