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
 * 30/03/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.seq;

import japsa.util.JapsaMath;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.GZIPOutputStream;


/**
 * This is a re-implementation of the standard BufferedOutputStream class.
 * Note this class is not synchronised.
 * Methods write() are used for native types (such as a byte or an int which 
 * is also convereted to byte. For writing the actual value (string, numbers, 
 * char etc, use methods print().
 *  
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
public class SequenceOutputStream  extends OutputStream {
	/**
	 * The underlying output stream to be filtered.
	 */
	protected OutputStream out;

	/**
	 * The internal buffer where data is stored.
	 */	
	protected byte buf[];

	/**
	 * The number of valid bytes in the buffer. This value is always
	 * in the range <tt>0</tt> through <tt>buf.length</tt>; elements
	 * <tt>buf[0]</tt> through <tt>buf[count-1]</tt> contain valid
	 * byte data.
	 */
	protected int count;

	/**
	 * Creates a new buffered output stream to write data to the
	 * specified underlying output stream.
	 *
	 * @param   out   the underlying output stream.
	 */
	public SequenceOutputStream(OutputStream out) {
		this(out, 8192);
	}

	/**
	 * Creates a new buffered output stream to write data to the
	 * specified underlying output stream with the specified buffer
	 * size.
	 *
	 * @param   out    the underlying output stream.
	 * @param   size   the buffer size.
	 * @exception IllegalArgumentException if size &lt;= 0.
	 */
	public SequenceOutputStream(OutputStream out, int size) {
		this.out = out;	        
		if (size <= 0) {
			throw new IllegalArgumentException("Buffer size <= 0");
		}
		buf = new byte[size];
	}
	
	/**
	 * Create an output stream to a file. If the file name ends with .gz, a gzip
	 * stream is created.
	 * If the string "-" is passed, the stream is bound to standard output
	 * @param fileName
	 * @return
	 * @throws IOException
	 */
	public static SequenceOutputStream makeOutputStream(String fileName)throws IOException{
		if (fileName.endsWith(".gz"))
			return new SequenceOutputStream 
					(new GZIPOutputStream 
					  (new FileOutputStream(fileName)));		
		else if (fileName.equals("-"))
			return 	new SequenceOutputStream(System.out);
		
		return new SequenceOutputStream (new FileOutputStream(fileName));
		
	}

	/**
	 * Writes <code>b.length</code> bytes to this output stream.
	 * <p>
	 * The <code>write</code> method of <code>FilterOutputStream</code>
	 * calls its <code>write</code> method of three arguments with the
	 * arguments <code>b</code>, <code>0</code>, and
	 * <code>b.length</code>.
	 * <p>
	 * Note that this method does not call the one-argument
	 * <code>write</code> method of its underlying stream with the single
	 * argument <code>b</code>.
	 *
	 * @param      b   the data to be written.
	 * @exception  IOException  if an I/O error occurs.
	 */
	public void write(byte b[]) throws IOException {
		write(b, 0, b.length);
	}

	/**
	 * Writes the specified byte to output stream.
	 *
	 * @param      b   the byte to be written.
	 * @exception  IOException  if an I/O error occurs.
	 */
	public void write(int b) throws IOException {
		if (count >= buf.length) {
			flushBuffer();
		}
		buf[count++] = (byte)b;
	}
	
	/**
	 * Write the specified byte to this stream
	 * @param b
	 * @throws IOException
	 */
	public void write(byte b) throws IOException {
		if (count >= buf.length) {
			flushBuffer();
		}
		buf[count++] = b;
	}


	/**
	 * Writes <code>len</code> bytes from the specified byte array
	 * starting at offset <code>off</code> to this buffered output stream.
	 *
	 * <p> Ordinarily this method stores bytes from the given array into this
	 * stream's buffer, flushing the buffer to the underlying output stream as
	 * needed.  If the requested length is at least as large as this stream's
	 * buffer, however, then this method will flush the buffer and write the
	 * bytes directly to the underlying output stream.  Thus redundant
	 * <code>BufferedOutputStream</code>s will not copy data unnecessarily.
	 *
	 * @param      b     the data.
	 * @param      off   the start offset in the data.
	 * @param      len   the number of bytes to write.
	 * @exception  IOException  if an I/O error occurs.
	 */
	public void write(byte b[], int off, int len) throws IOException {
		if (len >= buf.length) {
			/* If the request length exceeds the size of the output buffer,
	               flush the output buffer and then write the data directly.
	               In this way buffered streams will cascade harmlessly. */
			flushBuffer();
			out.write(b, off, len);
			return;
		}
		if (len > buf.length - count) {
			flushBuffer();
		}
		System.arraycopy(b, off, buf, count, len);
		count += len;
	}

	/**
	 * Flushes this buffered output stream. This forces any buffered
	 * output bytes to be written out to the underlying output stream.
	 *
	 * @exception  IOException  if an I/O error occurs.
	 */
	public void flush() throws IOException {
		flushBuffer();
		out.flush();
	}

	/**
	 * Closes this output stream and releases any system resources
	 * associated with the stream.
	 * <p>
	 * The <code>close</code> method of <code>FilterOutputStream</code>
	 * calls its <code>flush</code> method, and then calls the
	 * <code>close</code> method of its underlying output stream.
	 *
	 * @exception  IOException  if an I/O error occurs.
	 */
	public void close() throws IOException {
		try {
			flush();
		} catch (IOException ignored) {
		}
		out.close();
	}



	/** Flush the internal buffer */
	private void flushBuffer() throws IOException {
		if (count > 0) {
			out.write(buf, 0, count);
			count = 0;
		}
	}


	
	/**
	 * Write the character to the stream
	 * @param c
	 * @throws IOException
	 */
	
	public void print(char c) throws IOException {
		if (count >= buf.length) {
			flushBuffer();
		}
		buf[count++] = (byte) c;
	}
	
	public void println() throws IOException {
		print('\n');
	}
	
	public void printFormat(String str, int field) throws IOException{
		if (field > str.length()){
			print(str);
			for (int i = str.length();i< field;i++)
				print(' ');
		}else{
			print(str.substring(0,  field));
		}		
	}
	
	/**
	 * Write every byte of the string to the stream
	 * @param str
	 * @throws IOException
	 */
	public void print(String str)throws IOException {
		write(str.getBytes());
	}
	
	
	/**
	 * Write the string presenting the number to the stream
	 * @param num: the number to write
	 * @throws IOException
	 */
	public void print(int num)throws IOException {
		//FIXME: this is slow
		print(Integer.toString(num));		
	}
	
	public void print(long num)throws IOException {
		//FIXME: this is slow
		print(Long.toString(num));		
	}
	
	public void print(float num)throws IOException {
		//FIXME: this is slow
		print(Float.toString(num));		
	}
	
	public void print(double num)throws IOException {
		//FIXME: this is slow
		print(Double.toString(num));		
	}
	
	
	public void print(long num, int fm)throws IOException {
		write(JapsaMath.formatLong((int)num,fm));		
	}
	
	/**
	 * Write the string presenting the number in fm spaces, right justified, to
	 * the stream. 
	 * @param num: the number to write
	 * @param fm: the number of spaces
	 * @throws IOException
	 */
	public void print(int num, int fm)throws IOException {
		write(JapsaMath.formatInt(num,fm));		
	}


}
