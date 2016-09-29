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

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.zip.GZIPOutputStream;

import japsa.util.Logging;


/**
 * Represent a collection of output streams. An operation on this stream will 
 * be done on all the underlaying stream. When an underlaying stream fails to 
 * function, it will be removed from the collection.  
 *  
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
public class CompositeOutputStream  extends OutputStream {

	ArrayList< OutputStream> outs;

	/**
	 * Creates a new buffered output stream to write data to the
	 * specified underlying output stream.
	 *
	 * @param   out   the underlying output stream.
	 */
	public CompositeOutputStream() {
		outs = new ArrayList< OutputStream>();
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
	public CompositeOutputStream(OutputStream out){
		this();
		outs.add(out);
	}

	public void addStream(OutputStream out){
		outs.add(out);
	}

	/**
	 * Create an output stream to a file. If the file name ends with .gz, a gzip
	 * stream is created.
	 * If the string "-" is passed, the stream is bound to standard output
	 * @param fileName
	 * @return
	 * @throws IOException
	 */
	public static CompositeOutputStream makeOutputStream(String fileName)throws IOException{
		if (fileName.endsWith(".gz"))
			return new CompositeOutputStream 
					(new GZIPOutputStream 
							(new FileOutputStream(fileName)));		
		else if (fileName.equals("-"))
			return 	new CompositeOutputStream(System.out);

		return new CompositeOutputStream (new FileOutputStream(fileName));
	}
	
	/**
	 * Writes <code>b.length</code> bytes to this output stream.
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
		Iterator<OutputStream> it = outs.iterator();
		while (it.hasNext()){
			try{
				it.next().write(b);
			}catch(IOException e){
				Logging.warn(e.getMessage());
				it.remove();
			}			
		}		
	}

	/**
	 * Write the specified byte to this stream
	 * @param b
	 * @throws IOException
	 */
	public void write(byte b) throws IOException {
		Iterator<OutputStream> it = outs.iterator();
		while (it.hasNext()){
			try{
				it.next().write(b);
			}catch(IOException e){
				Logging.warn(e.getMessage());
				it.remove();
			}			
		}
		
		
	}


	/**
	 * Writes <code>len</code> bytes from the specified byte array
	 * starting at offset <code>off</code> to this buffered output stream.
	 *
	 *
	 * @param      b     the data.
	 * @param      off   the start offset in the data.
	 * @param      len   the number of bytes to write.
	 * @exception  IOException  if an I/O error occurs.
	 */
	public void write(byte b[], int off, int len) throws IOException {
		Iterator<OutputStream> it = outs.iterator();
		while (it.hasNext()){
			try{
				it.next().write(b, off, len);
			}catch(IOException e){
				Logging.warn(e.getMessage());
				it.remove();
			}			
		}
	}

	/**
	 * Flushes this buffered output stream. This simply flushes all the 
	 * underlying streams. If any stream throws an exception, the stream is removed 
	 *
	 * @exception  IOException  if an I/O error occurs.
	 */
	public void flush() throws IOException {
		Iterator<OutputStream> it = outs.iterator();
		while (it.hasNext()){
			try{
				it.next().flush();
			}catch(IOException e){
				Logging.warn(e.getMessage());
				it.remove();
			}			
		}
	}

	/**
	 * Closes this output stream by closing all the individual streams
	 *
	 * @exception  IOException  if an I/O error occurs.
	 */
	public void close() throws IOException {
		Iterator<OutputStream> it = outs.iterator();
		while (it.hasNext()){
			try{
				it.next().close();
			}catch(IOException e){
				Logging.warn(e.getMessage());
				it.remove();
			}			
		}	

	}
}
