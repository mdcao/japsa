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
 * 23/01/2014 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsa.util;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

/**
 * @author minhduc
 *
 * We should use slf4j instead
 */
@Deprecated
public class Logging {
	private static String prefix = "#";

	private static PrintStream infoStr = System.err;
	private static PrintStream errorStr = System.err;
	private static PrintStream warnStr = System.err;


	/**
	 * A simple logging system
	 */
	private Logging() {
		// TODO Auto-generated constructor stub
	}	

	/**
	 * @param args
	 * @throws FileNotFoundException 
	 */
	public static void setInfoFile(String filePath) throws FileNotFoundException{
		infoStr = new PrintStream(new FileOutputStream(filePath,true));
	}
	
	public static void setWarnFile(String filePath) throws FileNotFoundException{
		warnStr = new PrintStream(new FileOutputStream(filePath,true));
	}
	
	public static void setErrorFile(String filePath) throws FileNotFoundException{
		errorStr = new PrintStream(new FileOutputStream(filePath,true));
	}
	/**
	 * Force all three streams to the same file
	 * @param filePath
	 * @throws FileNotFoundException
	 */
	public static void setFile(String filePath) throws FileNotFoundException{
		warnStr = errorStr = infoStr = new PrintStream(new FileOutputStream(filePath,true));
	}	
	
	
	public static void info(String msg) {
		synchronized(infoStr){
			infoStr.println(prefix+msg);
		}

	}

	public static void warn(String msg) {
		synchronized(warnStr){
			warnStr.println(prefix+msg);
		}
	}

	public static void error(String msg) {
		synchronized(errorStr){
			errorStr.println(prefix+msg);
		}
	}

	public static void exit(String msg, int status) {
		error(msg);		
		System.exit(status);
	}
}
