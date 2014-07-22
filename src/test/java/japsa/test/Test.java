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
 * 01/01/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.test;


import java.io.IOException;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;



/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
public class Test {
	  private final static Logger LOGGER = Logger.getLogger(Test.class.getName());
	  public final static boolean DEBUG_MODE = false;
	  
	  public final static int DEBUG_LEVEL = 1;
	  
	  public final static int PRODUCTION = 0;
	  public final static int WARNING = 1;
	  public final static int DEBUG = 2;
	  
	  
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		String s = "45";
		String [] ss = s.split("h");
		for (int i = 0; i < ss.length;i++)
			System.out.println(ss[i]);
		
		LOGGER.log(Level.FINE, "fine");
		
		LOGGER.log(Level.INFO, "info");
		
		
		Timer timer = new Timer();		
		//bio.la.Sequence japsa.seq = bio.la.Sequence.fasta(bio.la.Alphabet.dna, args[0], new java.util.Random());
		//System.out.println(japsa.seq.length());		
		
		//FastaFileReader reader = new FastaFileReader(args[0]);
		//japsa.seq.Sequence japsa.seq = reader.readFasta(null);
		//System.out.println(japsa.seq.length());
		Random random = new Random();
		double x = 1;
		for (int i = 1; i <= 100000000;i++){
			x += ((double)i)/208;			
		}
		
		timer.mark("done");
		
		LOGGER.log(Level.INFO, "done");
	}
}
