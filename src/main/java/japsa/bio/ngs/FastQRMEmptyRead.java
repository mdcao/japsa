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
 * 02/10/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.bio.ngs;

import japsa.seq.Alphabet;
import japsa.seq.FastqReader;
import japsa.seq.FastqSequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.IOException;


/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 * Program to trim reads in a fastq file and split the file to smaller pieces
 */
@Deployable(scriptName = "jsa.hts.fqrmempty",
           scriptDesc = "Remove empty reads from a fastq file, especially after trimming")
public class FastQRMEmptyRead {
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		/*********************** Setting up script ****************************/
		Deployable annotation = FastQRMEmptyRead.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options]", annotation.scriptDesc());		
		/**********************************************************************/

		
		
		cmdLine.addString("1", null, "Name of the first fastq file", true);
		cmdLine.addString("2", null, "Name of the second fastq file for paired end");	
		
		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/		

		String f1 = cmdLine.getStringVal("1");
		String f2 = cmdLine.getStringVal("2");		
		
		
		SequenceOutputStream out1 = 
				 SequenceOutputStream.makeOutputStream("R"+f1);		
		FastqReader reader1 = new FastqReader(f1);
		
		Alphabet.DNA dna = Alphabet.DNA(); 
		FastqSequence seq;
		
		if (f2 == null){						
			while ((seq = reader1.nextSequence(dna)) != null){
				if (seq.length() > 0)
					seq.print(out1);
			}
			reader1.close();
			out1.close();			
		}else{
			SequenceOutputStream out2 = 
					 SequenceOutputStream.makeOutputStream("R"+f2);		
			FastqReader reader2 = new FastqReader(f2);
			
			while ((seq = reader1.nextSequence(dna)) != null){
				FastqSequence seq2 = reader2.nextSequence(dna);
				if (seq.length() > 0 && seq2.length() > 0){
					seq.print(out1);
					seq2.print(out2);
				}
			}
			reader1.close();
			out1.close();
			reader2.close();
			out2.close();	
			
		}
	}
}
