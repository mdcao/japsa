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

import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.IOException;


/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 * Program to trim reads in a fastq file and split the file to smaller pieces
 */
@Deployable(scriptName = "jsa.ngs.fqtrim",
           scriptDesc = "Trim reads from a fastq file and break the file to smaller ones")
public class FastQTrim {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		/*********************** Setting up script ****************************/
		Deployable annotation = FastQTrim.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options]", annotation.scriptDesc());		
		/**********************************************************************/

		cmdLine.addStdInputFile();			
		
		cmdLine.addString("output", null, "Name of output fastq file, output files will be added with suffix P<index>_", true);
		cmdLine.addInt("size", 0, "The number of reads per file, a negative number for not spliting");
		cmdLine.addInt("begin", 0, "Begin position of a read (1-index, inclusive) - 0 for not trimming");
		cmdLine.addBoolean("trim", false, "Whether to trim Ns at the 3' end. Note this will trim before begin/end trimming");
		cmdLine.addInt("end", 0, "End position of a read (1-index, inclusive)  - 0 for not trimming");	
		
		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/		

		String output = cmdLine.getStringVal("output");
		String inFile = cmdLine.getStringVal("input");
		
		int begin = cmdLine.getIntVal("begin");
		int end = cmdLine.getIntVal("end");			
		int size = cmdLine.getIntVal("size");
		boolean trim = cmdLine.getBooleanVal("trim");
		
		if (end > 0 && begin >= end) {
			Logging.exit("Begin "+(begin) + " must be smaller than end (" + end +")", -1);			
		}	

		if (size == 0)
			size = -1;
			
		int index = 1;
		SequenceOutputStream outStream = 
				 SequenceOutputStream.makeOutputStream("P"+index+"_" + output);

		BufferedReader reader = SequenceReader.openFile(inFile);
		String line = "";

		int count = 0;		
		int countAll = 0;
		
		while ( (line = reader.readLine()) != null){
			String name  = line.trim();
			if (name.charAt(0) != '@')
				throw new RuntimeException(name + " unexpected at read " + countAll);

			//standardise read name: replace ' ' space by a '_'
			//name = name.replace(' ', '_');
			String seq = reader.readLine();
			reader.readLine();//'+'
			String qual = reader.readLine();
			if (trim){
				int lastIndex = seq.length();
				while (lastIndex > 0 && seq.charAt(lastIndex - 1) =='N'){
					lastIndex --;
				}
				if (lastIndex < seq.length()){
					seq = seq.substring(0 , lastIndex);
					qual = qual.substring(0 , lastIndex);
				}
			}
			
			if (begin > 0) {
				if (begin > seq.length()){
					seq  = "";
					qual = "";
				}else if (seq.length() > end){
					seq = seq.substring(begin - 1 , end);
					qual = qual.substring(begin - 1 ,end);
				}else{
					seq = seq.substring(begin - 1);
					qual = qual.substring(begin - 1);
				}
			}else if (end > 0 && end < seq.length()){
				seq = seq.substring(0 , end);
				qual = qual.substring(0 ,end);
			}			
			
			if (count == size){
				//write this file
				outStream.close();

				///start a new one
				index ++;
				count = 0;
				outStream = //new SequenceOutputStream("P"+index+"_" + output);	
				SequenceOutputStream.makeOutputStream("P"+index+"_" + output);
			}
			
			count ++;
			countAll ++;
			
			outStream.print(name);
			outStream.print('\n');
			outStream.print(seq);
			outStream.print("\n+\n");		
			outStream.print(qual);
			outStream.print('\n');
		}//while
		outStream.close();
		Logging.info("Write " + countAll + " reads to " + index + " files" );		
	}
}
