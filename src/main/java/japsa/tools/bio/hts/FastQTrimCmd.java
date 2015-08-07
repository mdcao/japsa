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
package japsa.tools.bio.hts;

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
@Deployable(
	scriptName = "jsa.hts.fqtrim",
	scriptDesc = "Trim reads from a fastq file and break the file to smaller ones")
public class FastQTrimCmd extends CommandLine{	
	public FastQTrimCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addStdInputFile();			

		addString("output", null, "Name of output fastq file, output files will be added with suffix P<index>_", true);
		addInt("size", 0, "The number of reads per file, a negative number for not spliting");
		addInt("begin", 0, "Begin position of a read (1-index, inclusive) - 0 for not trimming");

		addBoolean("trim", false, "Whether to trim Ns at the 3' end. Note this will trim before begin/end trimming");
		addInt("qual", 0, "Minimum of qual to be trimed from the 3', 0 for not trimming");
		addInt("end", 0, "End position of a read (1-index, inclusive)  - 0 for not trimming");	

		addStdHelp();		
	} 

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{

		CommandLine cmdLine = new FastQTrimCmd();
		args = cmdLine.stdParseLine(args);			
		
		/**********************************************************************/
		String output = cmdLine.getStringVal("output");
		String inFile = cmdLine.getStringVal("input");

		int begin = cmdLine.getIntVal("begin");
		int end = cmdLine.getIntVal("end");			
		int size = cmdLine.getIntVal("size");
		int qual = cmdLine.getIntVal("qual");

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
			String qualStr = reader.readLine();

			int lastIndex = seq.length();
			if (trim){			
				while (lastIndex > 0 && seq.charAt(lastIndex - 1) =='N'){
					lastIndex --;
				}
			}

			while (lastIndex > 0 && qualStr.charAt(lastIndex - 1) - '!' < qual){
				lastIndex --;
			}			

			if (lastIndex < seq.length()){
				seq = seq.substring(0 , lastIndex);
				qualStr = qualStr.substring(0 , lastIndex);
			}

			if (begin > 0) {
				if (begin > seq.length()){
					seq  = "";
					qualStr = "";
				}else if (seq.length() > end){
					seq = seq.substring(begin - 1 , end);
					qualStr = qualStr.substring(begin - 1 ,end);
				}else{
					seq = seq.substring(begin - 1);
					qualStr = qualStr.substring(begin - 1);
				}
			}else if (end > 0 && end < seq.length()){
				seq = seq.substring(0 , end);
				qualStr = qualStr.substring(0 ,end);
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
			outStream.print(qualStr);
			outStream.print('\n');
		}//while
		outStream.close();
		Logging.info("Write " + countAll + " reads to " + index + " files" );		
	}
}
