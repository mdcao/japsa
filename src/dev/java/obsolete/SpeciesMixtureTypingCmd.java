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

/*****************************************************************************
 *                           Revision History                                
 * 7 Aug 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package obsolete;

import java.io.BufferedReader;
import java.io.IOException;

import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.IntArray;
import japsa.util.deploy.Deployable;


/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.np.speciesTyping", 
	scriptDesc = "Species typing using Nanopore Sequencing"
	)
@Deprecated
public class SpeciesMixtureTypingCmd extends CommandLine {

	public SpeciesMixtureTypingCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("output", "output.dat",  "Output file");		
		addString("bamFile", null,  "The bam file");		
		addString("indexFile", null,  "indexFile ");
		addString("hours", null,  "The file containging hours against yields, if set will output acording to tiime");
		addBoolean("GUI", false,  "Run on GUI");
		addInt("number", 50,  "Number of reads");
		addInt("timestamp", 0,  "Timestamp to check, if <=0 then use read number instead");
		addInt("sim", 0,  "Scale for simulation");
		addDouble("qual", 0,  "Minimum alignment quality");

		addStdHelp();		
	} 
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLine cmdLine = new SpeciesMixtureTypingCmd();		
		args = cmdLine.stdParseLine(args);		
		
		/**********************************************************************/

		String output    = cmdLine.getStringVal("output");
		String bamFile   = cmdLine.getStringVal("bamFile");			
		String indexFile = cmdLine.getStringVal("indexFile");
		String hours     = cmdLine.getStringVal("hours");
		boolean GUI      = cmdLine.getBooleanVal("GUI");
		int number       = cmdLine.getIntVal("number");
		double qual      = cmdLine.getDoubleVal("qual");
		
		SpeciesMixtureTyping paTyping = new SpeciesMixtureTyping(GUI);

		paTyping.simulation = cmdLine.getIntVal("sim");
		paTyping.qual = qual;

		if (hours !=null){
			BufferedReader bf = SequenceReader.openFile(hours);
			String line = bf.readLine();//first line -> ignore
			paTyping.hoursArray = new IntArray();
			paTyping.readCountArray = new IntArray();

			while ((line = bf.readLine())!= null){
				String [] tokens = line.split("\\s+");
				int hrs = Integer.parseInt(tokens[0]);
				int readCount = Integer.parseInt(tokens[2]);

				paTyping.hoursArray.add(hrs);
				paTyping.readCountArray.add(readCount);	
			}
			bf.close();
		}

		//	paTyping.prefix = prefix;
		paTyping.countsOS = SequenceOutputStream.makeOutputStream(output);
		paTyping.preTyping(indexFile);
		paTyping.typing(bamFile, number);
		paTyping.countsOS.close();
		paTyping.close();		

	}

}
