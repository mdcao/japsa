/*****************************************************************************
 * Copyright (c) 2017 Minh Duc Cao (minhduc.cao@gmail.com).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the names of the institutions nor the names of the contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/
/*                           Revision History                                
 * 12-03-2017 - Minh Duc Cao: Created       
 *                                
 ****************************************************************************/
package japsa.tools.bio.amra;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 * 
 */
@Deployable(
		scriptName = "jsa.amra.assppro",
		scriptDesc = "Extract subsequences"
		)
public class AssemblyPostProcessingCmd extends CommandLine {
	public AssemblyPostProcessingCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 
		///////////////////////////////////////////////////////////////

		addString("sample", null, "Sample ID", true);
		addString("input", null, "Name of the input file, - for standard input", true);
		addString("output", null, "Name of the output file, - for standard output", true);		
		addString("summary", null, "Name of the summary file", true);

		addStdHelp();
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new AssemblyPostProcessingCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String sample      = cmdLine.getStringVal("sample");
		String inputFile   = cmdLine.getStringVal("input");
		String outputFile  = cmdLine.getStringVal("output");
		String summaryFile = cmdLine.getStringVal("summary");
		
		String summaryStr = postProcessing(sample, inputFile, outputFile);
		if (summaryStr != null){
			Files.write(Paths.get(summaryFile), (summaryStr + "\n").getBytes());
		}else
			System.exit(1);
	}

	public static String postProcessing( String sampleID, String input, String output) throws IOException{				
		ArrayList<Sequence> seqList = SequenceReader.readAll(input, Alphabet.DNA());
		if (seqList.size() ==0)
			return null;		
		Collections.sort(seqList, 
				new Comparator<Sequence>(){
			public int compare(Sequence seq1, Sequence seq2) {
				return Integer.compare(seq2.length(), seq2.length());				//	      
			}
		});
		
		int sum = 0;				
		//c. Rename contig and write to file
		int ind = 0;
		//how many chars needed 
		int fieldSize = 3;//String.valueOf(seqList.size()).length();
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);
		for (Sequence seq:seqList){
			ind ++;
			String name = String.valueOf(ind);
			//pad in 0
			while (name.length() < fieldSize)
				name = "0" + name;

			seq.setDesc(seq.getName() +" " + seq.getDesc());
			seq.setName(sampleID + "_C" + name);
			seq.writeFasta(sos);
			sum += seq.length();
		}
		sos.close();
		
		//second round: compute N50
		int n50 = 0;
		int sumHalf = 0;
		for (Sequence seq:seqList){
			sumHalf += seq.length();
			if (sumHalf * 2 >= sum){
				n50 = seq.length();
				break;
			}
		}

		return "length " + sum + "\ncontig " + seqList.size() + "\nn50 " + n50;
	}

}