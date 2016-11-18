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
package japsadev.misc.obsolete;

import java.io.BufferedReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.IntArray;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deprecated
@Deployable(
	scriptName = "jsa.np.resistGenes", 
	scriptDesc = "Antibiotic resistance genes identification")
public class ResistanceGeneCmd extends CommandLine{	
	public ResistanceGeneCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("output", "output.dat",  "Output file");
		addString("bamFile", null,  "The bam file");
		//addString("geneFile", null,  "The gene file");
		//addString("figure", null,  "Figure file");
		addString("mcoordFile", null,  "Alignment with dnadiff between the genome (illumina) and genes database");

		addDouble("scoreThreshold", 1.5,  "The alignment score threshold");
		addString("msa", "kalign",
				"Name of the msa method, support poa, kalign, muscle and clustalo");
		addString("global", "needle",
				"Name of the global method, support needle and hmm");
		addString("tmp", "tmp/t",  "Temporary folder");
		addString("hours", null,  "The file containging hours against yields, if set will output acording to time");
		addInt("read", 500,  "Number of reads before a typing, NA if timestamp is set");		
		addInt("timestamp", 0,  "Number of seconds between internval");
		addDouble("il", 0.9,   "Threshold for Illumina");
		addBoolean("twodonly", false,  "Use only two dimentional reads");

		
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException, InterruptedException{
		CommandLine cmdLine = new ResistanceGeneCmd();
		args = cmdLine.stdParseLine(args);		

		String output = cmdLine.getStringVal("output");
		String bamFile = cmdLine.getStringVal("bam");
		String mcoordFile = cmdLine.getStringVal("mcoordFile");
		String msa = cmdLine.getStringVal("msa");
		String global = cmdLine.getStringVal("global");

		//String figure = cmdLine.getStringVal("figure");
		String tmp = cmdLine.getStringVal("tmp");
		String hours = cmdLine.getStringVal("hours");

		double scoreThreshold = cmdLine.getDoubleVal("scoreThreshold");				
		int read = cmdLine.getIntVal("read");
		int timestamp = cmdLine.getIntVal("timestamp");

		//double np = cmdLine.getDoubleVal("np");
		double il = cmdLine.getDoubleVal("il");

		boolean twodonly = cmdLine.getBooleanVal("twodonly");

		ResistanceGene paTyping = new ResistanceGene();		
		paTyping.msa = msa;
		paTyping.global = global;

		paTyping.prefix = tmp;
		paTyping.scoreThreshold = scoreThreshold;
		paTyping.twoDOnly = twodonly;
		paTyping.readNumber = read;
		DateFormat df = new SimpleDateFormat("dd/MM/yy HH:mm:ss");
		Logging.info("START : " + df.format(Calendar.getInstance().getTime()));

		if (paTyping.readNumber < 1)
			paTyping.readNumber = 1;

		paTyping.datOS = SequenceOutputStream.makeOutputStream(output);

		paTyping.getGeneClassInformation(mcoordFile);
		if (hours !=null){
			BufferedReader bf = SequenceReader.openFile(hours);
			String line = bf.readLine();//first line
			paTyping.hoursArray = new IntArray();
			paTyping.readCountArray = new IntArray();

			while ((line = bf.readLine())!= null){
				String [] tokens = line.split("\\s");
				int hrs = Integer.parseInt(tokens[0]);
				int readCount = Integer.parseInt(tokens[2]);

				paTyping.hoursArray.add(hrs);
				paTyping.readCountArray.add(readCount);	
			}
		}

		paTyping.ilThreshold = il;
		paTyping.typing(bamFile);
		paTyping.datOS.close();
	}
}
