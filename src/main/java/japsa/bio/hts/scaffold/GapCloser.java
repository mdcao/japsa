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
 * 30/06/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.hts.scaffold;

import japsa.seq.Alphabet;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.IOException;

/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.dev.gapcloser", scriptDesc = "Gap closer")
public class GapCloser {
	static Alphabet alphabet = Alphabet.DNA();
	//static boolean hardClip = false;

	public static void main(String[] args) throws IOException,
	InterruptedException {
		/*********************** Setting up script ****************************/
		Deployable annotation = GapCloser.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options]",
				annotation.scriptDesc());

		cmdLine.addString("bamFile", null, "Name of the bam file", true);
		cmdLine.addString("sequenceFile", null, "Name of the assembly file (sorted by length)",true);
		cmdLine.addString("output", "-",
				"Name of the output file, -  for stdout");		
		cmdLine.addInt("threshold", 500, "Threshold");
		cmdLine.addDouble("cov", 0, "Expected average coverage of Illumina, <=0 to estimate");
		cmdLine.addInt("qual", 1, "Minimum quality");
		cmdLine.addString("connect", null, "Name of the connection file");
		

		//String [] processArgs = 
				cmdLine.stdParseLine_old(args);		
		/**********************************************************************/
		String output = cmdLine.getStringVal("output");
		String bamFile = cmdLine.getStringVal("bamFile");
		String sequenceFile = cmdLine.getStringVal("sequenceFile");
		int threshold = cmdLine.getIntVal("threshold");
		double cov = cmdLine.getDoubleVal("cov");
		int qual = cmdLine.getIntVal("qual");
		/**********************************************************************/

		SequenceOutputStream outOS = null;
		if (cmdLine.getStringVal("connect") != null){
			outOS = SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("connect"));			
		}
		
		ScaffoldGraph graph = new ScaffoldGraphDFS(sequenceFile);
		if (cov <=0)
			cov = graph.estimatedCov;

		graph.makeConnections(bamFile, cov / 1.6,  cov * 1.46, threshold, qual,outOS);
		if (outOS != null)
			outOS.close();
		
		graph.connectBridges();
		
		outOS = SequenceOutputStream.makeOutputStream(output);
		graph.printSequences(outOS);
		//graph.printScaffoldSequence(outOS);
		outOS.close();		
	}
}
