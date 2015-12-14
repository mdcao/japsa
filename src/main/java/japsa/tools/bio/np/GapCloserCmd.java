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

package japsa.tools.bio.np;

import japsa.bio.hts.scaffold.RealtimeScaffolding;
import japsa.bio.hts.scaffold.ScaffoldGraph;
import japsa.bio.hts.scaffold.ScaffoldGraphDFS;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.IOException;

/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.np.gapcloser",
scriptDesc = "Scaffold and finish assemblies using Oxford Nanopore sequencing reads")
public class GapCloserCmd extends CommandLine{
	//static Alphabet alphabet = Alphabet.DNA();

	public GapCloserCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("bamFile", null, "Name of the bam file", true);
		addString("sequenceFile", null, "Name of the assembly file (sorted by length)",true);
		addString("output", "-", "Name of the output file, -  for stdout");		
		addInt("threshold", 0, "Margin threshold: to limit distance to the contig's ends of the alignment used in bridging."); 
		addInt("minContig", 300, "Minimum contigs length that are used in scaffolding (default 200)."); 
		//addInt("minMarker", 300, "Minimum length of markers that are used in scaffolding (default 1000)."); 
		
		addDouble("cov", 0, "Expected average coverage of Illumina, <=0 to estimate");
		addInt("qual", 1, "Minimum quality");
		addInt("support", 2, "Minimum supporting long read needed for a link between markers");
		addString("connect", null, "Name of the connection file");
		addString("stat", null, "Name of the stastistic file for Nanopore read alignment");
		addBoolean("realtime", false, "Real-time processing mode. Default is batch mode (false)");
		addInt("read", 50,  "Minimum number of reads between analyses");		
		addInt("time", 30,   "Minimum number of seconds between analyses");
		addBoolean("verbose", false, "Turn on debugging mode if true (default false)");

		addStdHelp();		
		
	} 	
	//static boolean hardClip = false;

	public static void main(String[] args) throws 
	IOException, InterruptedException {
		CommandLine cmdLine = new GapCloserCmd();		
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String output = cmdLine.getStringVal("output");
		String bamFile = cmdLine.getStringVal("bamFile");
		String sequenceFile = cmdLine.getStringVal("sequenceFile");
		int threshold = cmdLine.getIntVal("threshold");
		
		int minContig = cmdLine.getIntVal("minContig");
		if(minContig != 300)
			ScaffoldGraph.minContigLength = minContig;
		
		int minSupport = cmdLine.getIntVal("support");
		if(minSupport != 2)
			ScaffoldGraph.minSupportReads = minSupport;
		
		if(threshold != 0)
			ScaffoldGraph.marginThres = threshold;
		boolean verbose = cmdLine.getBooleanVal("verbose");
		if(verbose)
			ScaffoldGraph.verbose = verbose;
		
		double cov = cmdLine.getDoubleVal("cov");
		int qual = cmdLine.getIntVal("qual");
		boolean rt = cmdLine.getBooleanVal("realtime");
		int number = cmdLine.getIntVal("read");
		int time = cmdLine.getIntVal("time");	
		/**********************************************************************/

		SequenceOutputStream outOS = null, connectOS = null, statOS = null;
		if (cmdLine.getStringVal("connect") != null){
			connectOS = SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("connect"));			
		}
		if (cmdLine.getStringVal("stat") != null){
			statOS = SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("stat"));			
		}
		ScaffoldGraph graph;
		if(rt){
			RealtimeScaffolding rtScaffolding = new RealtimeScaffolding(sequenceFile, "-");
			rtScaffolding.scaffolding(bamFile, number, time, cov/1.6, qual);
			graph = rtScaffolding.graph;
		}
		else{
			graph = new ScaffoldGraphDFS(sequenceFile);
			if (cov <=0)
				cov = graph.estimatedCov;
			graph.makeConnections(bamFile, cov / 1.6, qual, connectOS, statOS);
			
			if (connectOS != null)
				connectOS.close();
			if(statOS != null)
				statOS.close();
	
			graph.connectBridges();
		}
		outOS = SequenceOutputStream.makeOutputStream(output);
		graph.printSequences(outOS);
		//graph.printScaffoldSequence(outOS);
		outOS.close();
		
	}
}
