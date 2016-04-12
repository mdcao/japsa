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
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.IOException;

/**
 * @author sonnguyen, minhduc
 * 
 */
@Deployable(
		scriptName = "jsa.np.gapcloser",
		scriptDesc = "Scaffold and finish assemblies using Oxford Nanopore sequencing reads",
		seeAlso = "jsa.np.f5reader, jsa.util.streamServer, jsa.util.streamClient"
)
public class GapCloserCmd extends CommandLine{

	public GapCloserCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("bamFile", null, "Name of the bam file", true);
		addString("sequenceFile", null, "Name of the assembly file (sorted by length)",true);
		addString("prefix", "out", "Prefix for the output files, default is out.*");	
		
		addString("resistGene", null , "Realtime annotation: name of antibiotic resistance gene fasta file");
		addString("insertSeq", null , "Realtime annotation: name of IS fasta file");
		addString("oriRep", null, "Realtime annotation: name of fasta file containing possible origin of replication");
		
		addInt("marginThres", 1000, "Margin threshold: to limit distance to the contig's ends of the alignment used in bridging."); 
		addInt("minContig", 300, "Minimum contigs length that are used in scaffolding (default 300)."); 
		addInt("maxRepeat", 7500, "Maximum length of repeat in considering species (default 7500 for bacteria)."); 
		
		addDouble("cov", 0, "Expected average coverage of Illumina, <=0 to estimate");
		addInt("qual", 1, "Minimum quality");
		addInt("support", 1, "Minimum supporting long read needed for a link between markers");
	
		addBoolean("realtime", false, "Process in real-time mode. Default is batch mode (false)");
		addInt("read", 500,  "Minimum number of reads between analyses");		
		addInt("time", 300,   "Minimum number of seconds between analyses");
		addBoolean("verbose", false, "Turn on debugging mode");

		addStdHelp();		
		
	} 	
	//static boolean hardClip = false;

	public static void main(String[] args) throws 
	IOException, InterruptedException {
		CommandLine cmdLine = new GapCloserCmd();		
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String prefix = cmdLine.getStringVal("prefix");
		String bamFile = cmdLine.getStringVal("bamFile");
		String sequenceFile = cmdLine.getStringVal("sequenceFile"),
				resistFile = cmdLine.getStringVal("resistGene"),
				isFile = cmdLine.getStringVal("insertSeq"),
				oriFile = cmdLine.getStringVal("oriRep");
		
		int 	marginThres = cmdLine.getIntVal("marginThres"),
				minContig = cmdLine.getIntVal("minContig"),
				minSupport = cmdLine.getIntVal("support"),
				maxRepeat = cmdLine.getIntVal("maxRepeat");
		if(marginThres < 0)
			Logging.exit("Marginal threshold must not be negative", 1);			
		if(minContig <= 0)
			Logging.exit("Minimum contig length has to be positive", 1);
		if(minSupport <= 0)
			Logging.exit("Minimum supporting reads has to be positive", 1);
		if(maxRepeat <= 0)
			Logging.exit("Maximal possible repeat length has to be positive", 1);
			
		
		ScaffoldGraph.minContigLength = minContig;
		ScaffoldGraph.minSupportReads = minSupport;	
		ScaffoldGraph.maxRepeatLength = maxRepeat;		
		ScaffoldGraph.marginThres = marginThres;
		ScaffoldGraph.verbose = cmdLine.getBooleanVal("verbose");
			
		
		double cov = cmdLine.getDoubleVal("cov");
		int qual = cmdLine.getIntVal("qual");
		if(qual < 0)
			Logging.exit("Phred score of quality has to be positive", 1);
		
		int number = cmdLine.getIntVal("read"),
			time = cmdLine.getIntVal("time");
		
		if(number <= 0)
			Logging.exit("Number of reads has to be positive", 1);			
		if(time < 0)
			Logging.exit("Sleeping time must not be negative", 1);	
		/**********************************************************************/

		ScaffoldGraph graph;
		boolean rt = cmdLine.getBooleanVal("realtime");
		if(rt){
			RealtimeScaffolding rtScaffolding = new RealtimeScaffolding(sequenceFile, resistFile, isFile, oriFile, "-");
			rtScaffolding.scaffolding(bamFile, number, time, cov/1.6, qual);
			graph = rtScaffolding.graph;
		}
		else{
			graph = new ScaffoldGraphDFS(sequenceFile, resistFile, isFile, oriFile);
			if (cov <=0)
				cov = graph.estimatedCov;
			graph.makeConnections(bamFile, cov / 1.6, qual);
	
			graph.connectBridges();
		}
		if(prefix != null)
			graph.prefix = prefix;
		graph.printSequences();
	}
}
/*RST*
-------------------------------------------------------------------------
*npScaffolder*: real-time scaffolder using SPAdes contigs and Nanopore sequencing reads
-------------------------------------------------------------------------
*npScaffolder* (jsa.np.gapcloser) is a program that connect contigs from a draft genomes 
to generate sequences that are closer to finish. These pipelines can run on a single laptop
for microbial datasets. In real-time mode, it can be integrated with simple structural 
analyses such as gene ordering, plasmid forming.
npScaffolder is included in the `Japsa package <http://mdcao.github.io/japsa/>`_.
<usage>
~~~~~~~~~~~~~~
Usage examples
~~~~~~~~~~~~~~
A summary of *npScaffolder* usage can be obtained by invoking the --help option::
   	jsa.np.gapcloser --help
Input
======
*npScaffolder* takes two files as required input::
	jsa.np.gapcloser -s <*draft*> -b <*bam*>
	
<*draft*> input is the FASTA file containing the pre-assemblies. Normally this 
is the output from running SPAdes on Illumina MiSeq paired end reads.
<*bam*> contains SAM/BAM formated alignments between <*draft*> file and <*nanopore*> 
FASTA/FASTQ file of long read data. We use BWA-MEM as the recommended aligner 
with the fixed parameter set as follow::
	bwa mem -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y <*draft*> <*nanopore*> > <*bam*>
	
Output
=======
*npScaffolder* output is specified by *-prefix* option. The default prefix is \'out\'.
Normally the tool generate two files: *prefix*.fin.fasta and *prefix*.fin.japsa which 
indicate the result scaffolders in FASTA and JAPSA format.
In realtime mode, if any annotation analysis is enabled, a file named 
*prefix*.anno.japsa is generated instead. This file contains features detected after
scaffolding.
Real-time scaffolding
===============
To run *npScaffolder* in streaming mode::
   	jsa.np.gapcloser -realtime [options]
In this mode, the <*bam*> file will be processed block by block. The size of block 
(number of BAM/SAM records) can be manipulated through option *-read* and *-time*.
The idea of streaming mode is when the input <*nanopore*> file is retrieved in stream.
npReader is the module that provides such data from fast5 files returned from the real-time
base-calling cloud service Metrichor. Ones can run::
jsa.np.f5reader -realtime -folder c:\Downloads\ -fail -output - | \
bwa mem -t 10 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y -K 3000 <*draft*> - 2> /dev/null | \ 
jsa.np.gapcloser --realtime -b - -seq <*draft*> > log.out 2>&1
or if you have the whole set of Nanopore long reads already and want to emulate the 
streaming mode::
jsa.np.timeEmulate -s 100 -i <*nanopore*> -output - | \
bwa mem -t 10 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y -K 3000 <*draft*> - 2> /dev/null | \ 
jsa.np.gapcloser --realtime -b - -seq <*draft*> > log.out 2>&1
Note that jsa.np.timeEmulate based on the field *timeStamp* located in the read name line to
decide the order of streaming data. So if your input <*nanopore*> already contains the field,
you have to sort it::
jsa.seq.sort -i <*nanopore*> -o <*nanopore-sorted*> -sortKey=timeStamp
or if your file does not have the *timeStamp* data yet, you can manually make ones. For example::
cat <*nanopore*> |awk 'BEGIN{time=0.0}NR%4==1{printf "%s timeStamp=%.2f\n", $0, time; time++}NR%4!=1{print}'
> <*nanopore-with-time*> 
Real-time annotation
====================
The tool includes usecase for streaming annotation. Ones can provides database of antibiotic
resistance genes and/or Origin of Replication in FASTA format for the analysis of gene ordering
and/or plasmid identifying respectively::
jsa.np.timeEmulate -s 100 -i <*nanopore*> -output - | \
bwa mem -t 10 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y -K 3000 <*draft*> - 2> /dev/null | \ 
jsa.np.gapcloser --realtime -b - -seq <*draft*> -resistGene <*resistDB*> -oriRep <*origDB*> > log.out 2>&1
*RST*/
