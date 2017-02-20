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

package japsadev.tools;

import japsadev.bio.hts.scaffold.ContigBridge;
import japsadev.bio.hts.scaffold.RealtimeScaffolding;
import japsadev.bio.hts.scaffold.ScaffoldGraph;
import japsadev.bio.hts.scaffold.ScaffoldGraphDFS;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;

/**
 * @author sonnguyen, minhduc
 * 
 */
@Deployable(
		scriptName = "jsa.dev.npscarf",
		scriptDesc = "Experimental Scaffold and finish assemblies using Oxford Nanopore sequencing reads",
		seeAlso = "jsa.np.npreader, jsa.util.streamServer, jsa.util.streamClient"
		)
public class GapCloserCmd extends CommandLine{

	public GapCloserCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("seqFile", null, "Name of the assembly file (sorted by length)",true);
		//addString("bamFile", null, "Name of the bam file", true);

		addString("input", "-", "Name of the input file, - for stdin", true);
		addString("format", "fastq/fasta", "format of the input fastq/fasta or sam/bam");
		addString("bwaExe", "bwa", "Path to bwa");
		addInt("bwaThread", 4, "Theads used by bwa");
		addBoolean("long", false, "Whether report all sequences, including short/repeat contigs (default) or only long/unique/completed sequences.");
		addBoolean("eukaryotic", false, "Whether eukaryotic or bacterial (default) genomes");

		addString("assembler", "spades", "Name of the assembler used for Illumina assembly: SPAdes (default) or ABySS.");
		addString("graphDir", null, "Name of the output folder by SPAdes/ABySS: assembly graph and paths will be used for better gap-filling.");
		addString("prefix", "out", "Prefix for the output files");	
		
		addString("genes", null , "Realtime annotation: name of annotated genes in GFF 3.0 format");
		addString("resistGene", null , "Realtime annotation: name of antibiotic resistance gene fasta file");
		addString("insertSeq", null , "Realtime annotation: name of IS fasta file");
		addString("oriRep", null, "Realtime annotation: name of fasta file containing possible origin of replication");
		addInt("minContig", 300, "Minimum contigs length that are used in scaffolding."); 
		addInt("maxRepeat", 7500, "Maximum length of repeat in considering species."); 

		addDouble("cov", 0, "Expected average coverage of Illumina, <=0 to estimate");
		addInt("qual", 1, "Minimum quality");
		addInt("support", 1, "Minimum supporting long read needed for a link between markers");

		addBoolean("realtime", false, "Process in real-time mode. Default is batch mode (false)");
		addInt("read", 50,  "Minimum number of reads between analyses");		
		addInt("time", 10,   "Minimum number of seconds between analyses");
		addBoolean("verbose", false, "Turn on debugging mode");

		addStdHelp();		

	} 	
	//static boolean hardClip = false;

	public static void main(String[] args) throws 
	IOException, InterruptedException {
		CommandLine cmdLine = new GapCloserCmd();		
		args = cmdLine.stdParseLine(args);

		/***********************************************************************/
		String prefix = cmdLine.getStringVal("prefix");
		//String bamFile = cmdLine.getStringVal("bamFile");

		String input = cmdLine.getStringVal("input");
		String bwaExe = cmdLine.getStringVal("bwaExe");
		int bwaThread = cmdLine.getIntVal("bwaThread");
		String format = cmdLine.getStringVal("format").toLowerCase().trim();
		


		String sequenceFile = cmdLine.getStringVal("seqFile"),
				graphFolder = cmdLine.getStringVal("graphDir"),

				genesFile = cmdLine.getStringVal("genes"),
				resistFile = cmdLine.getStringVal("resistGene"),
				isFile = cmdLine.getStringVal("insertSeq"),
				oriFile = cmdLine.getStringVal("oriRep");


		String assembler = cmdLine.getStringVal("assembler").toLowerCase();
		//TODO: driver for ABySS
		if(assembler.equals("abyss")){
			ScaffoldGraph.assembler=0b01;
			//graphFile = ;
			//pathFile = ;
		}
		
		if (format.startsWith("fastq") ||
				format.startsWith("fasta") ||
				format.startsWith("fq") ||
				format.startsWith("fa")){
			try{
				ProcessBuilder pb = new ProcessBuilder(bwaExe).redirectErrorStream(true);
				Process process =  pb.start();
				BufferedReader bf = SequenceReader.openFile(process.getInputStream());
				String line;
				String version = "";
				while ((line = bf.readLine())!=null){
					if (line.startsWith("Version: ")){
						version = line.substring(9).trim();
						break;//while
					}				
				}	
				bf.close();
				if (version.length() == 0){
					System.err.println(bwaExe + " is not the rith path to bwa. bwa is required");
					System.exit(1);
				}else{
					if (!version.startsWith("0.7.1")){
						System.err.println(" Require bwa of 0.7.11 or above");
						System.exit(1);
					}
				}
			}catch (IOException e){
				System.err.println(e.getMessage());
				System.exit(1);
			}

		}else if (format.startsWith("sam") || format.startsWith("bam")){
			// no problem
		}else{
			System.err.println("I dont understand format " + format);
			System.exit(1);
		}

		File 	graphFile = new File(graphFolder+"/assembly_graph.fastg"),
				pathFile = new File(graphFolder+"/contigs.paths");
		
		if(graphFolder !=null){
			if(ScaffoldGraph.assembler==0b00  && graphFile.exists() && pathFile.exists()){
				Logging.info("===> Use assembly graph from short-read assembler: SPAdes!");
				
			}else if (ScaffoldGraph.assembler==0b01){
				File f = new File(graphFolder);
				File[] matchingFiles = f.listFiles(new FilenameFilter() {
				    public boolean accept(File dir, String name) {
				        return name.endsWith("contigs.dot");
				    }
				});
				if(matchingFiles.length != 1){
					Logging.error("Failed to looking for an unique *-contigs.dot file in " + graphFolder + " . Proceeding without assembly graph...");
					graphFolder=null;
				} else{
					graphFile=matchingFiles[0];
					Logging.info("===> Use assembly graph from short-read assembler ABySS: " + graphFile);
				}
				
			}
			
		}
		else{
			Logging.warn("Not found any legal assembly output folder, assembly graph thus not included!");
			graphFolder=null;
		}


		int 	//marginThres = cmdLine.getIntVal("marginThres"),
		minContig = cmdLine.getIntVal("minContig"),
		minSupport = cmdLine.getIntVal("support"),
		maxRepeat = cmdLine.getIntVal("maxRepeat");
		//if(marginThres < 0)
		//	Logging.exit("Marginal threshold must not be negative", 1);			
		if(minContig <= 0)
			Logging.exit("Minimum contig length has to be positive", 1);
		if(minSupport <= 0)
			Logging.exit("Minimum supporting reads has to be positive", 1);
		if(maxRepeat <= 0)
			Logging.exit("Maximal possible repeat length has to be positive", 1);


		ScaffoldGraph.verbose = cmdLine.getBooleanVal("verbose");
		ScaffoldGraph.reportAll = !cmdLine.getBooleanVal("long");
		ScaffoldGraph.eukaryotic = cmdLine.getBooleanVal("eukaryotic");
		
		ScaffoldGraph.minContigLength = minContig;
		ScaffoldGraph.minSupportReads = minSupport;	
		ScaffoldGraph.maxRepeatLength = ScaffoldGraph.eukaryotic?Math.max(maxRepeat,10000):Math.max(maxRepeat, 7500);		
		//ScaffoldGraph.marginThres = marginThres;


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
		ContigBridge.relaxFilling();
		if(rt){
			RealtimeScaffolding rtScaffolding = new RealtimeScaffolding(sequenceFile, genesFile, resistFile, isFile, oriFile, "-");

			graph = rtScaffolding.graph;
			if(prefix != null)
				graph.prefix = prefix;
			if(graphFolder!=null){
				synchronized(graph){
					if(ScaffoldGraph.assembler==0b00)
						graph.readMore(graphFile.getAbsolutePath(),pathFile.getAbsolutePath());
					else if(ScaffoldGraph.assembler==0b01)
						graph.readMore(graphFile.getAbsolutePath(),"");
									
				}
			}
			if (cov <=0)
				cov = ScaffoldGraph.estimatedCov;

			rtScaffolding.scaffolding2(input, number, time, cov/1.6, qual, format, bwaExe, bwaThread, sequenceFile);

		}
		else{
			graph = new ScaffoldGraphDFS(sequenceFile, genesFile, resistFile, isFile, oriFile);
			if(graphFolder!=null){
				if(ScaffoldGraph.assembler==0b00)
					graph.readMore(graphFile.getAbsolutePath(),pathFile.getAbsolutePath());
				else if(ScaffoldGraph.assembler==0b01)
					graph.readMore(graphFile.getAbsolutePath(),"");
			}

			if (cov <=0)
				cov = ScaffoldGraph.estimatedCov;

			graph.makeConnections2(input, cov / 1.6, qual, format, bwaExe, bwaThread, sequenceFile);

			graph.connectBridges();
			if(prefix != null)
				graph.prefix = prefix;

			ContigBridge.forceFilling();
			graph.printSequences();
		}

	}
}

/*RST*
---------------------------------------------------------------------------------------
 *npScaffolder*: real-time scaffolder using SPAdes contigs and Nanopore sequencing reads
---------------------------------------------------------------------------------------

 *npScaffolder* (jsa.np.npscarf) is a program that connect contigs from a draft genomes 
to generate sequences that are closer to finish. These pipelines can run on a single laptop
for microbial datasets. In real-time mode, it can be integrated with simple structural 
analyses such as gene ordering, plasmid forming.

npScaffolder is included in the `Japsa package <http://mdcao.github.io/japsa/>`_.

<usage>

~~~~~~~~~~~~~~
Usage examples
~~~~~~~~~~~~~~

A summary of *npScarf* usage can be obtained by invoking the --help option::

    jsa.np.npscarf --help

Input
=====
 *npScarf* takes two files as required input::

	jsa.np.npscarf -seq <*draft*> -input <*nanopore*>

<*draft*> input is the FASTA file containing the pre-assemblies. We support outputs 
from running SPAdes (for small genomes, e.g. microbial) or ABySS (for large eukaryotic 
genomes) on Illumina MiSeq paired end reads. 
The assembler used to generate the pre-assemblies can be specified using option:

 	--assembler <*assembler*> 

The default <*assembler*> is SPAdes. 
Input from another non-supported short-read assemblers could be used with your own risk.

<*nanopore*> is either the long reads in FASTA/FASTQ file or SAM/BAM formated alignments 
between them to <*draft*> file. We use BWA-MEM as the recommended aligner 
with the fixed parameter set as follow::

	bwa mem -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y <*draft*> <*nanopore*> > <*bam*>
	
The input file format is specified by option --format. The default is FASTA/FASTQ in which 
the path to BWA version 0.7.11 or newer is required. If SAM/BAM is provided as input instead,
then do not worry about the aligner. 

	Note: Remember to always *INDEXING* the reference before running BWA::
	
	bwa index <*draft*>
	
	Missing this step would break down the whole pipeline.

Output
=======
 *npScarf* output is specified by *-prefix* option. The default prefix is \'out\'.
Normally the tool generates two files: *prefix*.fin.fasta and *prefix*.fin.japsa which 
indicate the result contigs in FASTA and JAPSA format.

In realtime mode, if any annotation analysis is enabled, a file named 
 *prefix*.anno.japsa is generated instead. This file contains features detected after
scaffolding.

Real-time scaffolding
=====================
To run *npScarf* in streaming mode::

   	jsa.np.npscarf -realtime [options]

In this mode, the <*bam*> file will be processed block by block. The size of block 
(number of BAM/SAM records) can be manipulated through option *-read* and *-time*.

The idea of streaming mode is when the input <*nanopore*> file is retrieved in stream.
npReader is the module that provides such data from fast5 files returned from the real-time
base-calling cloud service Metrichor. Ones can run::

    jsa.np.npreader -realtime -folder c:\Downloads\ -fail -output - | \
      jsa.np.npscarf --realtime -bwaExe=<path_to_BWA> -bwaThread=10 -input - -seq <*draft*> > log.out 2>&1
    
For the same purpose, you can also invoke BWA-MEM explicitly as in the old version of *npScarf*,
In this case, option --format=SAM must be presented as follow:
      
    jsa.np.npreader -realtime -folder c:\Downloads\ -fail -output - | \
      bwa mem -t 10 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y -K 3000 <*draft*> - 2> /dev/null | \ 
      jsa.np.npscarf --realtime -input - -format=SAM -seq <*draft*> > log.out 2>&1

or if you have the whole set of Nanopore long reads already and want to emulate the 
streaming mode::

    jsa.np.timeEmulate -s 100 -i <*nanopore*> -output - | \
      jsa.np.npscarf --realtime -bwaExe=<path_to_BWA> -bwaThread=10 -input - -seq <*draft*> > log.out 2>&1

Note that jsa.np.timeEmulate based on the field *timestamp* located in the read name line to
decide the order of streaming data. So if your input <*nanopore*> already contains the field,
you have to sort it::

    jsa.seq.sort -i <*nanopore*> -o <*nanopore-sorted*> -sortKey=timestamp

or if your file does not have the *timestamp* data yet, you can manually make ones. For example::

    cat <*nanopore*> | \
       awk 'BEGIN{time=0.0}NR%4==1{printf "%s timestamp=%.2f\n", $0, time; time++}NR%4!=1{print}' \
       > <*nanopore-with-time*> 

Real-time annotation
====================
The tool includes usecase for streaming annotation. Ones can provides database of antibiotic
resistance genes and/or Origin of Replication in FASTA format for the analysis of gene ordering
and/or plasmid identifying respectively::

    jsa.np.timeEmulate -s 100 -i <*nanopore*> -output - | \ 
      jsa.np.npscarf --realtime -bwaExe=<path_to_bwa> -input - -seq <*draft*> -resistGene <*resistDB*> -oriRep <*origDB*> > log.out 2>&1

Assembly graph
==============
 *npScarf* can read the assembly graph info from SPAdes/AbySS to make the results more precise.
The results might be slightly deviate from the old version in term of number of final contigs::

    jsa.np.npscarf --graphDir=<output_directory> <options...>

where output_directory indicates the result folder of SPAdes/ABySS, containing files such as contigs.fasta, 
contigs.paths, assembly_graph.fastg in case of SPAdes and *-contigs.fa, *-contigs.dot if from ABySS.
 *RST*/
