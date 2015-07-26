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
 * File: Deploy.java
 * 15/11/2013 - Minh Duc Cao: Created
 *
 ****************************************************************************/

package japsa.util.deploy;

import japsa.bio.alignment.AlignmentEM;
import japsa.bio.hts.BreakBam;
import japsa.bio.hts.CountReadInRegion;
import japsa.bio.hts.FastQRMEmptyRead;
import japsa.bio.hts.FastQTrim;
import japsa.bio.hts.HTSAlignmentParam;
import japsa.bio.hts.HTSErrorAnalysis;
import japsa.bio.hts.SelectReadIntersect;
import japsa.bio.hts.SelectReadSpan;
import japsa.bio.np.GeneStrainTyping;
import japsa.bio.np.MLSTStrainTyping;
import japsa.bio.np.ResistanceGene;
import japsa.bio.np.SpeciesMixtureTyping;
import japsa.bio.phylo.XMDistance;
import japsa.bio.phylo.XMDistance2;
import japsa.bio.phylo.tools.NormaliseTree;
import japsa.bio.sim.SimHTSWithFSM;
import japsa.bio.sim.SimProbFSM;
import japsa.bio.tr.Fragment2TRV;
import japsa.bio.tr.Japsa2TR;
import japsa.bio.tr.ParseTRF;
import japsa.bio.tr.TRV2Bed;
import japsa.bio.tr.TRV2VCF;
import japsa.bio.tr.Sam2FragmentSize;
import japsa.bio.tr.SortFragmentFile;
import japsa.bio.tr.VCF2TRV;
import japsa.bio.tr.VNTRDepth;
import japsa.seq.nanopore.NanoporeReadFilter;
//import japsa.seq.nanopore.NanoporeReader;
import japsa.seq.nanopore.NanoporeReaderStream;
import japsa.seq.tools.AddAnnotation;
import japsa.seq.tools.AnnotateRegions;
import japsa.seq.tools.AnnotateVCF;
import japsa.seq.tools.Bed2Japsa;
import japsa.seq.tools.ExtractGeneSequence;
import japsa.seq.tools.FileFormatConverter;

import japsa.tools.seq.JoinSequenceFile;
import japsa.tools.seq.SequenceExtract;
import japsa.tools.seq.SequenceReverseComplement;
import japsa.tools.seq.SequenceSort;
import japsa.tools.seq.SequenceStats;
import japsa.tools.seq.SplitSequenceFile;
import japsa.tools.util.StreamServer;
import japsa.tools.xm.ExpertModelDriver;
import japsa.util.CommandLine;
import japsa.util.StringSeparator;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Date;

/**
 * This class is used to deploy tools: create a makefile to generate scripts
 * 
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 */
public class Deploy {	
	private static ArrayList<Object> tools = new ArrayList<Object>();
	public static String VERSION = "1.5-8a";
	//private static String AUTHORS = "Minh Duc Cao";
	static{
		//jsa.seq.*
		tools.add("Sequence manipulation tools:");
		tools.add(SequenceStats.class);
		tools.add(SequenceSort.class);
		tools.add(SequenceExtract.class);
		tools.add(SequenceReverseComplement.class);		
		tools.add(Bed2Japsa.class);
		tools.add(SplitSequenceFile.class);
		tools.add(JoinSequenceFile.class);	
		tools.add(FileFormatConverter.class);
		tools.add(AddAnnotation.class);
		tools.add(AnnotateRegions.class);
		tools.add(AnnotateVCF.class);
		tools.add(ExtractGeneSequence.class);
		tools.add(AlignmentEM.class);
		
		//tools.add(MarkovCompress.class);

		//jsa.hts.*
		tools.add("HTS analysis tools:");
		tools.add(FastQTrim.class);
		tools.add(FastQRMEmptyRead.class);		
		tools.add(BreakBam.class);
		tools.add(SelectReadIntersect.class);
		tools.add(SelectReadSpan.class);				
		tools.add(CountReadInRegion.class);
		tools.add(HTSAlignmentParam.class);
		tools.add(HTSErrorAnalysis.class);
		

		//jsa.np.
		//tools.add(NanoporeReader.class);
		tools.add("Oxford Nanopore sequencing analysis tools:");
		tools.add(NanoporeReaderStream.class);
		tools.add(NanoporeReadFilter.class);		
		tools.add(SpeciesMixtureTyping.class);		
		tools.add(GeneStrainTyping.class);
		tools.add(MLSTStrainTyping.class);
		tools.add(ResistanceGene.class);		

		//jsa.trv.*
		tools.add("Tandem repeat variation analysis tools:");
		tools.add(ParseTRF.class);		
		tools.add(Sam2FragmentSize.class);
		tools.add(SortFragmentFile.class);
		tools.add(Fragment2TRV.class);
		tools.add(TRV2VCF.class);
		tools.add(VCF2TRV.class);		
		tools.add(Japsa2TR.class);
		tools.add(TRV2Bed.class);	
		tools.add(VNTRDepth.class);
		
		tools.add("Utilities:");
		tools.add(StreamServer.class);
		//tools.add(VNTRDepth.class);

		//jsa.phylo		
		tools.add("Phylogenetics analysis tools:");
		tools.add(XMDistance.class);
		tools.add(XMDistance2.class);
		tools.add(NormaliseTree.class);	
		
		//jsa.sim
		tools.add("Alignment with Finite State Machines");
		tools.add(SimProbFSM.class);		
		tools.add(SimHTSWithFSM.class);
		
		//jsa.xm
		tools.add("Export Model compression");
		tools.add(ExpertModelDriver.class);
		//tools.add(.class);		
	}	



	public static void main(String[] args) throws NoSuchFieldException,
	SecurityException, IOException {
		CommandLine cmdLine = new CommandLine();
		cmdLine.addString("mode", "install", "install or uinstall");
		cmdLine.addString("libs", "", "list of extenal libraries");
		cmdLine.addString("japsa", "japsa.jar", "name of the jar file");
		cmdLine.addString("prefix", ".", "the directory to install");
		cmdLine.addString("jlp", "", "java.library.path");
		cmdLine.addBoolean("version", false, "Get version and exit");
		
		cmdLine.addStdHelp();// help

		/********************** Standard processing ***************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")) {
			System.out.println(cmdLine.usage());
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usage());
			System.exit(-1);
		}
		/**********************************************************************/
		if (cmdLine.getBooleanVal("version")){
			System.out.println(VERSION);
			System.exit(0);
		}

		String dirPath = cmdLine.getStringVal("prefix").trim();
		String jlp = cmdLine.getStringVal("jlp");
		// Java doesnt understand ~ as home directory
		if (dirPath.startsWith("~/")) {
			dirPath = System.getProperty("user.home") + dirPath.substring(1);
		}

		File dir = new File(dirPath);
		String mode = cmdLine.getStringVal("mode");
		String libs = cmdLine.getStringVal("libs");
		String japsa = cmdLine.getStringVal("japsa");
		File jsa = new File(dir + "/bin/jsa");

		if ("install".equals(mode)) {
			String cp = dir.getCanonicalPath() + "/lib/japsa/" + japsa;
			StringSeparator ss = new StringSeparator(libs, ':');
			while (ss.hasNext()) {
				String l = ss.next();
				if (l.length() > 0)
					cp = cp + ":" + dir.getCanonicalPath() + "/lib/japsa/" + l;
			}

			PrintStream outJsa = new PrintStream(new FileOutputStream(jsa));
			outJsa.println("#!/bin/sh\n\ncat << EOF");

			outJsa.println("Japsa: A Java Package for Statistical Sequence Analysis\n"
					+ " Version " + VERSION + ", Built on " + (new Date()));
			outJsa.println("\nList of tools:");

			String JAVA_COMMAND ="java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -server"; 
			if (jlp.length() > 0){
				JAVA_COMMAND += " -Djava.library.path="+jlp;
			}

			System.out.println("\nInstalling Japsa");
			
			for (Object obj : tools) {		
				
				if (!(obj instanceof Class<?>)){
					outJsa.printf("\n%s\n",obj);
					continue;
				}				
				Class<?> tool = (Class<?>) obj;
				
				Deployable annotation = tool.getAnnotation(Deployable.class);
				File file = new File(dir + "/bin/" + annotation.scriptName());

				PrintStream out = new PrintStream(new FileOutputStream(file));
				out.println("#!/bin/sh");

				out.println("case $JSA_MEM in '')JSA_MEM=7000m;;*);;esac\n\n");

				out.println("case $JSA_CP in\n  '')JSA_CP="
						+ cp
						+ ";;\n  *)echo \"[INFO] Use ${JSA_CP} as path \" 1>&2;;\nesac\n\n");
				
				out.println("JSA_CMD=\"`basename $0` $@\"\n");
				
				out.println(JAVA_COMMAND + " -classpath ${JSA_CP} "
						+ tool.getCanonicalName() + " \"$@\"");
				out.close();

				Runtime.getRuntime().exec(
						"chmod a+x " + file.getCanonicalPath());
				System.out.println(" " + file.getCanonicalPath() + " created");
				outJsa.printf("  %-23s  %s\n", annotation.scriptName(),
						annotation.scriptDesc());
			}
			outJsa.println("\nEOF");
			outJsa.close();
			System.out.println(" " + jsa.getCanonicalPath() + " created\nDone");

			Runtime.getRuntime().exec("chmod a+x " + jsa.getCanonicalPath());

		} else if ("uninstall".equals(mode)) {
			// Delete all the scripts
			
			for (Object obj : tools) {				
				if (!(obj instanceof Class<?>)){			
					continue;
				}				
				Class<?> tool = (Class<?>) obj;
				Deployable annotation = tool.getAnnotation(Deployable.class);
				File file = new File(dir + "/bin/" + annotation.scriptName());
				System.out.println("rm " + file.getCanonicalPath());				
				file.delete();
			}
			System.out.println("rm " + jsa.getCanonicalPath());
			jsa.delete();
		} else {
			System.err.println("Mode " + mode + " not recognised");
			System.err.println(cmdLine.usage());
			System.exit(-1);
		}
	}
}
