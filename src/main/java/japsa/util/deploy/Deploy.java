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

import japsa.bio.ngs.AlignmentEvaluation;
import japsa.bio.ngs.BreakBam;
import japsa.bio.ngs.BreakFastQ;
import japsa.bio.ngs.CountReadInRegion;
import japsa.bio.ngs.SelectReadIntersect;
import japsa.bio.ngs.SelectReadSpan;
import japsa.bio.phylo.XMDistance;
import japsa.bio.phylo.tools.NormaliseTree;
import japsa.bio.tr.Fragment2TRV;
import japsa.bio.tr.Japsa2TR;
import japsa.bio.tr.ParseTRF;
import japsa.bio.tr.TRV2Bed;
import japsa.bio.tr.TRV2VCF;
import japsa.bio.tr.Sam2FragmentSize;
import japsa.bio.tr.SortFragmentFile;
import japsa.bio.tr.VCF2TRV;
import japsa.seq.tools.AddAnnotation;
import japsa.seq.tools.AnnotateRegions;
import japsa.seq.tools.Bed2Japsa;
import japsa.seq.tools.BreakSequenceFile;
import japsa.seq.tools.FileFormatConverter;
import japsa.seq.tools.SequenceStats;
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
	@SuppressWarnings("rawtypes")
	private static ArrayList<Class> tools = new ArrayList<Class>();
	private static String VERSION = "1.4";
	//private static String AUTHORS = "Minh Duc Cao";
	static{
		//jsa.seq.*
		tools.add(SequenceStats.class);
		tools.add(Bed2Japsa.class);
		tools.add(BreakSequenceFile.class);
		tools.add(FileFormatConverter.class);
		tools.add(AddAnnotation.class);
		tools.add(AnnotateRegions.class);
		//tools.add(MarkovCompress.class);
			

		//jsa.nsg.*
		tools.add(BreakBam.class);
		tools.add(BreakFastQ.class);
		tools.add(SelectReadIntersect.class);
		tools.add(SelectReadSpan.class);
		tools.add(AlignmentEvaluation.class);		
		tools.add(CountReadInRegion.class);
		
//		tools.add(CountReadInRegion.class);		

		//jsa.str.*
		tools.add(ParseTRF.class);		
		tools.add(Sam2FragmentSize.class);
		tools.add(SortFragmentFile.class);
		tools.add(Fragment2TRV.class);
		tools.add(TRV2VCF.class);
		tools.add(VCF2TRV.class);		
		tools.add(Japsa2TR.class);
		tools.add(TRV2Bed.class);	

		//jsa.phylo
		tools.add(XMDistance.class);
		tools.add(NormaliseTree.class);	
		
		//tools.add(.class);
		//tools.add(.class);		
	}	



	public static void main(String [] args) throws NoSuchFieldException, SecurityException, IOException{
		CommandLine cmdLine = new CommandLine();
		cmdLine.addString("mode", "install", "install or uinstall");
		cmdLine.addString("libs", "", "list of extenal libraries");
		cmdLine.addString("japsa", "japsa.jar", "name of the jar file");
		cmdLine.addString("prefix", ".", "the directory to install");
		cmdLine.addStdHelp();//help		
		
		/********************** Standard processing ***************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")){
			System.out.println(cmdLine.usage());			
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usage());
			System.exit(-1);
		}	
		/**********************************************************************/		

		String dirPath = cmdLine.getStringVal("prefix").trim();
		//Java doesnt understand ~ as home directory
		if (dirPath.startsWith("~/")){
			dirPath = System.getProperty("user.home") + dirPath.substring(1);
		}

		File dir = new File(dirPath);		
		String mode = cmdLine.getStringVal("mode");		
		String libs = cmdLine.getStringVal("libs");
		String japsa = cmdLine.getStringVal("japsa");
		File jsa = new File(dir+"/bin/jsa");			

		if ("install".equals(mode)){			
			String cp = dir.getCanonicalPath() + "/lib/japsa/"+japsa;			
			StringSeparator ss = new StringSeparator(libs,':');		
			while (ss.hasNext()){
				String l = ss.next();
				if (l.length() > 0)
					cp = cp + ":" + dir.getCanonicalPath() + "/lib/japsa/"+l;			
			}			

			PrintStream outJsa = new PrintStream(new FileOutputStream(jsa));
			outJsa.println("#!/bin/sh\n\ncat << EOF");		

			outJsa.println("Japsa: A Java Package for Statistical Sequence Analysis\n"
					+ " Version " + VERSION +  ", Built on " + (new Date()));
			outJsa.println("\nUsage:\n");

			System.out.println("\nInstalling Japsa");
			for (Class<?> tool: tools){
				Deployable annotation = tool.getAnnotation(Deployable.class);
				File file = new File(dir+"/bin/"+ annotation.scriptName());


				PrintStream out = new PrintStream(new FileOutputStream(file));
				out.println("#!/bin/sh");						


				out.println("case $MEM in '')MEM=7000m;;*);;esac\n");
				out.println("java -Xmx$MEM -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -server -classpath " + cp + ":${CLASS_PATH} " +  tool.getCanonicalName() + " $@");				
				out.close();

				Runtime.getRuntime().exec("chmod a+x "+file.getCanonicalPath());				
				System.out.println(" " +file.getCanonicalPath() + " created");							
				outJsa.printf("  %-25s  %s\n", annotation.scriptName(),annotation.scriptDesc());
			}
			outJsa.println("\nEOF");
			outJsa.close();
			System.out.println(" " +jsa.getCanonicalPath() + " created\nDone");

			Runtime.getRuntime().exec("chmod a+x "+jsa.getCanonicalPath());

		}else if ("uninstall".equals(mode)){
			//Delete all the scripts

			for (Class<?> tool: tools){				
				Deployable annotation = tool.getAnnotation(Deployable.class);
				File file = new File(dir+"/bin/"+ annotation.scriptName());				
				System.out.println("rm " + file.toPath());
				file.delete();
			}
			System.out.println("rm " + jsa.toPath());
			jsa.delete();
		}else{
			System.err.println("Mode " + mode + " not recognised");
			System.err.println(cmdLine.usage());
			System.exit(-1);
		}
	}
}
