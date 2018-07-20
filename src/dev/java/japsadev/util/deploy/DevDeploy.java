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

package japsadev.util.deploy;

import japsa.tools.bio.phylo.XMDistance2Cmd;
import japsa.util.CommandLine;
import japsa.util.deploy.Deploy;
import japsadev.tools.*;
//import japsadev.tools.ConsensusGenerateCmd;
import japsadev.tools.work.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * This class is used to deploy tools: create a makefile to generate scripts
 * 
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 */
public class DevDeploy {
	private static ArrayList<Object> tools = new ArrayList<Object>();
	static {
		//tools.add(new SampleCmd());

		tools.add(new CaptureVNTR());
		tools.add(new ResGeneGenomesCmd());
		tools.add(new GetFlankBlast());		

		tools.add(new String("Working commands"));
		tools.add(new TreeHetLinageCmd());		

		tools.add(new FilterPEConcordance());		

		tools.add(new VNTRLongReadsHmmerCmd());				

		tools.add(new RepeatPrimerCmd());		

		tools.add(new ProfileDPCmd());

		tools.add(new PlasmaAnalysisCmd());
		tools.add(new MethylationAnalysisCmd());
		tools.add(new MethylationAnalysis2Cmd());
		tools.add(new PlasmaAnalysisCrossCorrelationCmd());
		//tools.add(new RemoveNsCmd());
		tools.add(new XMDistance2Cmd());
		tools.add(new BuildMLSTTreeCmd());

		tools.add(new BuildXMTreeCmd());
		tools.add(new FixNamesTreeCmd());
		tools.add(new VNTRDepthAnalyserCmd());
		tools.add(new VNTRDepthSumCmd());
		tools.add(new FixFastqNameCmd());
		tools.add(new SampleCmd());
		tools.add(new AnalyseCaptureCmd());
		tools.add(new BreakPointAnalysisCmd());
		tools.add(new CaptureProbeDesignCmd());
		tools.add(new ConvertProbeCmd());
		tools.add(new VNTRSelectCmd());
		tools.add(new NewScarfCmd());
        tools.add(new GetCDHitCmd());

		tools.add(new GapCloserCmd());

		tools.add(new SelectReadsCmd());
		tools.add(new StructuralVariationCmd());
		
		tools.add(new VNTRClusteringCmd());
		tools.add(new VNTRClusteringHmmCmd());
        tools.add(new CheckInductionCmd());
        
        tools.add(new SelectPlasmidReadsCmd());
        tools.add(new FlankSeqsDetectorCmd());
		//new 
		//tools.add(new CheckInductionCmd());
		//

	}

	public static void main(String[] args) throws NoSuchFieldException,
	SecurityException, IOException {
		CommandLine cmdLine = new CommandLine();
		cmdLine.addString("mode", "install", "install or uinstall");
		cmdLine.addString("libs", "", "list of extenal libraries");
		cmdLine.addString("installDir", null, "the directory to install");
		cmdLine.addString("jlp", null, "Directories to libhdf5 and to jri");
		cmdLine.addString("xmx", null, "Set default maximum memory");
		cmdLine.addString("compiler", null, "Compiler version");
		cmdLine.addBoolean("version", false, "Get version and exit");	
		cmdLine.addString("server", "na", "Run on server: yes/true for yes; no/false for no");

		cmdLine.addStdHelp();// help
		/********************** Standard processing ***************************/
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		if (cmdLine.getBooleanVal("version")){
			System.out.println(Deploy.VERSION);
			System.exit(0);
		}		

		///Get command lines option
		String mode = cmdLine.getStringVal("mode");		


		Deploy.japsaPath = cmdLine.getStringVal("installDir");
		Deploy.compiler = cmdLine.getStringVal("compiler");		
		Deploy.jlp = cmdLine.getStringVal("jlp");		
		Deploy.libs = cmdLine.getStringVal("libs");		
		Deploy.maxMem = cmdLine.getStringVal("xmx");
		String serverOpt = cmdLine.getStringVal("server").toLowerCase();

		if(serverOpt.equals("yes") || serverOpt.equals("true"))
			Deploy.server = 1;		
		if(serverOpt.equals("no") || serverOpt.equals("false"))
			Deploy.server = 0;		

		Deploy.japsaJar = "japsa-dev.jar";
		/**********************************************************************/

		if ("install".equals(mode)) {		
			Deploy.setUpDirectory();

			Deploy.setUpScripts(Deploy.tools,"jsa");
			Deploy.setUpScripts(tools, "jsa.dev");	
			//System.out.println("Japsa-dev installtion complete\nFor your convenience, please add your PATH: " + Deploy.japsaPath + File.separator+"bin\n");
			System.out.println("Japsa-dev installtion complete\nFor your convenience, please add the following directory your PATH: " + Deploy.japsaPath + File.separator+"bin\n");
		} else if ("uninstall".equals(mode)) {
			if (Deploy.uninstallLibraries()){
				Deploy.uninstallScripts(Deploy.tools,"jsa");
				Deploy.uninstallScripts(tools, "jsa.dev");
			}
		}
		else if ("galaxy".equals(mode)) {
			ArrayList<Object> fullList = (ArrayList<Object>)Deploy.tools.clone();
			fullList.addAll(tools);			
			Deploy.setUpGalaxyScripts(fullList);			
		}
		else {
			System.err.println("Mode " + mode + " not recognised");
			System.err.println(cmdLine.usageString());
			System.exit(-1);
		}
	}

}
