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
package japsa.tools.bio.np;


import java.io.IOException;

import japsa.bio.np.RealtimeStrainTyping;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */

@Deployable(
	scriptName = "jsa.np.rtStrainTyping", 
	scriptDesc = "Realtime strain typing using Nanopore sequencing data")
public class RealtimeStrainTypingCmd extends CommandLine{	
	public RealtimeStrainTypingCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("output", "output.dat",  "Output file");
		addString("profile", null,  "Output file containing gene profile of all strains",true);		
		addString("bamFile", null,  "The bam file",true);
		addString("geneFile", null,  "The gene file",true);

		addInt("top", 10,  "The number of top strains");
		addString("hours", null,  "The file containging hours against yields, if set will output acording to tiime");
		addInt("timestamp", 0,  "Timestamp to check, if <=0 then use read number instead");
		addInt("read", 500,  "Number of reads before a typing, NA if timestamp is set");

		addBoolean("twodonly", false,  "Use only two dimentional reads");

		addStdHelp();		
	} 

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException, InterruptedException{
		CommandLine cmdLine = new RealtimeStrainTypingCmd();		
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/

		String output = cmdLine.getStringVal("output");
		String profile = cmdLine.getStringVal("profile");				
		String bamFile = cmdLine.getStringVal("bamFile");
		String geneFile = cmdLine.getStringVal("geneFile");
		
		int top = cmdLine.getIntVal("top");		

		RealtimeStrainTyping paTyping = new RealtimeStrainTyping(geneFile, profile, output);
		paTyping.typing(bamFile,  top);
		paTyping.close();
		
	}


}
