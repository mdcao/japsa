/*****************************************************************************
 * Copyright (c) Bhuvaneswari Thirugnanasambandham, buvan.suji@gmail.com
 * All rights reserved.         *
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

/*                           Revision History                                
 *
 ****************************************************************************/

package japsadev.tools;




import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import japsadev.bio.hts.clustering.KmeanClustering;



/**
 * @author Bhuvaneswari Thirugnanasambandham
 *
 */
@Deployable(
		scriptName = "jsa.dev.kcluster",
		scriptDesc = "Clustering reads"
		)
public class VNTRClusteringCmd extends CommandLine{
	//CommandLine cmdLine;
	public VNTRClusteringCmd(){
		
		
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", null,  "Input file");
		addString("output", "-",  "Output file");

		addStdHelp();
		
		
		
		
	}
	public static void main(String [] args) throws Exception{
		CommandLine cmdLine = new VNTRClusteringCmd();
		args = cmdLine.stdParseLine(args);	
		/*String[] s = args;
		System.out.println(s);*/
		
		


		/**********************************************************************/

		String input = cmdLine.getStringVal("input");
		String output= cmdLine.getStringVal("output");
		////YOUR CODE GOES HERE		
		/*System.out.println("Hello world input is " + input);
		System.out.println("Testing this statement");*/
		////////////
		
		KmeanClustering cluster = new KmeanClustering();
		cluster.Clustering();
		
		
	}
}
