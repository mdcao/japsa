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

/*                           Revision History                                
 * 28/05/2014 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsadev.tools.work;



import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import japsa.seq.nanopore.Fast5DetailReader;
import japsa.util.CommandLine;
import japsa.util.JapsaException;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.dev.npanalysis", 
	scriptDesc = "Analysis of an np read from a fast5 file"
	)
public class NpReadAnalysys extends CommandLine{
	//CommandLine cmdLine;
	public NpReadAnalysys(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 
				
		addString("input", null, "Name of the input bam file",true);
		addString("output", null, "Name of the output bam file");
		
		addStdHelp();
	}
	public static void main(String [] args) throws OutOfMemoryError, JapsaException, Exception{
		NpReadAnalysys cmdTool = new NpReadAnalysys ();
		args = cmdTool.stdParseLine(args);

		/**********************************************************************/
		String input = cmdTool.getStringVal("input");
		String output = cmdTool.getStringVal("output");

		addSequence(input, output);

	}
	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 * @param inFile
	 * @param pad
	 * @throws Exception 
	 * @throws JapsaException 
	 * @throws OutOfMemoryError 
	 */

	static void addSequence(String inFile, String outFile) throws OutOfMemoryError, JapsaException, Exception{
		Fast5DetailReader npReader = new Fast5DetailReader(inFile);
		npReader.readData();
		npReader.close();
		
		Fast5DetailReader.DetectionEvents events = npReader.getEvents();
		System.out.println("Events = " + events.getLength().length);
		
		DescriptiveStatistics desc = new DescriptiveStatistics();//events.getLength());
		
		System.out.println("========================== Events stats ==========================");		
		System.out.println("Event length");		
		desc.clear();
		for (int i = 0; i < events.getLength().length;i++)
			desc.addValue(events.getLength()[i] );
		
		System.out.println(" Mean = " + desc.getMean());
		System.out.println(" Std = " + desc.getStandardDeviation());
		System.out.println(" Max = " + desc.getMax());
		
		
		System.out.println("Events mean");
		desc.clear();
		for (int i = 0; i < events.getMean().length;i++)
			desc.addValue(events.getMean()[i] );
		
		System.out.println(" Mean = " + desc.getMean());
		System.out.println(" Std = " + desc.getStandardDeviation());
		System.out.println(" Max = " + desc.getMax());
		
		System.out.println("Events stdv");
		desc.clear();
		for (int i = 0; i < events.getStdv().length;i++)
			desc.addValue(events.getStdv()[i] );
		
		System.out.println(" Mean = " + desc.getMean());
		System.out.println(" Std = " + desc.getStandardDeviation());
		System.out.println(" Max = " + desc.getMax());
				
		
		System.out.println("========================== Base called stats ==========================");		
		Fast5DetailReader.BaseCallEvents bcEvents = npReader.getBcTempEvents();
		
		
		System.out.println("BC length");		
		desc.clear();
		for (int i = 0; i < bcEvents.length().length;i++)
			desc.addValue(bcEvents.length()[i] );
		
		System.out.println(" Mean = " + desc.getMean());
		System.out.println(" Std = " + desc.getStandardDeviation());
		System.out.println(" Max = " + desc.getMax());
		
		
		System.out.println("BC mean:");		
		desc.clear();
		for (int i = 0; i < bcEvents.mean().length;i++)
			desc.addValue(bcEvents.mean()[i] );
		
		System.out.println(" Mean = " + desc.getMean());
		System.out.println(" Std = " + desc.getStandardDeviation());
		System.out.println(" Max = " + desc.getMax());
		
		System.out.println("BC stdv:");		
		desc.clear();
		for (int i = 0; i < bcEvents.stdv().length;i++)
			desc.addValue(bcEvents.stdv()[i] );
		
		System.out.println(" Mean = " + desc.getMean());
		System.out.println(" Std = " + desc.getStandardDeviation());
		System.out.println(" Max = " + desc.getMax());
			
		
		System.out.println("BC move:");		
		desc.clear();
		for (int i = 0; i < bcEvents.getMove().length;i++)
			desc.addValue(bcEvents.getMove()[i] );
		
		System.out.println(" Mean = " + desc.getMean());
		System.out.println(" Std = " + desc.getStandardDeviation());
		System.out.println(" Max = " + desc.getMax());
	}

}
