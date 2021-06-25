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

package japsadev.tools;

import japsa.seq.SequenceOutputStream;
import japsa.seq.XAFReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;



/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.dev.vntrDepthSum", 
scriptDesc = "Sum read depth information")
public class VNTRDepthSumCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(VNTRDepthSumCmd.class);

	public VNTRDepthSumCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options] file1 [file2...]");
		setDesc(annotation.scriptDesc());

		addString("output", "-", "Name of output file, - for standard out");
		addString("sampleID", "ID", "Sample ID");
		
		addStdHelp();		
	}  



	public static void main(String [] args) throws IOException, InterruptedException{
		CommandLine cmdLine = new VNTRDepthSumCmd();		
		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/
		String output      =  cmdLine.getStringVal("output");
		String sampleID      =  cmdLine.getStringVal("sampleID");
		
		
		if (args.length < 1){
			LOG.error("Need to supply some files",1);
			System.exit(1);
		}
		
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);
		sos.print("#H:ID\tchrom\tstart\tend\trepLen\tseqLen\t"	
		        + sampleID + "_R1\t" 	+ sampleID + "_S1\t" 	
				+ sampleID + "_R2\t" 	+ sampleID + "_S2\t"				
				+ sampleID + "_R3\t" 	+ sampleID + "_S3\t"
				+ sampleID + "_B1\t" 	+ sampleID + "_B2\t"
				+ sampleID + "_B3\t"    + sampleID + "_B4\t"
				+ sampleID + "_B5\n");
		
		XAFReader [] sReaders = new XAFReader[args.length];
		for (int i = 0; i < args.length;i++){					
			sReaders[i] = new XAFReader(args[i]);
		}
		
				
		while (sReaders[0].next() != null){
			//Extract			
			String ID = sReaders[0].getField("ID");
			String chrom = sReaders[0].getField("chrom");
			int startRep  = Integer.parseInt(sReaders[0].getField("start"));
			int endRep = Integer.parseInt(sReaders[0].getField("end"));
			int repLen = Integer.parseInt(sReaders[0].getField("repLen"));
			int seqLen = Integer.parseInt(sReaders[0].getField("seqLen"));
			
			sos.print(ID+"\t" + chrom + "\t" + startRep + "\t" + endRep + "\t"
					+ repLen + "\t" + seqLen);
			
			for (int i = 1; i < sReaders.length;i++)
				sReaders[i].next();
						
			for (int j = 6; j < 17; j ++){
				int sum = 0;
				for (int i = 0; i < sReaders.length;i++)
					sum += Integer.parseInt(sReaders[i].getField(j));
				sos.print("\t" + sum);
			}
			sos.println();			
		}
		
		sos.close();
		
		//depthAnalysis(xafFile,resample, resAllele, args, output, readLength,  model);
	}
}
