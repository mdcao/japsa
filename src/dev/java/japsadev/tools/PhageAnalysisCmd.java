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
 * 05/01/2017 - Minh Duc Cao: Revised                                        
 ****************************************************************************/

package japsadev.tools;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.HTSUtilities;
import japsa.util.deploy.Deployable;
import japsadev.bio.np.phage.SequenceExtractor;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.phageAnalysis",
	scriptDesc = "Analysis phage dataset. Now only extract insert sequences..."
	)
public class PhageAnalysisCmd extends CommandLine{	
	public PhageAnalysisCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", "-", "Name of the input file, - for standard input", true);
		addString("format", "fasta", "Format of the input: SAM/BAM or FASTA/FASTQ");
		addString("plasmid", null, "Name of a sample plasmid file in FASTA format", true);
		addString("output", "out.fasta", "Name of the output file, - for standard input");		
		
		addString("bwaExe", "bwa", "Path to BWA mem.");	
				
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new PhageAnalysisCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		
		String 	input = cmdLine.getStringVal("input"),
				format = cmdLine.getStringVal("format"),
				plasmid = cmdLine.getStringVal("plasmid"),
				bwaExe = cmdLine.getStringVal("bwaExe"),
				output = cmdLine.getStringVal("output");
		try {
			SequenceExtractor vector = new SequenceExtractor(plasmid, bwaExe, 1658, 2735);
			
			vector.extractInsertSequence(input, 0, format, 4, output);
			
//			ArrayList<Sequence> list = SequenceReader.readAll("/home/s.hoangnguyen/Projects/Phage/test.fasta",Alphabet.DNA());
//			for(Sequence e:list)
//				vector.extractInsertSequence(e);
//				
				
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}
		
}



/*RST*



 
  
  
  
  
  
  
  
  
*RST*/
  
