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


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.phageAnalysis",
	scriptDesc = "Analysis phage dataset"
	)
public class PhageAnalysisCmd extends CommandLine{	
	public PhageAnalysisCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", null, "Name of the input file, - for standard input", true);
		addString("output", "out.fasta", "Name of the output file, - for standard input");		
				
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new PhageAnalysisCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		
		String input = cmdLine.getStringVal("input");
		String output = cmdLine.getStringVal("output");

		
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(input));
		
		
		int startFront = 1658;
		int endFront = 1758;
		
		int startBack = 2635;
		int endBack = 2735;
					
		SAMRecordIterator iter = reader.iterator();
		int count = 0;
		String currentName = "";
		boolean direction = true;
		int readFront = 0, readBack = 0;
		SequenceOutputStream outFile = SequenceOutputStream.makeOutputStream(output);
		
		while (iter.hasNext()){
			count ++;
			SAMRecord record = iter.next();
			if (record.getReadString().length() < 10){
				System.out.println("== " + record.getReadName());
				continue;//while
			}
			
			if (!record.getReadName().equals(currentName)){
				currentName = record.getReadName();
				direction = true;
				readFront = 0;
				readBack = 0;
			}
			//FIXME: this can be improved
			if (record.getAlignmentStart() <= startFront && record.getAlignmentEnd() >= endFront){
				int  [] refPositions = {startFront, endFront}; 
				int [] pos = HTSUtilities.positionsInRead(record, refPositions);
				if (pos[0] > 0 || pos[1] > 0){
					//System.out.printf("%5d %s %s %d %d %b\n", count, record.getReadName(), "FRONT", pos[0], pos[1], record.getReadNegativeStrandFlag());
					readFront = pos[0];
					if (readBack > 0 && direction == record.getReadNegativeStrandFlag()){
						if (readBack < readFront){
							System.err.printf("Bugger 1\n");							
						}else{
							String readSub = record.getReadString().substring(readFront,readBack);
							Sequence rs = new Sequence(Alphabet.DNA16(), readSub, record.getReadName());
							rs.writeFasta(outFile);													
						}
						direction = record.getReadNegativeStrandFlag();
					}
				}
					
			}
			
			if (record.getAlignmentStart() <= startBack && record.getAlignmentEnd() >= endBack){
				int  [] refPositions = {startBack, endBack}; 
				int [] pos = HTSUtilities.positionsInRead(record, refPositions);
				if (pos[0] > 0 || pos[1] > 0){
					//System.out.printf("%5d %s %s %d %d %b\n", count, record.getReadName(), "FRONT", pos[0], pos[1], record.getReadNegativeStrandFlag());
					readBack = pos[0];
					if (readFront > 0 && direction == record.getReadNegativeStrandFlag()){
						if (readBack < readFront){
							System.err.printf("Bugger 2\n");							
						}else{
							String readSub = record.getReadString().substring(readFront,readBack);
							Sequence rs = new Sequence(Alphabet.DNA16(), readSub, record.getReadName());
							rs.writeFasta(outFile);													
						}
						direction = record.getReadNegativeStrandFlag();
					}
				}				
			}		
		}
		iter.close();
		outFile.close();
		reader.close();
	}	
}



/*RST*



 
  
  
  
  
  
  
  
  
*RST*/
  
