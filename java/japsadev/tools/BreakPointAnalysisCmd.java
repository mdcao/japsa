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
 * 02/09/2016 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsadev.tools;


import java.io.File;
import java.io.IOException;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Analysis of break-point patterns
 * 
 * 
 * @author minhduc
 *
 */
@Deployable(
		scriptName = "jsa.dev.breakpoint", 
		scriptDesc = "Analysis of break point pattern"
		)
public class BreakPointAnalysisCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(BreakPointAnalysisCmd.class);

	public BreakPointAnalysisCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("bamFile", null, "Bam file",true);		
		addString("output", null, "output",true);

		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		BreakPointAnalysisCmd cmdLine = new BreakPointAnalysisCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String bamFile  = cmdLine.getStringVal("bamFile");		
		String output = cmdLine.getStringVal("output");



		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFile));


		List<SAMSequenceRecord> seqList =  
				samReader.getFileHeader().getSequenceDictionary().getSequences();


		int halfWindow = 5;
		int norm = 100;
		for (SAMSequenceRecord sequenceRecord:seqList){
			String chrom = sequenceRecord.getSequenceName();			
			int [] countThrough  = new int[sequenceRecord.getSequenceLength()];
			int [] countBreaks = new int[sequenceRecord.getSequenceLength()];			
			//Arrays.fill(countReads, 0);
			//Arrays.fill(countBreaks, 0);

			SAMRecordIterator samIter = samReader.query(chrom, 0, 0, false);
			while (samIter.hasNext()){
				SAMRecord samRecord = samIter.next();
				Cigar cigar = samRecord.getCigar();

				int clipLeft = 0;
				int clipRight = 0;

				//check the left clip
				CigarElement endElement = cigar.getCigarElement(0);
				if (endElement.getOperator() == CigarOperator.S 
						|| endElement.getOperator() == CigarOperator.H
						|| endElement.getOperator() == CigarOperator.P){
					clipLeft = endElement.getLength();

				}

				//Now look at the right clip
				endElement = cigar.getCigarElement(cigar.numCigarElements() - 1);				
				if (endElement.getOperator() == CigarOperator.S 
						|| endElement.getOperator() == CigarOperator.H
						|| endElement.getOperator() == CigarOperator.P){
					clipRight = endElement.getLength();				
				}

				int start = samRecord.getAlignmentStart();
				int end = samRecord.getAlignmentEnd();

				if (clipLeft > norm){
					int i = start - 1 - halfWindow;
					if (i < 0) 
						i=0;
					for (; i < start - 1 + halfWindow && i < countBreaks.length;i++){
						countBreaks[i] ++;	
					}					
				}

				if (clipRight > norm){
					int i = end - 1 - halfWindow;
					if (i < 0) 
						i=0;
					for (; i < end - 1 + halfWindow && i < countBreaks.length;i++){
						countBreaks[i] ++;	
					}
				}				
				for (int i = start  + halfWindow; i < end - halfWindow - 2;i++){
					countThrough[i]++;
				}				
			}//while
			samIter.close();

			LOG.info("Write");
			SequenceOutputStream tCount = SequenceOutputStream.makeOutputStream(output + "_" + chrom +"_through.bedgraph");
			SequenceOutputStream bCount = SequenceOutputStream.makeOutputStream(output + "_" + chrom +"_breaks.bedgraph");
			tCount.print("track type=bedGraph\n");		
			bCount.print("track type=bedGraph\n");
			char sep = '\t';		
			for (int i = 0; i < countThrough.length;i++){
				tCount.print(chrom);
				tCount.print(sep);
				tCount.print(i);
				tCount.print(sep);						
				tCount.print(i+1);
				tCount.print(sep);			
				tCount.print(countThrough[i]);
				tCount.print('\n');

				bCount.print(chrom);
				bCount.print(sep);
				bCount.print(i);
				bCount.print(sep);						
				bCount.print(i+1);
				bCount.print(sep);			
				bCount.print(countBreaks[i]);
				bCount.print('\n');

			}
			tCount.close();
			bCount.close();			
		}//while

	}
}
