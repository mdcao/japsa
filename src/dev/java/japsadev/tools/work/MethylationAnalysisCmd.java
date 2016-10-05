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


import java.io.File;
import java.io.IOException;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * Filter a bam filem based on some criteria. Input file in bam format assumed
 * to be sorted and indexed
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.dev.methyC", 
	scriptDesc = "Count fraction of methylated"
	)
public class MethylationAnalysisCmd extends CommandLine{
	//CommandLine cmdLine;
	public MethylationAnalysisCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("inputBam", null, "Name of the input bam",true);		
		addString("chrom", null, "Name of the chromosome",true);		
		addInt("qual", 0, "Minimum mapping quality");
		
		addString("output", null, "name of the output file",true);				
		//addInt("filterBits", 0, "Filter reads based on flag. Common values:\n 0    no filter\n 256  exclude secondary alignment \n 1024 exclude PCR/optical duplicates\n 2048 exclude supplementary alignments");


		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		MethylationAnalysisCmd cmdLine = new MethylationAnalysisCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String inputBam  = cmdLine.getStringVal("inputBam");
		String chrom = cmdLine.getStringVal("chrom");
		String output = cmdLine.getStringVal("output");		
		int qual = cmdLine.getIntVal("qual");		

		analyse(inputBam, chrom, qual, output);
	}

	static void analyse(String inFile, String chrom, int qual, String output) throws IOException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(inFile));						
		SAMRecordIterator samIter = samReader.query(chrom, 0, 0, false);		
		SAMSequenceRecord refSequence = samReader.getFileHeader().getSequence(chrom);

		int seqLength = refSequence.getSequenceLength();
		int [] countC = new int[seqLength], 
			countTot = new int[seqLength];

		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();
			
			if (sam.getMappingQuality() < qual)
				continue;			
			////////////////////////////////////////////////////////////////int readPos = 0;//start from 0
			
			Sequence readSeq = new Sequence(Alphabet.DNA(), sam.getReadString(), sam.getReadName());
			
			int readPos = 0;//start from 0					
			int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index

			for (final CigarElement e : sam.getCigar().getCigarElements()) {
				final int  length = e.getLength();
				switch (e.getOperator()) {
				case H :
					//nothing todo
					break; // ignore hard clips
				case P : 
					//pad is a kind of hard clipped ?? 					
					break; // ignore pads	                
				case S :
					//advance on the read
					readPos += length;
					break; // soft clip read bases	                	
				case N : 
					refPos += length; 
					break;  // reference skip

				case D ://deletion      	
					refPos += length;
					break; 	

				case I :	                	
					readPos += length;
					break;
				case M :
				case EQ:
				case X :
					for (int i = 0; i < length; i++){
						countTot[refPos + i]++;
						if (readSeq.getBase(readPos + i) == Alphabet.DNA.C)
							countC[refPos + i]++;
					}							
					readPos += length;
					refPos  += length;
					break;

				default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
				}//case
			}//for			

			////////////////////////////////////////////////////////////////

		}
		samIter.close();
		samReader.close();

		SequenceOutputStream fCount = SequenceOutputStream.makeOutputStream(output);
		fCount.print("#pos\tstart\tend\n");
		for (int i = 0; i< countTot.length;i++){
			fCount.print(i+1);
			fCount.print('\t');
			fCount.print(countC[i]);
			fCount.print('\t');
			fCount.print(countTot[i]);
			fCount.print('\n');
		}
		fCount.close();
	}
}
