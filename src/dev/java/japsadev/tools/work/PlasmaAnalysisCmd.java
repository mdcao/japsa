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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
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
	scriptName = "jsa.dev.plasma", 
	scriptDesc = "Analysis of plasma sequencing"
	)
public class PlasmaAnalysisCmd extends CommandLine{
	//CommandLine cmdLine;
	public PlasmaAnalysisCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("inputBam", null, "Name of the input bam",true);		
		addString("chrom", null, "Name of the chromosome",true);
		addInt("qual", 0, "Minimum mapping quality");
		addInt("min", 0, "Minimum fragment length");
		addInt("max", 10000, "Maximum fragment length");
		addInt("window", 10, "Half window sise");


		addString("output", null, "Name of the output file",true);

		//addInt("filterBits", 0, "Filter reads based on flag. Common values:\n 0    no filter\n 256  exclude secondary alignment \n 1024 exclude PCR/optical duplicates\n 2048 exclude supplementary alignments");


		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		PlasmaAnalysisCmd cmdLine = new PlasmaAnalysisCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String inputBam  = cmdLine.getStringVal("inputBam");
		//String controlBam = cmdLine.getStringVal("controlBam");
		String chrom = cmdLine.getStringVal("chrom");
		String output = cmdLine.getStringVal("output");
		int qual = cmdLine.getIntVal("qual");
		int min = cmdLine.getIntVal("min");
		int max = cmdLine.getIntVal("max");
		int window = cmdLine.getIntVal("window");

		analyse(inputBam, chrom, qual, min, max, window, output);

	}

	static void analyse(String inFile, String chrom, int qual, int min, int max, int window, String output) throws IOException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(inFile));

		SAMRecordIterator samIter = samReader.query(chrom, 0, 0, false);

		SAMSequenceRecord refSequence = samReader.getFileHeader().getSequence(chrom);

		int length = refSequence.getSequenceLength();
		int [] countStart = new int[length], 
			countEnd = new int[length];

		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();
			if (sam.getMappingQuality() < qual)
				continue;			

			int insertSize = sam.getInferredInsertSize();

			if (insertSize <= 0)//this step may be redundant 
				continue;

			if (insertSize < min)
				continue;

			if (insertSize > max)
				continue;		


			int posStart = sam.getAlignmentStart() - 1;//-1 for 0-index
			int posEnd = posStart + insertSize;		

			countStart[posStart] ++;

			if (posEnd < length)
				countEnd[posEnd] ++;	
		}
		samIter.close();
		samReader.close();

		
		writeBedGraph(chrom, countStart, window, output + "start.bedgraph");
		writeBedGraph(chrom, countEnd, window, output + "end.bedgraph");
		for (int i = 0; i< countStart.length;i++){
			countStart[i] += countEnd[i];			
		}
		writeBedGraph(chrom, countStart, window, output + "tot.bedgraph");		
	}

	static void writeBedGraph(String chrom, int [] countStart, int halfWindowSize, String fileName) throws IOException{
		SequenceOutputStream fCount = SequenceOutputStream.makeOutputStream(fileName);
		fCount.print("track type=bedGraph\n");		
		char sep = '\t';		
		int sum = 0;		
		for (int i = 0; i < Math.min(halfWindowSize, countStart.length);i++)
			sum += countStart[i];

		//fCount.print("#pos\tstart\tend\n");
		for (int i = 0; i< countStart.length;i++){
			int newPos = i + halfWindowSize; 
			if (newPos < countStart.length)
				sum += countStart[newPos];			

			int oldPos = i - halfWindowSize - 1;
			if (oldPos >= 0)
				sum -= countStart[oldPos];


			fCount.print(chrom);
			fCount.print(sep);
			fCount.print(i);
			fCount.print(sep);						
			fCount.print(i+1);
			fCount.print(sep);			
			fCount.print(sum);
			fCount.print('\n');	
		}		
		fCount.close();
	}
}
