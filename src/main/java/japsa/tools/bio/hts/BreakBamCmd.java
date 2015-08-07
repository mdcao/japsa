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
 * 02/10/2013 - Minh Duc Cao: Created                                        
 * 16/11/2013 - Minh Duc Cai: Revised 
 ****************************************************************************/
package japsa.tools.bio.hts;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTextWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.io.IOException;



/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 * Program to break a fastq file to smaller pieces
 */
@Deployable(
	scriptName = "jsa.hts.breakbam",
	scriptDesc = "Break a sam/bam file to smaller ones")
public class BreakBamCmd extends CommandLine{	
	public BreakBamCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addStdInputFile();			
		addString("output", null, "Name of output s/bam file. If output file is .bam, bam format is outputed", true);
		addInt("size", 50000000, "The number of reads per file, a negative number for not spliting");	
		
		addStdHelp();		
	} 
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new BreakBamCmd();
		args = cmdLine.stdParseLine(args);
		

		String output = cmdLine.getStringVal("output");
		String inFile = cmdLine.getStringVal("input");

		int size = cmdLine.getIntVal("size");
		
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(inFile));	
		SAMFileHeader samHeader = samReader.getFileHeader();
		
		//TODO: should it be sorted?
		samHeader.setSortOrder(SortOrder.unsorted);
		
		System.out.println(samHeader.getSortOrder());
		if (size == 0)
			size = -1;
		
		int index = 1;
		
		SAMTextWriter samWriter = new SAMTextWriter(new File("P"+index+"_" + output));
		samWriter.setSortOrder(SortOrder.unsorted, false);		
		samWriter.writeHeader(samHeader.getTextHeader());
		
		
		
		int count = 0;	
		int countAll = 0;

		SAMRecordIterator samIter = samReader.iterator();
		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();			

			if (count == size){
				//write this samRecord
				samWriter.close();

				///start a new one
				index ++;
				count = 0;
				samWriter = new SAMTextWriter(new File("P"+index+"_" + output));
				samWriter.setSortOrder(SortOrder.unsorted, false);		
				samWriter.writeHeader(samHeader.getTextHeader());				
			}

			count ++;
			countAll ++;
			samWriter.addAlignment(sam);
		}//while
		samReader.close();
		System.out.println("Write " + countAll + " reads to " + index + " files" );
		samWriter.close();
	}
}
