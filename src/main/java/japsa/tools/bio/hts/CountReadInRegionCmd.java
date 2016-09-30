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

/****************************************************************************
 *                           Revision History                                
 * 15/05/2014 - Minh Duc Cao: Started
 *  
 ****************************************************************************/
package japsa.tools.bio.hts;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 *A program to count reads overlapping or containted with regions from a bam/sam file 
 * 
 */

@Deployable(
		scriptName = "jsa.hts.countReads",
		scriptDesc = "Count the number of reads in some regions from a sorted, indexed bam file"
		)
public class CountReadInRegionCmd extends CommandLine{	
	public CountReadInRegionCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options] <s1.bam> <s2.bam> <s3.bam> ...");
		setDesc(annotation.scriptDesc());	
		addString("bedFile", null, "Name of the regions file in bed format",true);
		addString("output", "-", "Name of output file, - for from standard out.");
		addInt("flanking", 0, "Size of the flanking regions");
		addInt("qual", 0, "Minimum quality");		
		addInt("filterBits", 0, "Filter reads based on flag. Common values:\n 0    no filter\n 256  exclude secondary alignment \n 1024 exclude PCR/optical duplicates\n 2048 exclude supplementary alignments");
		addBoolean("contained", false, "true: Reads contained in the region; false: reads overlap with the region");

		addStdHelp();		
	}  

	public static void main(String[] args) throws IOException {

		CommandLine cmdLine = new CountReadInRegionCmd();
		args = cmdLine.stdParseLine(args);


		String output = cmdLine.getStringVal("output");
		int flanking = cmdLine.getIntVal("flanking");		
		if (flanking < 0)
			flanking = 0;		


		int qual = cmdLine.getIntVal("qual");
		int filter = cmdLine.getIntVal("filterBits");		
		boolean contained = cmdLine.getBooleanVal("contained");		
		String bedFile = cmdLine.getStringVal("bedFile");


		/**********************************************************************/
		ArrayList<JapsaFeature> myList = JapsaFeature.readBED(bedFile);		

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SequenceOutputStream os = SequenceOutputStream.makeOutputStream(output);

		char sep = '\t';

		int notCount = 0;
		os.print("#H:chr\tID\tstart\tend");

		SamReader [] readers = new SamReader[args.length];
		for (int i = 0; i < readers.length; i++){
			File file = new File(args[i]);		
			os.print("\t" + file.getName().replace(".sam", "").replace(".bam",""));
			readers[i] = SamReaderFactory.makeDefault().open(file);				
		}
		os.print("\n");


		for (JapsaFeature str:myList){
			int start = str.getStart() - flanking;
			int end   = str.getEnd() +   flanking;

			if (start < 0 ) 
				start = 0;
			//TODO: check if end > chr.length
			os.print(str.getParent() + sep + str.getID() + sep + str.getStart() + sep + str.getEnd());

			for (int i = 0; i < readers.length; i++){
				SAMRecordIterator iter = readers[i].query(str.getParent(), start, end, contained);
				int count = 0;
				while (iter.hasNext()){
					SAMRecord rec = iter.next();

					//Check qualilty
					if(rec.getMappingQuality() < qual){
						notCount ++;
						continue;
					}				
					if ((filter & rec.getFlags()) != 0){
						notCount ++;
						continue;
					}				
					count ++;
				}//while			
				iter.close();			
				os.print(sep);
				os.print(count);
			}//for
			os.print("\n");
		}//for
		for (int i = 0; i < readers.length; i++){
			readers[i].close();
		}
		os.close();
		Logging.info("Ignore " + notCount + " reads");
	}

}



