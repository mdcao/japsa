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
package japsa.bio.hts;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


/**
 *A program to count reads overlapping or containted with regions from a bam/sam file
 *
 *FIXME: Generalise to any kinds of regions, not just STR  
 */

@Deployable(scriptName = "jsa.hts.countReads",
scriptDesc = "Count the number of reads in some regions from a sorted, indexed bam file")
public class CountReadInRegion {	
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/
		Deployable annotation = CountReadInRegion.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options] <s1.bam> <s2.bam> <s3.bam> ..." , annotation.scriptDesc());

		cmdLine.addString("xafFile", null, "Name of the regions file in xaf");	
		cmdLine.addString("bedFile", null, "Name of the regions file in bed\n"+
		                                   "Either xafFile or bedFile has to be specified");
				
		cmdLine.addString("output", "-", "Name of output file, - for from standard out.");

		cmdLine.addInt("flanking", 0, "Size of the flanking regions");
		cmdLine.addInt("qual", 0, "Minimum quality");		
		cmdLine.addInt("filterBits", 0, "Filter reads based on flag. Common values:\n 0    no filter\n 256  exclude secondary alignment \n 1024 exclude PCR/optical duplicates\n 2048 exclude supplementary alignments");

		cmdLine.addBoolean("contained", false, "true: Reads contained in the region; false: reads overlap with the region");

		args = cmdLine.stdParseLine_old(args);
		/**********************************************************************/		
		//Get options
		
		String output = cmdLine.getStringVal("output");
		int flanking = cmdLine.getIntVal("flanking");		
		if (flanking < 0)
			flanking = 0;		
				

		int qual = cmdLine.getIntVal("qual");
		int filter = cmdLine.getIntVal("filterBits");		
		boolean contained = cmdLine.getBooleanVal("contained");
		
		String strFile = cmdLine.getStringVal("xafFile");
		String bedFile = cmdLine.getStringVal("bedFile");
		
		if (strFile!= null &&  bedFile != null){
			System.err.println("##ERROR: only one of bedFile and strFile is specified");
			System.err.println(cmdLine.usageMessage());
			System.exit(-1);			
		}
		if (strFile== null &&  bedFile == null){
			System.err.println("##ERROR: one of bedFile and xafFile has to be specified");
			System.err.println(cmdLine.usageMessage());
			System.exit(-1);			
		}
		/**********************************************************************/
		ArrayList<JapsaFeature>  myList;
		
		if(bedFile != null)
			myList = JapsaFeature.readBED(bedFile);
		else{
			ArrayList<TandemRepeat> list = TandemRepeat.readFromFile(SequenceReader.openFile(strFile), new ArrayList<String>());			
			myList = new ArrayList<JapsaFeature>(list.size());
			
			for (TandemRepeat str:list){
				myList.add(str);
			}
		}

		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);		
		SequenceOutputStream os = SequenceOutputStream.makeOutputStream(output);
 
		char sep = '\t';
		
		int notCount = 0;
		os.print("#H:chr\tID\tstart\tend");
		
		SAMFileReader [] readers = new SAMFileReader[args.length];
		for (int i = 0; i < readers.length; i++){
			File file = new File(args[i]);		
			os.print("\t" + file.getName().replace(".sam", "").replace(".bam",""));
			readers[i] = new SAMFileReader(file);
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


