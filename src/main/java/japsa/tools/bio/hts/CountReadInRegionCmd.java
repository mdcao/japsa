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

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 *A program to count reads overlapping or containted with regions from a bam/sam file 
 * 
 */

@Deployable(
		scriptName = "jsa.hts.countReads",
		scriptDesc = "Count the number of reads in some regions from a sorted, indexed bam file"
		)
public class CountReadInRegionCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(CountReadInRegionCmd.class);

	public CountReadInRegionCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());


		CommandLine.Option bamFileOpt = 
				addString("bamFile", null, "Name of the bam file",true);
		bamFileOpt.setGalaxySetting(new GalaxySetting("data", "bam",false));
		
		//addExtraGalaxyCmd("ln -s $bamFile bamFile.bam &amp;&amp; ln -s $bamFile.metadata.bam_index bamFile.bai &amp;&amp;");

		CommandLine.Option bedFileOpt = 
				addString("bedFile", null, "Name of the regions file in bed format",true);		
		bedFileOpt.setGalaxySetting(new GalaxySetting("data", "bed",false));

				
		CommandLine.Option outputOpt =
				addString("output", "-", "Name of output file, - for from standard out.");
		
		GalaxySetting galaxyOutput = new GalaxySetting("data", "tabular",true); 
		galaxyOutput.setLabel("countRead.txt");
		outputOpt.setGalaxySetting(galaxyOutput);

		addInt("flanking", 0, "Size of the flanking regions, effectively expand the region by flanking")
		.setGalaxySetting(new GalaxySetting("integer", null,false));

		addInt("qual", 0, "Minimum quality")
		.setGalaxySetting(new GalaxySetting("integer", null,false));		

		addInt("filterBits", 0, "Filter reads based on flag. Common values:\n 0    no filter\n 256  exclude secondary alignment \n 1024 exclude PCR/optical duplicates\n 2048 exclude supplementary alignments");

		addBoolean("contained", false, "Count reads contained in the region")
		.setGalaxySetting(new GalaxySetting("boolean", null,false));

		addBoolean("overlap", false, "Count number of read overlap with the region")
		.setGalaxySetting(new GalaxySetting("boolean", null,false));

		addBoolean("span", false, "Count reads span the region")
		.setGalaxySetting(new GalaxySetting("boolean", null,false));

		addStdHelp();		

		setGalaxy(annotation.scriptName());
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
		boolean overlap = cmdLine.getBooleanVal("overlap");
		boolean span = cmdLine.getBooleanVal("span");

		String bedFile = cmdLine.getStringVal("bedFile");
		if (!(contained || overlap || span)){
			System.err.print("ERROR: At least one of contained, overlap and span switched on\n");
			System.err.println(cmdLine.usageString());

			System.exit(1);
		}


		/**********************************************************************/
		ArrayList<JapsaFeature> myList = JapsaFeature.readBED(bedFile);		

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SequenceOutputStream os = SequenceOutputStream.makeOutputStream(output);

		String sep = "\t";

		int notCount = 0;
		//os.print("#CMD:" + cmdLine.fullCmd() + '\n');
		os.print("#H:chrom\tID\tstart\tend");

		/******************************************************************
		 * This is to be removed later, will get back to list of bam files
		 */

		args = new String[1];
		args[0] = cmdLine.getStringVal("bamFile");		
		/******************  End of short cut  ***********************/

		SamReader [] readers = new SamReader[args.length];
		for (int i = 0; i < readers.length; i++){
			File file = new File(args[i]);

			readers[i] = SamReaderFactory.makeDefault().open(file);
			String sampleID = null;
			//SAMReadGroupRecord groupID = 

			List<SAMReadGroupRecord> readGroups = readers[i].getFileHeader().getReadGroups();

			if (readGroups != null && readGroups.size() > 0){
				sampleID = readGroups.get(0).getSample();
			}


			if (sampleID == null){
				sampleID = file.getName().replace(".sam", "").replace(".bam","");
			}

			if (overlap)
				os.print(sep + sampleID + "_overlap");

			if (contained)
				os.print(sep + sampleID + "_contained");			

			if (span)
				os.print(sep + sampleID + "_span");
		}
		os.print("\n");

		os.flush();
		for (JapsaFeature str:myList){
			int start = str.getStart() - flanking;
			int end   = str.getEnd() +   flanking;

			if (start <= 0 ) 
				start = 1;

			//TODO: check if end > chr.length
			os.print(str.getParent() + sep + str.getID() + sep + (str.getStart() -1) + sep + str.getEnd());

			for (int i = 0; i < readers.length; i++){
				SAMRecordIterator iter = readers[i].query(str.getParent(), start, end, false);

				int countOverlap = 0, countContained = 0, countSpan = 0;

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
					countOverlap ++;

					int alignmentStart = rec.getAlignmentStart();
					int alignmentEnd = rec.getAlignmentEnd();

					if (alignmentStart >= start && alignmentEnd <= end)
						countContained ++;

					if (alignmentStart < start && alignmentEnd > end)
						countSpan ++;					

				}//while			
				iter.close();			

				if (overlap)
					os.print(sep + countOverlap);

				if (contained)
					os.print(sep + countContained);

				if (span)
					os.print(sep + countSpan);

				//os.print(sep);
				//os.print(countOverlap);
			}//for
			os.print("\n");
		}//for
		for (int i = 0; i < readers.length; i++){
			readers[i].close();
		}
		os.close();
		LOG.info("Ignore " + notCount + " reads");
	}
}
/*RST*
------------------------------------------------
*jsa.hts.countReads*: Count reads from bam files 
------------------------------------------------


<usage> 

*RST*/


