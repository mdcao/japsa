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

package japsa.bio.tr;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import japsa.seq.SequenceOutputStream;
import japsa.seq.XAFReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * Program to compute sequencing depth of repeats and the flanking regions
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.tr.trdepth", 
scriptDesc = "Compute read depth in repeats and flanking regions")
public class VNTRDepth {
	public static void main(String [] args) throws IOException, InterruptedException{
		/*********************** Setting up script ****************************/
		Deployable annotation = VNTRDepth.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options] samFile1 [samFile2 ..]", annotation.scriptDesc());		
		/**********************************************************************/

		cmdLine.addString("xafFile", "VNTR.xaf",  "XAF file containing repeat information");		
		cmdLine.addInt("qual", 0, "Minimum mapping quality");
		cmdLine.addBoolean("depth", false, "Include depth coverage (R3 and S3)");
		cmdLine.addInt("filterBits", 0, "Filter reads based on flag. Common values:\n 0    no filter\n 256  exclude secondary alignment \n 1024 exclude PCR/optical duplicates\n 2048 exclude supplementary alignments");
		cmdLine.addString("output", "-", "Name of output file, - for standard out");


		String[] bamFiles = cmdLine.stdParseLine(args);			
		/**********************************************************************/
		String xafFile     =  cmdLine.getStringVal("xafFile");	
		String output      =  cmdLine.getStringVal("output");
		int qual           =  cmdLine.getIntVal("qual");
		int filter = cmdLine.getIntVal("filterBits");
		boolean depth      =  cmdLine.getBooleanVal("depth");		

		if (bamFiles.length == 0)		
			return;		

		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);

		SAMFileReader [] samReaders = new SAMFileReader[bamFiles.length];
		String [] sampleID = new String[bamFiles.length];

		for (int i = 0; i < samReaders.length;i++){
			File file = new File(bamFiles[i]);
			sampleID[i] = file.getName();
			sampleID[i] = sampleID[i].replaceAll(".sort.bam", "");
			sampleID[i] = sampleID[i].replaceAll(".bam", "");
			samReaders[i] = new  SAMFileReader(file);
		}


		XAFReader xafReader = new XAFReader(xafFile);

		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);	
		sos.print("#H:ID\tchrom\tstart\tend\trepLen\tseqLen");
		for (int i = 0; i < bamFiles.length;i++){
			sos.print("\t" + sampleID[i] + "_R1\t" + sampleID[i] + "_S1\t" + sampleID[i] + "_R2\t" + sampleID[i] + "_S2");

			if (depth){
				sos.print("\t"+ sampleID[i] + "_R3\t" + sampleID[i] + "_S3");
			}
		}		
		sos.print('\n');

		while (xafReader.next() != null){
			//Extract 

			String chrom = xafReader.getField("chrom");
			int startRep  = Integer.parseInt(xafReader.getField("start"));
			int endRep = Integer.parseInt(xafReader.getField("end"));

			int rflank = 1000, lflank = 1000;

			//Get flanking information
			String field = xafReader.getField("flank");
			if (field != null)
				lflank = rflank = Integer.parseInt(field);

			field = xafReader.getField("rflank");
			if (field != null)
				rflank = Integer.parseInt(field);

			field = xafReader.getField("lflank");
			lflank = Integer.parseInt(field);


			int startSeq = startRep - lflank;
			int endSeq = endRep + rflank;

			//Make sure not a valid start
			if (startSeq <1)
				startSeq = 1;

			sos.print(xafReader.getField("ID"));
			sos.print('\t');
			sos.print(chrom);
			sos.print('\t');
			sos.print(startRep);
			sos.print('\t');
			sos.print(endRep);
			sos.print('\t');						
			sos.print(endRep - startRep + 1);
			sos.print('\t');
			sos.print(endSeq - startSeq + 1);			

			for (int i = 0; i < samReaders.length;i++){	
				int     countRepInt = 0,       //R1: reads intersect with repeat
						countRepContained = 0, //R2: reads conained within repeats
						countSeqInt = 0,       //S1: reads intersect with repeat + flank
						countSeqContained = 0; //S2: reads contained within repeat + flank


				SAMRecordIterator iter = samReaders[i].query(chrom, startSeq, endSeq, false);				
				while (iter.hasNext()){
					SAMRecord record = iter.next();

					if (record.getMappingQuality() < qual)
						continue;
					
					if ((filter & record.getFlags()) != 0){						
						continue;
					}	

					countSeqInt ++;

					int alignmentStart = record.getAlignmentStart();
					int alignmentEnd   = record.getAlignmentEnd();

					if (alignmentStart >= startSeq && alignmentEnd <= endSeq)
						countSeqContained ++;

					if (alignmentStart >= startRep && alignmentEnd <= endRep)
						countRepContained ++;

					if (alignmentStart <= endRep && alignmentEnd >= startRep)
						countRepInt ++;
				}
				iter.close();

				sos.print('\t');
				sos.print(countRepInt);
				sos.print('\t');
				sos.print(countSeqInt);

				sos.print('\t');
				sos.print(countRepContained);
				sos.print('\t');
				sos.print(countSeqContained);

				if (depth){
					String cmd = "samtools depth -q "+qual+" -r " + chrom + ":" + startSeq + "-" + endRep +"  " + bamFiles[i];

					//Logging.info("Run " + cmd);
					Process process = Runtime.getRuntime().exec(cmd);
					BufferedReader depthReader = new BufferedReader (new InputStreamReader(process.getInputStream()));
					String depthLine = "";
					long
					countDepthSeq = 0,
					countDepthRep = 0;

					while ((depthLine = depthReader.readLine()) != null ){
						String [] toks = depthLine.trim().split("\t");
						long depthCount = Long.parseLong(toks[2]);					
						countDepthSeq += depthCount;

						int pos = Integer.parseInt(toks[1]);
						if (pos >= startRep && pos <= endRep)
							countDepthRep += depthCount;					
					}

					depthReader.close();				
					process.waitFor();

					sos.print('\t');
					sos.print(countDepthRep);
					sos.print('\t');
					sos.print(countDepthSeq);
				}
			}//for i


			sos.print('\n');
		}//while iter
		sos.close();
		xafReader.close();
	}

}
