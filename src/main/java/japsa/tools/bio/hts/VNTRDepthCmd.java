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

package japsa.tools.bio.hts;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.List;


import japsa.seq.SequenceOutputStream;
import japsa.seq.XAFReader;
import japsa.util.CommandLine;
import japsa.util.IntArray;
import japsa.util.deploy.Deployable;

/**
 * Program to compute sequencing depth of repeats and the flanking regions
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.tr.trdepth", 
scriptDesc = "Compute read depth in repeats and flanking regions")
public class VNTRDepthCmd extends CommandLine{	
	public VNTRDepthCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options] samFile1 [samFile2 ..]");
		setDesc(annotation.scriptDesc());

		addString("xafFile", "VNTR.xaf",  "XAF file containing repeat information");		
		addInt("qual", 0, "Minimum mapping quality");
		addBoolean("depth", false, "Include depth coverage (R3 and S3)");
		addInt("filterBits", 0, "Filter reads based on flag. Common values:\n 0    no filter\n 256  exclude secondary alignment \n 1024 exclude PCR/optical duplicates\n 2048 exclude supplementary alignments");
		addString("output", "-", "Name of output file, - for standard out");

		addInt("maxFragment", 0, "Fragment size, 0 if single end");
		addInt("readLength", 250, "Read length");

		addStdHelp();		
	} 

	public static void main(String [] args) throws IOException, InterruptedException{		 		
		CommandLine cmdLine = new VNTRDepthCmd();		
		String[] bamFiles = cmdLine.stdParseLine(args);
		
		if (bamFiles.length == 0){		
			System.err.println("Please enter one or more bam file \n" + cmdLine.usageString());
			System.exit(1);
		}
		

		/**********************************************************************/
		String xafFile     =  cmdLine.getStringVal("xafFile");	
		String output      =  cmdLine.getStringVal("output");
		int qual           =  cmdLine.getIntVal("qual");
		int filter = cmdLine.getIntVal("filterBits");
		boolean depth      =  cmdLine.getBooleanVal("depth");		
		int maxFragment           =  cmdLine.getIntVal("maxFragment");
		int readLength           =  cmdLine.getIntVal("readLength");

		if (bamFiles.length == 0)		
			return;		

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);		
		SamReader [] samReaders = new SamReader[bamFiles.length];		
		String [] sampleID = new String[bamFiles.length];

		for (int i = 0; i < samReaders.length;i++){
			File file = new File(bamFiles[i]);			
			samReaders[i] = SamReaderFactory.makeDefault().open(file);			
			sampleID[i] = null;
			List<SAMReadGroupRecord> readGroups = samReaders[i].getFileHeader().getReadGroups();			
			if (readGroups != null){
				for (SAMReadGroupRecord rgRec:readGroups){					
					sampleID[i] = rgRec.getSample();
					if (sampleID[i] != null)
						break;//for
				}
			}

			if (sampleID[i] == null){
				sampleID[i] = file.getName().replaceAll(".sort.bam", "").replaceAll(".bam", "");
			}


		}


		XAFReader xafReader = new XAFReader(xafFile);

		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);	
		sos.print("#H:ID\tchrom\tstart\tend\trepLen\tseqLen");
		for (int i = 0; i < bamFiles.length;i++){
			sos.print("\t" + sampleID[i] + "_R1\t" + sampleID[i] + "_S1\t" + sampleID[i] + "_R2\t" + sampleID[i] + "_S2");

			if (depth){
				sos.print("\t"+ sampleID[i] + "_R3\t" + sampleID[i] + "_S3");
			}

			sos.print("\t"+ sampleID[i] + "_B1\t"+ sampleID[i] + "_B2\t"+ sampleID[i] + "_B3\t"+ sampleID[i] + "_B4\t"+ sampleID[i] + "_B5");

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
			if (field != null)
				lflank = Integer.parseInt(field);


			int startSeq = startRep - lflank;
			int endSeq = endRep + rflank;

			//Make sure not a valid start
			if (startSeq <1)
				startSeq = 1;

			String ID = xafReader.getField("ID");
			if (ID == null){
				ID = chrom + ":" + startSeq + "-" + endRep;
			}
			sos.print(ID);
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
				int countRepInt = 0,       //R1: reads intersect with repeat
					countRepContained = 0, //R2: reads contained within repeats
					countSeqInt = 0,       //S1: reads intersect with repeat + flank
					countSeqContained = 0; //S2: reads contained within repeat + flank

				///////////////////////////////////////////////////////////////
				int leftBound = startRep - 2 * readLength;			
				int countBin1 = 0;
				int rightBound1 = startRep - readLength;
				//bin1 = (start-2*readlength,start-readlength)

				int countBin2 = 0;
				//bin2 = (start-readlength,start)
				int rightBound2 = startRep;

				int countBin3 = 0;
				//bin3 = (start, end - readLength)
				int rightBound3 = endRep - readLength;

				int countBin4 = 0;
				//bin4 = (end - readLength, end)
				int rightBound4 = endRep;

				int countBin5 = 0;
				//bin5 = (end, end+readlength)
				int rightBound5 = endRep + readLength;
				///////////////////////////////////////////////////////////////



				SAMRecordIterator iter;
				HashMap<String, IntArray> pairs = null;				

				if (maxFragment > 0){
					pairs =  new HashMap<String, IntArray>();
					iter = samReaders[i].query(chrom, startSeq - maxFragment, endSeq + maxFragment, false);					 
					while (iter.hasNext()){
						SAMRecord record = iter.next();
						String readName = record.getReadName();

						if (record.getFirstOfPairFlag())
							readName = readName + "_1";							
						else
							readName = readName + "_2";

						IntArray array = pairs.get(readName);
						if (array == null){
							array = new IntArray(16);
							pairs.put(readName,array);
						}

						array.add(record.getAlignmentStart());
					}
					iter.close();
				}

				//Reopen				
				iter = samReaders[i].query(chrom, startSeq, endSeq, false);

				while (iter.hasNext()){
					SAMRecord record = iter.next();
					if (record.getMappingQuality() < qual)
						continue;

					if ((filter & record.getFlags()) != 0){						
						continue;
					}
					
					int alignmentStart = record.getAlignmentStart();

					if (maxFragment > 0){
						boolean pass = false;
						String mateName = record.getReadName();
						if (record.getFirstOfPairFlag())
							mateName = mateName + "_2";
						else
							mateName = mateName + "_1";
						
						IntArray array = pairs.get(mateName);
						if (array !=null){
							for (int x =0; x < array.size();x++){
								int fragment = array.get(x) - alignmentStart;
								if (fragment < maxFragment && fragment > -maxFragment){
									pass = true;
									break;//for
								}								
							}							
						}//if		
						
						if (!pass)
							continue;//while
					}

					countSeqInt ++;					
					int alignmentEnd   = record.getAlignmentEnd();

					if (alignmentStart >= startSeq && alignmentEnd <= endSeq)
						countSeqContained ++;

					if (alignmentStart >= startRep && alignmentEnd <= endRep)
						countRepContained ++;

					if (alignmentStart <= endRep && alignmentEnd >= startRep)
						countRepInt ++;

					if (alignmentStart < leftBound || alignmentStart >= rightBound5)
						continue;

					if (alignmentStart < rightBound1)
						countBin1 ++;
					else if (alignmentStart < rightBound2)
						countBin2 ++;
					else if (alignmentStart < rightBound3)
						countBin3 ++;
					else if (alignmentStart < rightBound4)
						countBin4 ++;
					else if (alignmentStart < rightBound5)
						countBin5 ++;

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
				sos.print('\t');
				sos.print(countBin1);
				sos.print('\t');
				sos.print(countBin2);
				sos.print('\t');
				sos.print(countBin3);
				sos.print('\t');
				sos.print(countBin4);
				sos.print('\t');
				sos.print(countBin5);

			}//for i


			sos.print('\n');
		}//while iter
		sos.close();
		xafReader.close();
	}

}
