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
 * 08/03/2013 - Minh Duc Cao: Started
 * 24/06/2013 - Minh Duc Cao: updated                   
 *  
 ****************************************************************************/

package japsa.bio.hts;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;

import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


/**
 *A program to filter reads that intersect with a STR regions from a bam/sam file
 *
 *FIXME: Generalise to any kinds of regions, not just STR 
 */

@Deployable(scriptName = "jsa.hts.selectIntesect",
            scriptDesc = "Filter reads that intersect with some regions from a sorted b/sam file")

public class SelectReadIntersect {	
	public static void main(String[] args) throws Exception {	
		/*********************** Setting up script ****************************/		 
		String scriptName = "jsa.ngs.selectIntesect";
		String desc = "Filter reads that intersect with some regions from a sorted b/sam file\n";	
		//condition: 
		// 1. s/bam file sorted, the sequence in str must be the same order as in the header file of the sam file
		// 2. strs in each sequence are sorted
		// 3. no strs overlap	

		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [params]");
		/**********************************************************************/

		cmdLine.addStdInputFile();		
		cmdLine.addString("strFile", null, "Name of the regions file",true);		
		cmdLine.addString("output", "-", "Name of output sam file, - for from standard out.");
		cmdLine.addInt("flanking", 120, "Flanking regions");
		cmdLine.addInt("qual", 10, "Minimum quality");
		cmdLine.addBoolean("pair", false, "Paired-end data");
		cmdLine.addStdHelp();		

		/**********************************************************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")){
			System.out.println(desc + cmdLine.usage());			
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usage());
			System.exit(-1);
		}	
		/**********************************************************************/		

		String output = cmdLine.getStringVal("output");
		String samFile = cmdLine.getStringVal("input");
		String strFile = cmdLine.getStringVal("strFile");
		int gaps = cmdLine.getIntVal("flanking");
		qual = cmdLine.getIntVal("qual");

		if (cmdLine.getBooleanVal("pair"))			
			filterSamPair(samFile, output, strFile, gaps);
		else
			filterSam(samFile, output, strFile, gaps);
	}
	

	/**
	 * Filter out aligned reads from inFile if reads intersect with one of the 
	 * regions in strFiles
	 * 
	 * @param inFile
	 * @param outFile
	 * @param strFile
	 * @throws IOException
	 */
	static void filterSam(String inFile, String outFile, String strFile, int gaps)
			throws IOException {				
		/////////////////////////////////////////////////////////////////////////////				

		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader samReader = new  SAMFileReader(new File(inFile));
		SAMFileHeader samHeader = samReader.getFileHeader();		
		//SAMTextWriter samWriter = new SAMTextWriter(new File(outFile));
		
		boolean preOrder = false;
		SAMFileWriterFactory factory = new SAMFileWriterFactory();		
		SAMFileWriter bamWriter = factory.makeSAMOrBAMWriter(samHeader, preOrder, new File(outFile));	
		
		//bamWriter.writeHeader(samHeader.getTextHeader());
		
		ArrayList<TandemRepeat>  myList = TandemRepeat.readFromFile(SequenceReader.openFile(strFile), new ArrayList<String>());
		TandemRepeat str = myList.get(0);

		///System.err.print(str.toString()+" : ");
		int strIndex = 0;
		//int count = 0;

		int strSeqIndex = samHeader.getSequenceIndex(str.getChr());
		if (strSeqIndex < 0){
			samReader.close();		
			throw new RuntimeException("Sequence " + str.getChr() + " not found in the header of b/sam file " + inFile);
		}

		SAMRecordIterator samIter = samReader.iterator();
		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();	
			if (sam.getMappingQuality() < qual)
				continue;

			int seqIndex = sam.getReferenceIndex();

			//the samrecod is in an ealier sequence
			if (seqIndex < strSeqIndex)
				continue;

			int posStart = sam.getAlignmentStart();
			int posEnd = sam.getAlignmentEnd();

			//the sam is before the next region
			if (seqIndex == strSeqIndex && posEnd < str.getStart() - gaps)
				continue;

			//assert: seqIndex > strSeqIndex || postEnd >= str.start - gaps						
			if (seqIndex == strSeqIndex && posStart <= str.getEnd() + gaps){

				if (sam.getReadPairedFlag()){
					if (sam.getFirstOfPairFlag())
						sam.setReadName(sam.getReadName()+"-1");
					else
						sam.setReadName(sam.getReadName()+"-2");
					sam.setReadPairedFlag(false);
					sam.setProperPairFlag(false);
					sam.setFirstOfPairFlag(false);
					sam.setSecondOfPairFlag(false);
					sam.setMateAlignmentStart(0);
					sam.setMateReferenceIndex(-1);
					sam.setInferredInsertSize(0);				
				}
				bamWriter.addAlignment(sam);

				//count ++;
				continue;
			}
			//assert: seqIndex > strSeqIndex || postStart > str.End


			while (seqIndex > strSeqIndex || (seqIndex == strSeqIndex && posStart > str.getEnd() + gaps)){
				strIndex ++;
				//System.err.println(count);
			//	count = 0;
				if (strIndex < myList.size()){
					str = myList.get(strIndex);
					strSeqIndex = samHeader.getSequenceIndex(str.getChr());
					if (strSeqIndex < 0){
						samReader.close();
						throw new RuntimeException("Sequence " + str.getChr() + " not found in the header of b/sam file " + inFile);
					}
					//System.err.print(str.toString()+" : ");					
				}else{
					str = null;
					break;//while
				}
			}

			if (str == null)
				break;

			if (seqIndex == strSeqIndex && posStart <= str.getEnd() + gaps && posEnd >=  str.getStart() - gaps){
				if (sam.getReadPairedFlag()){
					if (sam.getFirstOfPairFlag())
						sam.setReadName(sam.getReadName()+"-1");
					else
						sam.setReadName(sam.getReadName()+"-2");
					sam.setReadPairedFlag(false);
					sam.setProperPairFlag(false);
					sam.setFirstOfPairFlag(false);
					sam.setSecondOfPairFlag(false);
					sam.setMateAlignmentStart(0);
					sam.setMateReferenceIndex(-1);
					sam.setInferredInsertSize(0);				
				}
				bamWriter.addAlignment(sam);
				//count ++;
				continue;
			}			
		}
		samReader.close();
		bamWriter.close();
	}


	static void filterSamPair(String inFile, String outFile, String strFile, int gaps)
			throws IOException {				
		/////////////////////////////////////////////////////////////////////////////				

		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);

		SAMFileReader samReader = new  SAMFileReader(new File(inFile));
		SAMFileHeader samHeader = samReader.getFileHeader();
		
		//FIXME: check with preOrder
		boolean preOrder = false;
		SAMFileWriterFactory factory = new SAMFileWriterFactory();		
		SAMFileWriter bamWriter = factory.makeSAMOrBAMWriter(samHeader, preOrder, new File(outFile));	

		//SAMTextWriter samWriter = new SAMTextWriter(new File(outFile));
		//samWriter.writeHeader(samHeader.getTextHeader());

		ArrayList<TandemRepeat>  myList = TandemRepeat.readFromFile(SequenceReader.openFile(strFile), new ArrayList<String>());
		TandemRepeat str = myList.get(0);

		//		System.err.print(str.toString()+" : ");
		int strIndex = 0;		

		HashSet<String> set = new HashSet<String>();

		int strSeqIndex = samHeader.getSequenceIndex(str.getChr());
		if (strSeqIndex < 0){
			samReader.close();
			throw new RuntimeException("Sequence " + str.getChr() + " not found in the header of b/sam file " + inFile);
		}

		SAMRecordIterator samIter = samReader.iterator();
		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();

			//if the pair is already in the set 
			if (set.remove(sam.getReadName())){
				bamWriter.addAlignment(sam);
				continue;
			}		

			if (sam.getMappingQuality() < qual)
				continue;

			int seqIndex = sam.getReferenceIndex();

			//the samrecod is in an ealier sequence
			if (seqIndex < strSeqIndex)
				continue;

			int posStart = sam.getAlignmentStart();
			int posEnd = sam.getAlignmentEnd();

			//the sam is before the next region
			if (seqIndex == strSeqIndex && posEnd < str.getStart() - gaps)
				continue;

			//assert: seqIndex > strSeqIndex || postEnd >= str.start - gaps						
			if (seqIndex == strSeqIndex && posStart <= str.getEnd() + gaps){
				set.add(sam.getReadName());
				bamWriter.addAlignment(sam);
				continue;
			}

			//assert: seqIndex > strSeqIndex || postStart > str.End		

			while (seqIndex > strSeqIndex || (seqIndex == strSeqIndex && posStart > str.getEnd() + gaps)){
				strIndex ++;
				//System.err.println(count);
				//count = 0;
				if (strIndex < myList.size()){
					str = myList.get(strIndex);

					int newSTRIndex = samHeader.getSequenceIndex(str.getChr());

					if (newSTRIndex != strSeqIndex){
						System.out.println(samHeader.getSequence(strSeqIndex).getSequenceName() + "   " + set.size() + " at " + System.currentTimeMillis());
						strSeqIndex = newSTRIndex; 
					}

					if (strSeqIndex < 0){
						samReader.close();
						throw new RuntimeException("Sequence " + str.getChr() + " not found in the header of b/sam file " + inFile);
					}
					//System.err.print(str.toString()+" : ");					
				}else{
					str = null;
					break;//while
				}
			}

			if (str == null)
				break;

			if (seqIndex == strSeqIndex && posStart <= str.getEnd() + gaps && posEnd >=  str.getStart() - gaps){
				bamWriter.addAlignment(sam);
				set.add(sam.getReadName());				
				continue;
			}//if			
		}//while

		bamWriter.close();
		System.out.println("Writing out " + set.size() * 2 + " at " + System.currentTimeMillis());
		//close samReader
		samReader.close();


		//Iterate the second time
		if (set.size() > 0){
			samReader = new  SAMFileReader(new File(inFile));				
			bamWriter = factory.makeSAMOrBAMWriter(samHeader, preOrder, new File("2_" + outFile));		
			
			
			//samWriter = new SAMTextWriter(new File(outFile+"_2"));		
			//samWriter.writeHeader(samHeader.getTextHeader());

			samIter = samReader.iterator();
			while (samIter.hasNext()){
				SAMRecord sam = samIter.next();
				if (set.remove(sam.getReadName())){
					bamWriter.addAlignment(sam);				
				}
			}		
			bamWriter.close();
			samReader.close();
			System.out.println("Done " + set.size() * 2 + " at " + System.currentTimeMillis());
		}
	}
	static int qual;	
}


