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

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;

import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

/**
 * Select reads that span an STR from a bam/sam file
 * FIXME: 1. Generalise to any features, not just STR
 * FIXME: 2. Test outputing to sam/bam file
 */

@Deployable(scriptName = "jsa.hts.selectSpan",
scriptDesc = "Filter reads that span some regions from a sorted b/sam file")
public class SelectReadSpan {	
	static int pad = 3;

	public static void main(String[] args) throws Exception {		
		/*********************** Setting up script ****************************/		 
		String scriptName = "jsa.ngs.selectSpan";
		String desc = "Filter reads that span some STR regions from a sorted b/sam file\n";	
//condition: 
// 1. s/bam file sorted, the sequence in tr must be the same order as in the header file of the sam file
// 2. trs in each sequence are sorted
// 3. no trs overlap		
		
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [options]");
		/**********************************************************************/		
		cmdLine.addStdInputFile();		
		cmdLine.addString("trFile", null, "Name of the tr file",true);		
		cmdLine.addString("output", "-", "Name of output sam file, - for from standard out.");
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
		String trFile = cmdLine.getStringVal("trFile");
			
		filterSam(samFile, output, trFile);		
	}
	
	static void filterSam(String inFile, String outFile, String trFile)
			throws IOException {				
		/////////////////////////////////////////////////////////////////////////////				
		
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader samReader = new  SAMFileReader(new File(inFile));
		SAMFileHeader samHeader = samReader.getFileHeader();		
				
		boolean preOrder = false;
		SAMFileWriterFactory factory = new SAMFileWriterFactory();		
		SAMFileWriter bamWriter = factory.makeSAMOrBAMWriter(samHeader, preOrder, new File(outFile));	
			
		
		ArrayList<TandemRepeat>  myList = TandemRepeat.readFromFile(SequenceReader.openFile(trFile), new ArrayList<String>());
		TandemRepeat tr = myList.get(0);
		
		System.err.print(tr.toString()+" : ");
		int trIndex = 0;
		int count = 0;
		
		int trSeqIndex = samHeader.getSequenceIndex(tr.getChr());
		if (trSeqIndex < 0){
			samReader.close();
			throw new RuntimeException("Sequence " + tr.getChr() + " not found in the header of b/sam file " + inFile);
		}
		
		SAMRecordIterator samIter = samReader.iterator();
		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();			
			

			int seqIndex = sam.getReferenceIndex();
			
			//the samrecod is in an ealier sequence
			if (seqIndex < trSeqIndex)
				continue;
			
			int posStart = sam.getAlignmentStart();
			int posEnd = sam.getAlignmentEnd();
			
			if (seqIndex == trSeqIndex && posEnd <= tr.getEnd() + pad)
				continue;
			
			if (seqIndex == trSeqIndex && posStart < tr.getStart() - pad){
				bamWriter.addAlignment(sam);
				count ++;
				continue;
			}
			
			while (seqIndex > trSeqIndex || (seqIndex == trSeqIndex && posStart > tr.getStart() - pad)){
				trIndex ++;
				System.err.println(count);
				count = 0;
				if (trIndex < myList.size()){
					tr = myList.get(trIndex);
					trSeqIndex = samHeader.getSequenceIndex(tr.getChr());
					if (trSeqIndex < 0){
						samReader.close();
						throw new RuntimeException("Sequence " + tr.getChr() + " not found in the header of b/sam file " + inFile);
					}
					System.err.print(tr.toString()+" : ");					
				}else{
					tr = null;
					break;//while
				}
			}
			
			if (tr == null)
				break;
			
			if (seqIndex == trSeqIndex && posStart > tr.getStart() - pad && posEnd < tr.getEnd() + pad){
				bamWriter.addAlignment(sam);
				count ++;
				continue;
			}
		}
		bamWriter.close();
		samReader.close();
	}
}


