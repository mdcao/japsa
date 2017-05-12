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
 * 10/01/2017 - Minh Duc Cao: Created
 ****************************************************************************/
package japsa.tools.bio.hts;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTextWriter;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;



/**
 * Tools such as bwa does not store read sequence in the secondary alignments.
 * This tool correct this. It is useful for select a region later on
 * 
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 */
@Deployable(
		scriptName = "jsa.hts.fixsam",
		scriptDesc = "Add read sequences to secondary alignment, applied only for"
				+ "\nsam files by bwa without sorting."
				+ "\nNote it does not support paired-end at this version")
public class AddReadSequence2SamCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(AddReadSequence2SamCmd.class);
	public AddReadSequence2SamCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());		

		addString("input", null, "Name of the input file, - for standard input", true);
		addString("output", null, "Name of output s/bam file. If output file is .bam, bam format is outputed", true);

		addStdHelp();		
	} 
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new AddReadSequence2SamCmd();
		args = cmdLine.stdParseLine(args);		

		String output = cmdLine.getStringVal("output");
		String inFile = cmdLine.getStringVal("input");

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader 
		samReader = "-".equals(inFile)?SamReaderFactory.makeDefault().open(SamInputResource.of(System.in)):
			SamReaderFactory.makeDefault().open(new File(inFile));


		SAMFileHeader samHeader = samReader.getFileHeader();

		SAMTextWriter samWriter = null;
		if ("-".equals(output)){
			samWriter = new SAMTextWriter(System.out);
		}else{
			samWriter = new SAMTextWriter(new File(output));			
		}			

		samWriter.setSortOrder(SortOrder.unsorted, false);		
		samWriter.writeHeader(samHeader.getTextHeader());

		String readID = "";
		String readSequence  = null;
		String revSequence = null;
		boolean firstFlag = true;

		SAMRecordIterator samIter = samReader.iterator();
		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();
			if (!readID.equals(sam.getReadName())){				
				readSequence = sam.getReadString();
				if (readSequence.length() < 2){
					LOG.error("Some thing wrong " + sam.getReadName());
					continue;
				}
				readID = sam.getReadName();
				revSequence = null;
				firstFlag = sam.getReadNegativeStrandFlag();
				readID = sam.getReadName();
			}else if (sam.getReadString().length() < 2){
				if (sam.getReadNegativeStrandFlag() == firstFlag)
					sam.setReadString(readSequence);
				else{
					if (revSequence == null){
						Sequence seq = new Sequence(Alphabet.DNA6(), readSequence, "somename");
						revSequence = Alphabet.DNA16.complement(seq).toString();
					}
					sam.setReadString(revSequence);
				}
			}			
			samWriter.writeAlignment(sam);
		}//while
		samReader.close();
		samWriter.close();
	}
}
