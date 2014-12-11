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
package japsa.bio.hts;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;


/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 * Program to break a fastq file to smaller pieces
 */
@Deployable(scriptName = "jsa.hts.breakbam",
            scriptDesc = "Break a sam/bam file to smaller ones")
public class BreakBam{
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{

		/*********************** Setting up script ****************************/		 
		String scriptName = "jsa.ngs.breakbam";
		String desc = "Break a sam/bam file to smaller ones\n";
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [options]");
		/**********************************************************************/

		cmdLine.addStdInputFile();			
		cmdLine.addString("output", null, "Name of output s/bam file. If output file is .bam, bam format is outputed", true);
		cmdLine.addInt("size", 50000000, "The number of reads per file, a negative number for not spliting");

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
		String inFile = cmdLine.getStringVal("input");
		
		int size = cmdLine.getIntVal("size");
		
		//example
		//@ERR091787.1 HSQ955_155:1:1101:1266:2037/1
		//TGCAGNGGTAAATTGACCCAAGAAACTTATTTAAGACTATCAGCT
		//+
		//@BCFF#2B>CFHHIJIJJJJJIJJIJJJGHJIIIJIJJIJJIJJJ
		
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader samReader = new  SAMFileReader(new File(inFile));
		SAMFileHeader samHeader = samReader.getFileHeader();
		
		samHeader.setSortOrder(SortOrder.unsorted);
		System.out.println(samHeader.getSortOrder());
		if (size == 0)
			size = -1;
		boolean preOrder = false;
			
		int index = 1;
		SAMFileWriterFactory factory = new SAMFileWriterFactory(); 
		SAMFileWriter bamWriter = factory.makeSAMOrBAMWriter(samHeader, preOrder, new File("P"+index+"_" + output));
		
		
		//SequenceOutputStream outStream = 
		//		new  SequenceOutputStream("P"+index+"_" + output);
		int count = 0;	
		int countAll = 0;
		
		SAMRecordIterator samIter = samReader.iterator();
		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();			
			//if (sam.getFlags() >= 256)
			//	continue;
			
			if (count == size){
				//write this samRecord
				bamWriter.close();

				///start a new one
				index ++;
				count = 0;
				bamWriter = factory.makeSAMOrBAMWriter(samHeader, preOrder, new File("P"+index+"_" + output));
			}
			
			count ++;
			countAll ++;
			bamWriter.addAlignment(sam);
		}//while
		samReader.close();
		System.out.println("Write " + countAll + " reads to " + index + " files" );
		bamWriter.close();
	}
}
