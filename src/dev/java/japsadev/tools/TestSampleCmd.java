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
 * 11/01/2012 - Minh Duc Cao: Revised 
 * 01/01/2013 - Minh Duc Cao, revised                                       
 ****************************************************************************/

package japsadev.tools;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.XAFReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
		scriptName = "jsa.dev.test",
		scriptDesc = "A test program: generate log normal distribution"
		)
public class TestSampleCmd extends CommandLine{	
	public TestSampleCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options] samFile1 [samFile2 ..]");
		setDesc(annotation.scriptDesc());

		addString("xafFile", "VNTR.xaf",  "XAF file containing repeat information");		
		addString("bamFile", null, "Name of output file, - for standard out", true);

		addStdHelp();		
	} 

	public static void main(String [] args) throws IOException, InterruptedException{		 		
		CommandLine cmdLine = new TestSampleCmd();		
		cmdLine.stdParseLine(args);
		/**********************************************************************/
		String xafFile     =  cmdLine.getStringVal("xafFile");	
		String bamFile      =  cmdLine.getStringVal("bamFile");


		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);		
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFile));
		XAFReader xafReader = new XAFReader(xafFile);

		while (xafReader.next() != null){
			//Extract 

			String chrom = xafReader.getField("chrom");			
			int endRep = Integer.parseInt(xafReader.getField("end"));			
			int rflank = 1000;

			//Get flanking information
			String field = xafReader.getField("rflank");
			if (field != null)
				rflank = Integer.parseInt(field);			

			SAMRecordIterator iter = samReader.query(chrom, endRep, endRep + rflank, true);
			while (iter.hasNext()){
				SAMRecord sam = iter.next();
				int insertSize = sam.getInferredInsertSize();

				if (insertSize > 0 && insertSize < 2000){
					System.out.println(insertSize + "\t0\t0");
				}
			}
			iter.close();

		}//for i

		xafReader.close();
	}		
}
/*RST*












 *RST*/

