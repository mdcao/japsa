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
 * 09/09/2016 - Minh Duc Cao started 
 *                                        
 ****************************************************************************/

package japsadev.tools;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
		scriptName = "jsa.dev.capDesign",
		scriptDesc = "Capture probe design"
		)
public class CaptureProbeDesignCmd extends CommandLine{	
	public CaptureProbeDesignCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());


		addString("bamFile",null,"Bam file",true);
		addString("output","-","Output file");
		addInt("distance",2000, "Acceptable distance");
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new CaptureProbeDesignCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String bamFile = cmdLine.getStringVal("bamFile");
		String output = cmdLine.getStringVal("output");
		int distance = cmdLine.getIntVal("distance");

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));

		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);



		SAMRecordIterator iter = reader.iterator();
		String currentName = "";
		String currentSeq = "";
		String currentQual = "";

		String myChrom = "";
		int myStart = 0;

		//int distance = 1000;

		int counterGood = 0;
		int counterBad = 0;

		while (iter.hasNext()){
			SAMRecord record = iter.next();
			String readName = record.getReadName();			

			if (record.getReadUnmappedFlag()){
				sos.print(readName + "\t0\t0\t" + record.getReadString() + "\t" + record.getBaseQualityString() + "\n");
				//odd
				continue;				
			}			

			if (!readName.equals(currentName)){
				if (currentName.length() > 0)
					sos.print(currentName + "\t" + counterGood + "\t"+ counterBad + "\t" + currentSeq + "\t" + currentQual + "\n");

				currentName = readName;
				currentSeq  = record.getReadString();
				currentQual  = record.getBaseQualityString();
				counterGood = 0;
				counterBad = 0;
				String [] toks = readName.split("_");
				myChrom = toks[0];
				myStart = Integer.parseInt(toks[1]);
			}

			String chrom = record.getReferenceName();
			if (!chrom.equals(myChrom)){
				counterBad ++;
				continue;
			}

			if (Math.abs(myStart - record.getAlignmentStart()) > distance){
				counterBad ++;				
			}else
				counterGood ++;
		}
		sos.print(currentName + "\t" + counterGood + "\t"+ counterBad + "\t" + currentSeq + "\t" + currentQual + "\n");

		iter.close();
		sos.close();

	}	
}
/*RST*



asdfds
sfasf 







 *RST*/

