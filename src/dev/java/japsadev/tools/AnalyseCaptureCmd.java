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
 * 25/08/2016 - Minh Duc Cao: Created                                        
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
 * @author minhduc
 *
 */
@Deployable(
		scriptName = "jsa.dev.capFragmentAnalysis", 
		scriptDesc = "Analyse the fragment size"
		)
public class AnalyseCaptureCmd extends CommandLine{
	//CommandLine cmdLine;
	public AnalyseCaptureCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("probe", null,  "File containing probes in bam format");
		addString("logFile", "-",  "Log file");

		addString("miseq", null, "Name of read file if miseq is simulated");
		addString("pacbio", null, "Name of read file if pacbio is simulated");

		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		CommandLine cmdLine = new AnalyseCaptureCmd ();
		args = cmdLine.stdParseLine(args);	


		/**********************************************************************/

		String logFile = cmdLine.getStringVal("logFile");
		String probe       =  cmdLine.getStringVal("probe");

		String miseq       =  cmdLine.getStringVal("miseq");
		String pacbio       =  cmdLine.getStringVal("pacbio");

		if (miseq == null && pacbio == null){
			System.err.println("At least one of fragment, miseq and pacbio has to be set\n" + cmdLine.usageString());			
			System.exit(-1);
		}


		SequenceOutputStream 
		logOS =	logFile.equals("-")? (new SequenceOutputStream(System.err)):(SequenceOutputStream.makeOutputStream(logFile));


		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader probeReader = SamReaderFactory.makeDefault().open(new File(probe));		


		if (pacbio != null){
			SamReader pbReader = SamReaderFactory.makeDefault().open(new File(pacbio));		 

			SAMRecordIterator iter = pbReader.iterator();
			while (iter.hasNext()){
				SAMRecord sam = iter.next();

				int fragmentSize = sam.getReadLength();
				if (sam.getReadUnmappedFlag()){
					System.out.println(fragmentSize + "\t" + 0 + "\t" + 0);
					continue;
				}


				int fragmentStart = sam.getAlignmentStart();
				int fragmentEnd = sam.getAlignmentEnd();
				String refName = sam.getReferenceName();			
				int countProbe = countProbe(probeReader, refName, fragmentStart, fragmentEnd);
				System.out.println(fragmentSize + "\t" + countProbe + "\t" + 1);			
			}

			iter.close();
			pbReader.close();
		}
		if (miseq != null){
			SamReader msReader = SamReaderFactory.makeDefault().open(new File(miseq));		 

			SAMRecordIterator iter = msReader.iterator();
			while (iter.hasNext()){
				SAMRecord sam = iter.next();
				
				if (sam.getReadUnmappedFlag()){					
					continue;
				}

				int fragmentSize = sam.getInferredInsertSize();
				//if mate is not aligned, or mate aligned to the left, ignore
				if ( fragmentSize <= 0)
					continue;			

				int fragmentStart = sam.getAlignmentStart();
				int fragmentEnd = fragmentStart + fragmentSize;
				String refName = sam.getReferenceName();			
				int countProbe = countProbe(probeReader, refName, fragmentStart, fragmentEnd);
				System.out.println(fragmentSize + "\t" + countProbe + "\t" + 1);			
			}

			iter.close();
			msReader.close();

		}

		logOS.close();
	}

	static int countProbe(SamReader probeReader, String refName, int fragmentStart, int fragmentEnd){
		SAMRecordIterator iter = probeReader.query(refName, fragmentStart, fragmentEnd, true);

		int countProbe = 0;
		int myEnd = 0;
		while (iter.hasNext()){				
			SAMRecord sam = iter.next();

			int start = sam.getAlignmentStart();
			if (start < myEnd)
				continue;

			int end = sam.getAlignmentEnd();				
			//probe can only bind if > 90%
			if ((end - start) < 0.9 * sam.getReadLength())
				continue;

			countProbe += (end - start + 1);
			myEnd = end;
		}
		iter.close();
		return countProbe;

	}

}
