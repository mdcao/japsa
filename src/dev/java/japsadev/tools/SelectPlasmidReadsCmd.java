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
 * 05/01/2017 - Minh Duc Cao: Revised                                        
 ****************************************************************************/

package japsadev.tools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.HTSUtilities;
import japsa.util.deploy.Deployable;



/**
 * @author Son H Nguyen, MD Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.selectPlasmidReads",
	scriptDesc = "Select reads mapping to a plasmid reference"
	)
public class SelectPlasmidReadsCmd extends CommandLine{	
	public SelectPlasmidReadsCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", null, "Name of the input SAM/BAM file, - for standard input", true);
		addString("output", "-", "Name of the output file, - for standard output");
		
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new SelectPlasmidReadsCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		
		String input = cmdLine.getStringVal("input");
		String output = cmdLine.getStringVal("output");
		
		
		
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(input));
		
		

		SequenceOutputStream outFile = SequenceOutputStream.makeOutputStream(output);
		
		SAMRecordIterator iter = reader.iterator();
		String readID="";
		int[] map = null;
		ArrayList<String> readsList = new ArrayList<String>(); //store names of mapped reads
		
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();			
			
			if (rec.getReadUnmappedFlag())
				continue;
			if (rec.getMappingQuality() < 10)
				continue;

//			Arrays.fill(tmp.isMapped, myRec.refStart, myRec.refEnd, 1);
			

			//not the first occurrance				
			if (!readID.equals(rec.getReadName())){
				if(map!=null){
					int sum = IntStream.of(map).sum();
					if(((double)sum/map.length) > .9)
						readsList.add(readID);
				}
				
				readID=rec.getReadName();
				map = new int[rec.getReadLength()];
				Arrays.fill(map, 0);
			}
			
			int  [] refPositions = {rec.getStart(), rec.getEnd()}; 
			int [] pos = HTSUtilities.positionsInRead(rec, refPositions);
			if (pos[0] == 0 || pos[1] == 0)
				continue;
			Arrays.fill(map, pos[0], pos[1], 1);		

		}// while
		for(String read:readsList)
			outFile.print(read);
		
		iter.close();			
		outFile.close();


		reader.close();		
	
	}
	
}



/*RST*



 
  
  
  
  
  
  
  
  
*RST*/
  
