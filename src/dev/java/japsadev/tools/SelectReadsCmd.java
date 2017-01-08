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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.util.CommandLine;
import japsa.util.HTSUtilities;
import japsa.util.deploy.Deployable;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.selectReadsMapToPosition",
	scriptDesc = "Sample script description"
	)
public class SelectReadsCmd extends CommandLine{	
	public SelectReadsCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", null, "Name of the input file, - for standard input", true);
		addString("regions", null, "The regions to extract format chr1:s1-e1,chr2:s2-e2 no spaces", true);
		
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new SelectReadsCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		
		String input = cmdLine.getStringVal("input");
		String regions = cmdLine.getStringVal("regions");
		
		
		String [] regionArray =  regions.split(",");
		
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(input));
		
		
		for (String region:regionArray){
			String[] toks = region.split(":");
			
			
			if (toks.length < 1){
				System.err.println("region need to be in format chrX:start-end");
				System.exit(1);
			}
			String chrom = toks[0];
			toks = toks[1].split("-");
			if (toks.length < 1){
				System.err.println("region need to be in format chrX:start-end");
				System.exit(1);
			}
			int start = Integer.parseInt(toks[0]);
			int end = Integer.parseInt(toks[1]);
			
			
			System.out.println(region + ":" + (end - start) + ":"); 
			SAMRecordIterator iter = reader.query(chrom, start, end,false);
			while (iter.hasNext()){
				SAMRecord record = iter.next();
				int  [] refPositions = {start, end}; 
				int [] pos = HTSUtilities.positionsInRead(record, refPositions);
				if (pos[0] == 0 || pos[1] == 0)
					continue;
				
				String readSub = record.getReadString().substring(pos[0],pos[1]-1);
				
				System.out.printf("%5d %s %s %s\n", readSub.length(),readSub.substring(0, 22),readSub.substring(readSub.length() - 24), readSub);			
				
			}
			iter.close();			
		}

		reader.close();		
	
	}
	
}



/*RST*



 
  
  
  
  
  
  
  
  
*RST*/
  
