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
 * 11/05/2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.tools.bio.np;

import java.io.IOException;

import japsa.seq.Alphabet;
import japsa.seq.FastqReader;
import japsa.seq.FastqSequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * 
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.np.filter", 
	scriptDesc = "Filter nanopore reads data from fastq file")
public class NanoporeReadFilter// implements Closeable
{
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/
		Deployable annotation = NanoporeReadFilter.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine(annotation.scriptName() + " [options]",
				annotation.scriptDesc());
		
		cmdLine.addStdInputFile();
		cmdLine.addStdOutputFile();
		cmdLine.addInt("lenMin", 0, 
				"Minimum sequence length");
		
		cmdLine.addInt("lenMax", Integer.MAX_VALUE, 
				"Minimum sequence length");
		
		cmdLine.addDouble("qualMin", 0, 
				"Minimum average quality");		

		cmdLine.addDouble("qualMax", 1000, 
				"Maximum average quality");
		
		cmdLine.addBoolean("excl2D", false, "Exclude 2D reads");
		cmdLine.addBoolean("exclTemp", false, "Exclude template reads");
		cmdLine.addBoolean("exclComp", false, "Exclude complement reads");
		
				
		cmdLine.addString("format", "fastq", "Format of the output file");
		args = cmdLine.stdParseLine_old(args);
		/**********************************************************************/

		String output = cmdLine.getStringVal("output");
		String input = cmdLine.getStringVal("input");
		
		
		int lenMin  = cmdLine.getIntVal("lenMin");
		int lenMax  = cmdLine.getIntVal("lenMax");
		
		double qualMin  = cmdLine.getDoubleVal("qualMin");
		double qualMax  = cmdLine.getDoubleVal("qualMax");
		
		boolean exclude2D =  cmdLine.getBooleanVal("exclude2D");
		boolean excludeTemplate =  cmdLine.getBooleanVal("excludeTemplate");
		boolean excludeComplement =  cmdLine.getBooleanVal("excludeComplement");
		
						
		String format = cmdLine.getStringVal("format");
		
		boolean fastaOutput = "fasta".equals(format.trim().toLowerCase());
		
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);
		
		FastqReader reader = "-".equals(input)? (new FastqReader(System.in) ) 
				: (new FastqReader(input));
		FastqSequence seq;
		
		
		
		while ( (seq = reader.nextSequence(Alphabet.DNA())) != null){
			
			//Min length
			if (seq.length() < lenMin)
				continue;
			
			//max length
			if (seq.length() > lenMax)
				continue;
			
			double qual = -1;
			
			//min quality
			if (qualMin > 0){
				qual = NanoporeReaderStream.averageQuality(seq);
				if (qual < qualMin)
					continue;				
			}
			
			//max quality			
			if (qualMax < 1000){
				if (qual < 0)
					qual = NanoporeReaderStream.averageQuality(seq);
				if (qual >= qualMax)
					continue;		
			}
			
			if (excludeComplement && seq.getName().contains("complement"))
				continue;
			
			if (excludeTemplate && seq.getName().contains("template"))
				continue;
			
			if (exclude2D && seq.getName().contains("twodim"))
				continue;
			
			//done all the fitlering
			if (fastaOutput)
				seq.writeFasta(sos);
			else
				seq.print(sos);
			
		}
		
		reader.close();		
		sos.close();
	}//main

}
