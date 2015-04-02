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
 * 21/12/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.seq.tools;


import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.IOException;



/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
@Deployable(scriptName = "jsa.seq.join",
            scriptDesc = "Break a multiple sequence files to each sequence per file")
public class JoinSequenceFile {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		/*********************** Setting up script ****************************/
		Deployable annotation = JoinSequenceFile.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options] file1 file2 ...",
				annotation.scriptDesc());
				
		cmdLine.addStdAlphabet();
		
		cmdLine.addString("output", "-", "Name of the output file");
		cmdLine.addString("name", "name", "Name of the combined sequence");
		cmdLine.addBoolean("removeN", false, "Remove wildcards");
		
		//cmdLine.addString("format", "fasta", "Format of output files. Options : japsa or fasta");
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		//Get dna 		
		String alphabetOption = cmdLine.getStringVal("alphabet");		
		Alphabet alphabet = Alphabet.getAlphabet(alphabetOption);
		if (alphabet == null)
			alphabet = Alphabet.DNA();
		
		String output = cmdLine.getStringVal("output");
		String name = cmdLine.getStringVal("name");
		boolean removeN = cmdLine.getBooleanVal("removeN");
		
		//String format = cmdLine.getStringVal("format");		

		SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA(), 1000000, name);
		Sequence seq;
		for (String arg:args){
			SequenceReader reader = SequenceReader.getReader(arg);
			while ((seq = reader.nextSequence(alphabet)) != null){
				for (int i = 0; i < seq.length();i++){
					byte base =seq.getBase(i); 
					if ((!removeN) || base <4) 
						sb.append(base);
				}
				sb.setDesc(sb.getDesc() + ";" + seq.getDesc());
			}
			reader.close();			
		}				
		sb.writeFasta(output);		
	}
}
