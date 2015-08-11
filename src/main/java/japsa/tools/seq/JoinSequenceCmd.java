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
package japsa.tools.seq;


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
@Deployable(
	scriptName = "jsa.seq.join",
	scriptDesc = "Join multiple sequences into one"
	)

public class JoinSequenceCmd extends CommandLine{	
	public JoinSequenceCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options] file1 file2 ...");
		setDesc(annotation.scriptDesc()); 
		
		addStdAlphabet();
		addString("output", "-", "Name of the output file, - for standard output");
		addString("name", "newseq", "Name of the new sequence");
		addBoolean("removeN", false, "Remove wildcards (N)");
		
		addStdHelp();		
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		/*********************** Setting up script ****************************/
		JoinSequenceCmd cmdLine = new JoinSequenceCmd();		
		args = cmdLine.stdParseLine(args);		
		/*********************************************************************/
		
		String alphabetOption = cmdLine.getStringVal("alphabet");		
		Alphabet alphabet = Alphabet.getAlphabet(alphabetOption);
		if (alphabet == null)
			alphabet = Alphabet.DNA();

		String output = cmdLine.getStringVal("output");
		String name = cmdLine.getStringVal("name");
		boolean removeN = cmdLine.getBooleanVal("removeN");

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

/*RST*
-----------------------------------------------------
*jsa.seq.join*: Join multiple sequences into one file 
-----------------------------------------------------


<usage> 

*RST*/