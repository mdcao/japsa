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

package japsa.tools.seq;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.IOException;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(scriptName = "jsa.seq.rev",
           scriptDesc = "Reverse complement sequences (must be DNA)")
public class SequenceReverseComplement extends CommandLine {	
	public SequenceReverseComplement(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 
		
		addStdInputFile();
		addStdOutputFile();		
		addStdAlphabet();//aphabet
		
		addStdHelp();
	}
	
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new SequenceReverseComplement();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/	
		
		Alphabet alphabet = Alphabet.getAlphabet(cmdLine.getStringVal("alphabet"));
		if (alphabet == null)
			alphabet = Alphabet.DNA16();

		String input  = cmdLine.getStringVal("input");
		String output = cmdLine.getStringVal("output");
		/**********************************************************************/
		
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);
		SequenceReader reader = SequenceReader.getReader(input);		
		
		
		Sequence seq;
		while ((seq = reader.nextSequence(alphabet))!= null){		
			Sequence rseq = Alphabet.DNA.complement(seq);
			rseq.setName(seq.getName()+"_rev");
			rseq.writeFasta(sos);
		}
		sos.close();
		reader.close();		
		
	}
}
