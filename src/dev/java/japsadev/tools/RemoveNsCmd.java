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
 * 28/05/2014 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsadev.tools;



import java.io.IOException;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * Filter a bam filem based on some criteria. Input file in bam format assumed
 * to be sorted and indexed
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.dev.removeNs", 
	scriptDesc = "Removes from fasta file"
	)
public class RemoveNsCmd extends CommandLine{
	//CommandLine cmdLine;
	public RemoveNsCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("input", "-", "Name of the input file");
		addString("output", "-", "Name of the output file");

		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		RemoveNsCmd cmdLine = new RemoveNsCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String input  = cmdLine.getStringVal("input");
		String output = cmdLine.getStringVal("output");

		SequenceReader reader = SequenceReader.getReader(input);
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);		
		Sequence seq;		
		while ((seq = reader.nextSequence(Alphabet.DNA16()))!=null){
			SequenceBuilder sb = removeNs(seq);
			sb.writeFasta(sos);

		}
		reader.close();
		sos.close();
	}


	public static SequenceBuilder removeNs(Sequence inSeq){
		SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA4(), inSeq.length(), inSeq.getName());
		for (int i = 0; i < inSeq.length();i++){
			byte base =inSeq.getBase(i); 
			if (base < 4) 
				sb.append(base);
		}
		sb.setDesc(inSeq.getDesc());
		return sb;
	}

}
