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

package japsa.seq.tools;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(scriptName = "jsa.seq.sort",
           scriptDesc = "Sort sequences based on their lengths")
public class SequenceSort {	
	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/
		Deployable annotation = SequenceSort.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options] ",
				annotation.scriptDesc());
		
		cmdLine.addStdInputFile();
		cmdLine.addStdOutputFile();		
		cmdLine.addStdAlphabet();
		cmdLine.addBoolean("number",false,"Add the order number to the beginning of contig name");
		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		//Get dna 		
		String alphabetOption = cmdLine.getStringVal("alphabet");		
		Alphabet alphabet = Alphabet.getAlphabet(alphabetOption);
		if (alphabet == null)
			alphabet = Alphabet.DNA16();

		String input = cmdLine.getStringVal("input");
		SequenceReader reader = SequenceReader.getReader(input);
		ArrayList<SequenceLength> seqList = new ArrayList<SequenceLength>(); 
		
		Sequence seq;
		while ((seq = reader.nextSequence(alphabet))!= null){
			seqList.add(new SequenceLength(seq));			
		}
		reader.close();		
		Collections.sort(seqList);
		Collections.reverse(seqList);
		
		String output = cmdLine.getStringVal("output");
		SequenceOutputStream sos = 	SequenceOutputStream.makeOutputStream(output);
		
		if (cmdLine.getBooleanVal("number")){
			for (int i = 0; i < seqList.size();  i++){
				String name = seqList.get(i).seq.getName();
				seqList.get(i).seq.setName(i + "-" + name);
			}	
		}
		for (int i = 0; i < seqList.size();  i++){
			seqList.get(i).seq.writeFasta(sos);
		}
		sos.close();
	}


	
	static class SequenceLength implements Comparable<SequenceLength>{
		Sequence seq;
		SequenceLength(Sequence seq){
			this.seq = seq;
		}
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(SequenceLength o) {
			// TODO Auto-generated method stub
			return seq.length() - o.seq.length();
		}
	}
}
