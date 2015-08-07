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
 * 26/07/2015 - minhduc reorganised to tools                                  
 ****************************************************************************/

package japsa.tools.seq;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.IOException;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.seq.stats",
	scriptDesc = "Show statistical composition of sequences stored in a file (or from STDIN)",
	scriptDocs = ""
		+ ""
	)
public class SequenceStats extends CommandLine{	
	public SequenceStats(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addStdInputFile();
		addStdAlphabet();//aphabet
		
		addStdHelp();		
	}
	public static void main(String[] args) throws IOException {				
		CommandLine cmdLine = new SequenceStats();
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		//Get dna 		
		String alphabetOption = cmdLine.getStringVal("alphabet");		
		Alphabet alphabet = Alphabet.getAlphabet(alphabetOption);
		if (alphabet == null)
			alphabet = Alphabet.DNA16();

		String input = cmdLine.getStringVal("input");
		/**********************************************************************/	


		SequenceReader reader = SequenceReader.getReader(input);
		long total = 0;
		int numSeq = 0;
		Sequence seq;
		while ((seq = reader.nextSequence(alphabet))!= null){
			total += seq.length();
			numSeq ++;
			System.out.println(seq.getName() + " :  " + seq.length() + " bases");
			System.out.println(seq.getDesc());
			getComposition(seq);
		}
		reader.close();
		System.out.println("Total = " + total + " bases in " + numSeq + " sequences.");
	}


	/**
	 * Get properties of the sequences in the file
	 * 
	 * @param args	 
	 */
	static private void getComposition(Sequence seq) {
		if (seq.length() == 0) {
			System.out.println("Sequence " + seq + " contains 0 base/residue");
			return;
		}
		Alphabet alphabet = seq.alphabet();
		int[] counts = new int[alphabet.size()];
		int others = 0;

		for (int i = 0; i < seq.length(); i++) {
			int index = seq.symbolAt(i);
			if (index < 0 || index >= counts.length)
				others++;
			else
				counts[index]++;
		}
		for (int index = 0; index < counts.length; index++) {
			if (counts[index] > 0)
				System.out.printf("%10d  %c : %5.2f%%\n", counts[index], alphabet
					.int2char(index), (counts[index] * 100.0 / seq.length()));
		}
		if (others > 0)
			System.out.printf("%10d  %c : %5.2f%%\n", others, 'X',
				(others * 100.0 / seq.length()));

	}

}
