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
import japsa.seq.FastqSequence;
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
public class SequenceSortCmd extends CommandLine{	
	public SequenceSortCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addStdInputFile();
		addStdOutputFile();		
		addStdAlphabet();
		addBoolean("number",false,"Add the order number to the beginning of contig name");
		addBoolean("reverse",false,"Reverse sort order");
		addString("sortKey","length","Sort key");

		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new SequenceSortCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		String inputOption = cmdLine.getStringVal("input");
		String sortKeyOption = cmdLine.getStringVal("sortKey");
		String outputOption = cmdLine.getStringVal("output");
		boolean numberOption =  cmdLine.getBooleanVal("number");
		boolean reverseOption =  cmdLine.getBooleanVal("reverse");


		//Get dna 		
		String alphabetOption = cmdLine.getStringVal("alphabet");		
		Alphabet alphabet = Alphabet.getAlphabet(alphabetOption);
		if (alphabet == null)
			alphabet = Alphabet.DNA16();



		SequenceReader reader = SequenceReader.getReader(inputOption);
		ArrayList<SequenceLength> seqList = new ArrayList<SequenceLength>(); 

		Sequence seq;
		String sortKeyOptionPrefix = sortKeyOption + "=";
		int sortKeyOptionIndex = sortKeyOptionPrefix.length();
		while ((seq = reader.nextSequence(alphabet))!= null){
			String [] toks = seq.getName().split(" ");

			SequenceLength seqL = new SequenceLength(seq);
			if (sortKeyOption.equals("length"))
				seqL.keyCompare = seq.length();
			else{
				for (int i = 0; i < toks.length;i++){
					if (toks[i].startsWith(sortKeyOptionPrefix)){
						double t = Double.parseDouble(toks[i].substring(sortKeyOptionIndex));
						seqL.keyCompare = (long) t;						
					}	
				}
			}

			seqList.add(seqL);			
		}
		reader.close();		
		Collections.sort(seqList);

		if (reverseOption)
			Collections.reverse(seqList);


		SequenceOutputStream sos = 	SequenceOutputStream.makeOutputStream(outputOption);


		for (int i = 0; i < seqList.size();  i++){
			seq = seqList.get(i).seq;
			if (numberOption){
				seq.setName(i + "-" + seq.getName());
			}
			if (seq instanceof FastqSequence){
				FastqSequence fq = ((FastqSequence) seq);
				fq.print(sos);				
			}else				
				seq.writeFasta(sos);
		}	
		sos.close();
	}



	static class SequenceLength implements Comparable<SequenceLength>{
		Sequence seq;
		long keyCompare;
		SequenceLength(Sequence seq){
			this.seq = seq;
		}
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(SequenceLength o) {
			// TODO Auto-generated method stub
			return (int) (keyCompare - o.keyCompare);
		}
	}
}


/*RST*
--------------------------------------------
*jsa.seq.sort*: Sort the sequences in a file
--------------------------------------------

*jsa.seq.sort* sort the sequences from a file or from a standard input into
some order.

*jsa.seq.sort* is included in the 
`Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

<usage> 

*RST*/
