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

package japsa.tools.bio.hts;



import java.io.IOException;
import java.util.Arrays;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.IntArray;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(
		scriptName = "jsa.hts.n50", 
		scriptDesc = "Compute N50 of an assembly"
		)
public class GetN50Cmd extends CommandLine{
	//CommandLine cmdLine;
	public GetN50Cmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("input", null, "Name of the file",true);

		addStdHelp();
	}
	public static void main(String [] args) throws IOException, InterruptedException{
		GetN50Cmd cmdTool = new GetN50Cmd ();
		args = cmdTool.stdParseLine(args);

		/**********************************************************************/
		String input = cmdTool.getStringVal("input");


		//ArrayList<Sequence> seqs = SequenceReader.readAll(input, Alphabet.DNA());
		Alphabet dna = Alphabet.DNA();
		//double n50 = HTSUtilities.n50(seqs);
		//System.out.println(n50 + "\t" + seqs.size());

		IntArray lengthArray = new IntArray();
		SequenceReader reader = SequenceReader.getReader(input);

		Sequence seq = null;
		while ((seq = reader.nextSequence(dna))!= null){
			lengthArray.add(seq.length());
		}

		int [] lengths = lengthArray.toArray();

		double sum = 0;
		for (int i = 0;i < lengths.length;i++){			
			sum += lengths[i];
		}		
		Arrays.sort(lengths);

		int index = lengths.length;
		double contains = 0;
		while (contains < sum/2){
			index --;
			contains += lengths[index];
		}

		int n50 = lengths[index]; 
		System.out.println(n50 + "\t" + lengths.length + "\t" + sum);

	}	
}

/*RST*
-----------------------------------------
 *jsa.hts.n50*: Compute N50 of an assembly 
-----------------------------------------

<usage> 

 *RST*/

