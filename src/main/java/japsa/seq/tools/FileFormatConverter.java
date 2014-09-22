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
 * 08/01/2012 - Minh Duc Cao: Revised                                        
 *  
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



/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.seq.format",
            scriptDesc = "Convert sequence(s) from a format to another")
public class FileFormatConverter {
	/**
	 * 
	 * @param args
	 */

	public static void main(String[] args) throws Exception {
		/*********************** Setting up script ****************************/
		Deployable annotation = FileFormatConverter.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options] ",
				annotation.scriptDesc());
		
		cmdLine.addStdInputFile();
		cmdLine.addString("output", "-", "Name of the output file,  - for standard output");
		cmdLine.addString("format", "fasta", "Format of output file. Options : japsa, fasta and phylip");
		cmdLine.addStdAlphabet();
		
		args = cmdLine.stdParseLine(args);	
		/********************************************************************/
		String output = cmdLine.getStringVal("output");
		String format = cmdLine.getStringVal("format");
		String input = cmdLine.getStringVal("input");

		Alphabet alphabet = Alphabet.getAlphabet(cmdLine.getStringVal("alphabet"));
		if (alphabet == null)
			alphabet = Alphabet.DNA5();
		/********************************************************************/

		SequenceReader reader = SequenceReader.getReader(input);
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream((output));			

		Sequence seq;
		//if (format.equals("nexus")){
		//	Sequence [] seqHash = new Sequence[array.size()];
		//	for (int i = 0; i < seqHash.length; i++)
		//		seqHash[i] = array.get(i);
		//	PhylogenyTree.writeNexus(seqHash,out);
		//}else	
		if (format.equals("phylip")){
			ArrayList<Sequence> seqs = new ArrayList<Sequence>(); 
			while ((seq = reader.nextSequence(alphabet)) != null){
				seqs.add(seq);
			}
			writePhylip(seqs,out);
		}else{ 
			if (format.equals("fasta")){
				while ((seq = reader.nextSequence(alphabet)) != null){
					seq.writeFasta(out);
				}
			}else {//if (outType.equals("jsa")){
				while ((seq = reader.nextSequence(alphabet)) != null){
					seq.writeJSA(out);
				}
			}  
		}

		out.close();
	}


	public static void writePhylip(ArrayList<Sequence> seqs, SequenceOutputStream out)
			throws IOException {

		Sequence seq = seqs.get(0);
		int length = seq.length();
		int charPerLine = length + 12;

		out.print(seqs.size()+ "   " + length+"\n");
		int count = 0;

		while (true) {
			for (int i = 0; i < seqs.size(); i++) {
				seq = seqs.get(i);
				if (count == 0) {
					out.print((seq.getName() + "             ").substring(0, 10) + "  ");
				}

				for (int x = count; x < count + charPerLine && x < length; x++) {
					if (x % 10 == 0 && x > count)
						out.print(' ');
					out.print(seq.charAt(x));
				}
				out.print('\n');
			}

			out.print('\n');

			count += charPerLine;
			if (count >= length)
				break;
		}
		//out.flush();		
	}
}
