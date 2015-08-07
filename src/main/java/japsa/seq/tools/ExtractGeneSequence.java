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
import japsa.seq.JapsaAnnotation;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.IOException;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(scriptName = "jsa.seq.gff2fasta",
           scriptDesc = "Extract sequences from a gff annotation")
public class ExtractGeneSequence {	
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/
		Deployable annotation = ExtractGeneSequence.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options] ",
				annotation.scriptDesc());
		
		cmdLine.addString("sequence", null, "The sequence (whole chromosome)",true);
		cmdLine.addString("gff", null, "Annotation file in gff format",true);
		cmdLine.addString("type", "gene", "types of features to be extracted (all, gene, CDS etc)");
		cmdLine.addStdOutputFile();		
		
		cmdLine.addStdAlphabet();//aphabet		
		
		args = cmdLine.stdParseLine_old(args);
		/**********************************************************************/
		//Get dna 		
		String alphabetOption = cmdLine.getStringVal("alphabet");		
		Alphabet alphabet = Alphabet.getAlphabet(alphabetOption);
		if (alphabet == null)
			alphabet = Alphabet.DNA16();

		String sequence = cmdLine.getStringVal("sequence");
		String gff = cmdLine.getStringVal("gff");
		String type = cmdLine.getStringVal("type");
		String output = cmdLine.getStringVal("output");
		/**********************************************************************/	
		

		SequenceReader reader = SequenceReader.getReader(sequence);
		Sequence seq = reader.nextSequence(alphabet);
		reader.close();
		String name  = seq.getName();
		
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(output);
		
		//a quick hack
		String [] toks = name.split("\\|");
		if (toks.length>3 && toks[2].equals("ref"))
			seq.setName(toks[3]);
		
		
		//Read annotation, without upstream and downstream
		BufferedReader aReader = SequenceReader.openFile(gff);
		JapsaAnnotation anno = JapsaAnnotation.readGFF(aReader,0,0,type);
		aReader.close();
		
		if (!anno.getAnnotationID().equals(seq.getName())){
			Logging.exit( "IDs dont match : " + seq.getName() + " vs " + anno.getAnnotationID(),1);
		}
		anno.setSequence(seq);		
		anno.writeFeatureSequence(out);
		out.close();
	}
}
