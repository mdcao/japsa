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
 * 02/04/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.tools.work;

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFileFormat;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;

import java.io.IOException;

/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 * 
 */
public class CompareAnnotation {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/
		String scriptName = "japsa.seq compareAnno";
		String desc = "Compare an annotation against a model one\n";
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName
				+ " [params]");
		/**********************************************************************/

		cmdLine.addString("anno", "-",
				"Name of the annotation file,  - for standard output");
		cmdLine.addString("model", null, "Name of the model annotation", true);

		cmdLine.addString("japsa.seq", null,
				"Name of the sequence file, if needed");
		cmdLine.addStdHelp();

		/**********************************************************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")) {
			System.out.println(desc + cmdLine.usageMessage());
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usageMessage());
			System.exit(-1);
		}
		/**********************************************************************/

		String annoFile = cmdLine.getStringVal("anno");
		String modelFile = cmdLine.getStringVal("model");
		String seqFile = cmdLine.getStringVal("japsa.seq");

		// Alphabet dna = Alphabet.getAlphabet(alphabetName);
		JapsaFileFormat reader;
		if (annoFile.equals("-"))
			reader = new JapsaFileFormat(System.in);
		else
			reader = new JapsaFileFormat(annoFile);

		JapsaAnnotation anno = reader.readAnnotation();
		reader.close();

		reader = new JapsaFileFormat(modelFile);
		JapsaAnnotation model = reader.readAnnotation();
		reader.close();

		if (model.getSequence() != null && anno.getSequence() == null)
			anno.setSequence(model.getSequence());

		if (model.getSequence() == null && anno.getSequence() != null)
			model.setSequence(anno.getSequence());

		if (seqFile != null) {
			// FIXME: Make sure seqFile is not -
			SequenceReader seqReader = SequenceReader.getReader(seqFile);
			// FIXME: checl dna
			Sequence seq = seqReader.nextSequence(null);
			if (seq != null) {
				model.setSequence(seq);
				anno.setSequence(seq);
			}
		}

		System.out.println("At feature level ");
		anno.compareAnnotationAtFeatureLevel(model, 1);

		System.out.println("\nAt base level ");
		anno.compareAnnotationAtBaseLevel(model);
	}
}
