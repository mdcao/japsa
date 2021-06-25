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
 * 15/03/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.xm;

/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsadev.xm.expert.TandemExpert;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;

/**
 * A Tool for finding tandem repeats.
 * 
 * 
 */

public class TandemRepeatFinder {
	/**
	 * Determine if a sequence has a tandem repeat in it.
	 * 
	 * @param seqArray
	 *            the sequence
	 * @return an empty set if no repeat, else a set of the possible lengths of
	 *         the repeat unit.
	 */

	public static void main(String[] args) throws Exception {
		// Get params from users
		CommandLine cmdLine = new CommandLine();
		cmdLine.addInt("max", 8, "Maximum period size");
		cmdLine.addString("output", "-", "Output file");
		cmdLine.addInt("context", 10, "Length of the context");
		cmdLine.addDouble("threshold", 0.9, "Listen threshold");

		args = cmdLine.parseLine(args);

		if (args == null || args.length <= 0) {
			System.err
					.println("Usage: java bio.xm.TandemFinder [options] file\n"
							+ cmdLine.usageMessage() + "\n");
			System.exit(1);
		}

		int max = cmdLine.getIntVal("max");
		double threshold = cmdLine.getDoubleVal("threshold");

		ExpertModelTandem tModelF = new ExpertModelTandem(Alphabet.DNA5(),
				cmdLine.getIntVal("context"), threshold, max);

		// SequenceFileReader fReader
		SequenceReader fReader = SequenceReader.getReader(args[0]);

		Sequence s;

		SequenceOutputStream outBFF;

		String output = cmdLine.getStringVal("output");
		if (!output.equals("-")) {
			outBFF = new SequenceOutputStream(new FileOutputStream(output));
		} else {
			outBFF = new SequenceOutputStream(System.out);
		}
		while ((s = fReader.nextSequence(Alphabet.DNA5())) != null) {

			System.out.println(s.getName());
			tModelF.calculateTandemCosts(s);
			ExpertModelTandem.anno.writeAnnotation(outBFF);

			// tModelB.calculateTandemCosts(s.reverseComplement());

			/***************************************************/

			Iterator<TandemExpert> iter = tModelF.panel.iterator();
			for (int m = 0; m < max; m++) {
				TandemExpert exp = iter.next();
				PrintStream out = new PrintStream(new FileOutputStream(
						s.getName() + (m + 1) + ".info"));
				PrintStream out2 = new PrintStream(new FileOutputStream(
						s.getName() + (m + 1) + ".inf"));
				for (int i = 0; i < s.length(); i++) {
					out.println(exp.tandemCosts[i]);
					out2.println(exp.positive[i]);
				}
				out.close();
				out2.close();
			}
			/***************************************************/
		}
		outBFF.close();
	}
}
