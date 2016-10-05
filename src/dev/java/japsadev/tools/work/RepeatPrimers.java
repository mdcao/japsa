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
 * 08/04/2012 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsadev.tools.work;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.dev.primerTR", scriptDesc = "Design primers for repeats. Each primer is added 30 bases for mapping")
public class RepeatPrimers {
	public static void main(String[] args) throws Exception {

		/*********************** Setting up script ****************************/
		Deployable annotation = RepeatPrimers.class
				.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options]",
				annotation.scriptDesc());

		cmdLine.addStdInputFile();
		cmdLine.addString("tr", null, "Name of the tandem repeat file", true);
		cmdLine.addInt("flanking", 300, "Length of flanking regions");
		cmdLine.addInt("pad", 30,
				"Pad to the primers for more specific mappings");

		args = cmdLine.stdParseLine_old(args);
		/**********************************************************************/

		String input = cmdLine.getStringVal("input");
		int flanking = cmdLine.getIntVal("flanking");
		int pad = cmdLine.getIntVal("pad");
		// int pad = 30;
		int gap = 5 + pad;

		// Reader sequences to the hash
		SequenceReader reader = SequenceReader.getReader(input);
		HashMap<String, Sequence> genomes = new HashMap<String, Sequence>();
		Sequence seq;

		String primerFile = "pinput";
		SequenceOutputStream pos = SequenceOutputStream
				.makeOutputStream(primerFile);
		Sequence wSeq = null;
		while ((seq = reader.nextSequence(Alphabet.DNA16())) != null) {
			genomes.put(seq.getName(), seq);
			if (wSeq == null)
				wSeq = seq;
		}
		reader.close();
		seq = wSeq;

		// Reader in tr
		BufferedReader bf = SequenceReader.openFile(cmdLine.getStringVal("tr"));
		ArrayList<TandemRepeat> trList = TandemRepeat.readFromFile(bf, null);
		bf.close();

		// Get sequence information
		for (TandemRepeat tr : trList) {
			if (!seq.getName().equals(tr.getChr())) {
				seq = genomes.get(tr.getChr());
			}
			if (seq == null) {
				throw new RuntimeException("Sequence " + tr.getChr()
						+ " not found");
			}

			int start = tr.getStart() - flanking;
			int end = tr.getEnd() + flanking;

			if (start <= 0)
				start = 1;

			if (end > seq.length())
				end = seq.length();

			int primerStart = tr.getStart() - start - gap;// to make sure
			int primerLength = tr.getLength() + gap * 2;

			pos.print("SEQUENCE_ID=" + tr.getID() + "#" + tr.getParent() + "#"
					+ start + "\n");
			pos.print("SEQUENCE_TEMPLATE=");
			for (int i = start - 1; i < end; i++) {
				pos.print(seq.charAt(i));
			}
			pos.print("\nSEQUENCE_TARGET=" + primerStart + "," + primerLength
					+ "\n=\n");
		}
		pos.close();
		// Run primer 3

		Process process = Runtime
				.getRuntime()
				.exec("primer3_core -p3_settings_file=settings -output=primer3.out -error=primer.error pinput");
		process.waitFor();
		// Reat the output
		bf = SequenceReader.openFile("primer3.out");

		SequenceOutputStream fq1 = SequenceOutputStream
				.makeOutputStream("1.fq");
		SequenceOutputStream fq2 = SequenceOutputStream
				.makeOutputStream("2.fq");

		String name = "";
		int index = 0;
		int seqOffset = 0;
		int leftPos = 0, rightPos = 0;
		String left = "", right = "";

		String line = "";
		while ((line = bf.readLine()) != null) {

			String[] toks = line.trim().split("=");
			// System.out.println(line + toks.length);
			if (toks.length < 2)
				System.out.println();// end of record
			else if (toks[0].equals("SEQUENCE_ID")) {
				System.out.print(toks[1]);
				int x = toks[1].lastIndexOf('#');
				name = toks[1].substring(0, x);

				seqOffset = Integer.parseInt(toks[1].substring(x + 1));
				index = 0;
			} else if (toks[0].equals("PRIMER_LEFT_NUM_RETURNED")
					|| toks[0].equals("PRIMER_RIGHT_NUM_RETURNED")
					|| toks[0].equals("PRIMER_PAIR_NUM_RETURNED")) {
				System.out.print("\t" + toks[1]);
			} else if (toks[0].equals("PRIMER_LEFT_" + index + "_SEQUENCE")) {
				left = toks[1];
				System.out.print("\t" + toks[1]);
			} else if (toks[0].equals("PRIMER_RIGHT_" + index + "_SEQUENCE")) {
				right = toks[1];
				System.out.print("\t" + toks[1]);
			} else if (toks[0].equals("PRIMER_LEFT_" + index)) {
				int x = toks[1].indexOf(',');
				leftPos = seqOffset + Integer.parseInt(toks[1].substring(0, x))
						- 1;
			} else if (toks[0].equals("PRIMER_RIGHT_" + index)) {
				int x = toks[1].indexOf(',');
				rightPos = seqOffset
						+ Integer.parseInt(toks[1].substring(0, x))
						- right.length();
			} else if (toks[0].equals("PRIMER_PAIR_" + index + "_PRODUCT_SIZE")) {
				/****************************************************************************
				 * String readName = name.replaceFirst("#", "#"+index+"#") + "#"
				 * + leftPos + "#"+rightPos + "#"+toks[1];
				 * 
				 * fq1.print("@"+readName+"\n"+left+"\n+\n"); for (int i = 0; i<
				 * left.length(); i++) fq1.print('I'); fq1.print('\n');
				 * 
				 * fq2.print("@"+readName+"\n"+right+"\n+\n");
				 * //fq2.print("@"+name
				 * +"#"+index+"#"+toks[1]+"\n"+right+"\n+\n"); for (int i = 0;
				 * i< right.length(); i++) fq2.print('I'); fq2.print('\n'); /
				 ****************************************************************************/
				int x = name.indexOf('#');
				String seqID = name.substring(x + 1);
				if (!seq.getName().equals(seqID)) {
					seq = genomes.get(seqID);
				}
				if (seq == null) {
					throw new RuntimeException("Sequence " + seqID
							+ " not found");
				}
				rightPos -= pad;
				String readName = name.replaceFirst("#", "#" + index + "#")
						+ "#" + leftPos + "#" + rightPos + "#" + toks[1];
				fq1.print("@" + readName + "\n");

				for (int i = 0; i < left.length() + pad; i++) {
					fq1.print(seq.charAt(i + leftPos - 1));
				}
				fq1.print("\n+\n");

				for (int i = 0; i < left.length() + pad; i++) {
					fq1.print('I');
				}
				fq1.print('\n');

				Alphabet.DNA alphabet = (Alphabet.DNA) seq.alphabet();
				fq2.print("@" + readName + "\n");
				for (int i = right.length() + pad - 1; i >= 0; i--) {
					fq2.print(alphabet.int2char(alphabet.complement(seq
							.getBase(i + rightPos - 1))));
				}
				fq2.print("\n+\n");

				for (int i = right.length() + pad - 1; i >= 0; i--) {
					fq2.print('I');
				}
				fq2.print('\n');

				System.out.print("\t" + toks[1]);
				index++;
			}
		}// while
		bf.close();
		fq1.close();
		fq2.close();
		// Run bwa to make sure
		/*******************************************************************
		 * String cmd = "bwa aln " + input + " 1.fq > 1.sai"; process =
		 * Runtime.getRuntime().exec(cmd);
		 * 
		 * process.getErrorStream().close(); process.getInputStream().close();
		 * 
		 * process.waitFor(); Logging.info("Running " + cmd);
		 * 
		 * 
		 * cmd = "bwa aln " + input + " 2.fq > 2.sai"; process =
		 * Runtime.getRuntime().exec(cmd); process.getErrorStream().close();
		 * process.getInputStream().close(); process.waitFor();
		 * Logging.info("Running " + cmd);
		 * 
		 * cmd = "bwa sampe -f a.sam " + input + " 1.sai 2.sai 1.fq 2.fq";
		 * process = Runtime.getRuntime().exec(cmd);
		 * process.getErrorStream().close(); process.getInputStream().close();
		 * process.waitFor(); Logging.info("Running " + cmd); /
		 *******************************************************************/
	}
}
