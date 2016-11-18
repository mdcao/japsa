/*****************************************************************************
 * Copyright (c) 2010 Minh Duc Cao, Monash University.  All rights reserved. *
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
 * 3. Neither the name of Monash University nor the names of its contributors*
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

package japsadev.misc.dnaPlatform.compModel;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import japsadev.misc.dnaPlatform.OptionsHandle;
import japsadev.misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: MarkovModel
 * </p>
 * 
 * <p>
 * Description: Class to use Markov model in DNAPlatform
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */

public class MarkovModel implements CompressionModel {

	public MarkovModel() {
	}

	/**
	 * Return all options needed to run Markov model in OptionsHandle object.
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		OptionsHandle myOps = new OptionsHandle(this, 2);
		myOps.addStringOption("dna", "atgc",
				"Alphabet used by the sequence.");
		myOps.addIntOption("order", 0, "Order of Markov Model");

		return myOps;
	}

	/**
	 * Returns a Class array with one element, DNA_SEQUENCE indicating this
	 * model can be applied to DNA sequences and any other type of char
	 * sequences.
	 * 
	 * @return Class[]
	 */
	@SuppressWarnings("rawtypes")
	public Class[] getTypeSequenceData() {
		Class[] types = { CharSequenceData.class };
		return types;
	}

	/**
	 * Runs markov model with inputSequenceData to create outputSequenceData.
	 * The output data for this model is a InfoContentSequenceData of type
	 * INFORMATION_CONTENT
	 * 
	 * @throws IOException
	 * @throws RuntimeException
	 */

	public SequenceData execute(OptionsHandle options, SequenceData data)
			throws RuntimeException {

		DoubleSequenceData outputData = new DoubleSequenceData(data);
		if (data instanceof CharSequenceData) {
			// get corresponding fields according to options
			String alphabet = options.getStringValue("dna");
			int order = options.getIntValue("order");

			Markov mark = new Markov(alphabet.length());

			// add this model to the history of outputData
			outputData.addHistory(options);

			// create double data in outputData by reading information content
			// file
			// created with markov
			outputData.setDoubleData(mark.calcInfoContent(
					((CharSequenceData) data).getCharData(), order));
		}
		return outputData;
	}

	/**
	 * Returns the name of the file where the information content calculated by
	 * model is stored. The name of the file returned is the filename given in
	 * MarkovDriver
	 * 
	 * @return String
	 */
	/*
	 * public String getInfoContentFile() { String prefix = (new
	 * File(sequenceFile)).getName(); return prefix + "_" +
	 * myOptions.getIntValue("order") + "MM-msgLen.txt"; }
	 */

	/**
	 * Returns null as Markov model does not create other files.
	 * 
	 * @return String[]
	 */

	public String[] getOtherGraphs() {
		return null;
	}

	public String toString() {
		return "Markov";
	}

}

/**
 * Markov implements a k-order Markov Model, the constructor for this model
 * takes as arguments the number of characters in sequence.
 * 
 * @author Julie Bernal
 * @version 1.0
 */
class Markov {

	private int myCharSet; // number of possible characters in myString
	private String nl;

	public Markov(int numCharSet) {
		myCharSet = numCharSet;
		nl = System.getProperty("line.separator");
	}

	/*
	 * printInfoContent: this class iterates through string s predicting the
	 * probability of every character according to a <order> Markov Model. The
	 * parameter order indicates the number of characters to consider before
	 * current character.
	 */
	public void printInfoContent(String s, final int order, String fileName) {
		try {

			BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
			HashMap<String, Integer> charCount = new HashMap<String, Integer>();
			String seq, totalSeq;
			int count = 1, totalCount = 1;

			double probability;// , bitsToCode = 0;

			// give probability of 1/charSet to all chars that don't have enough
			// <order> characters before them
			for (int i = 0; i < order; i++) {
				out.write(s.charAt(i) + "\t" + (-1 * log2(1.0 / myCharSet))
						+ nl);
				// bitsToCode += (-1 * log2(1.0/myCharSet));
			}

			for (int index = 0; index < s.length() - order; index++) {

				// 1. Look at current subsequence in s
				seq = s.substring(index, index + order + 1);
				totalSeq = seq.substring(0, seq.length() - 1) + "-";

				// Get counters to predict probability:
				// counter for all the characters
				if (charCount.containsKey(totalSeq))
					totalCount = ((Integer) charCount.get(totalSeq)).intValue();
				else
					// assume all chars seen once at the start
					totalCount = myCharSet;

				// conunter for current character
				if (charCount.containsKey(seq))
					count = ((Integer) charCount.get(seq)).intValue();
				else
					// assume all chars seen once at start
					count = 1;

				// 2. Predict probability of current character according to
				// preceeding sequence (InfoContent)
				probability = (double) count / (double) totalCount;

				// output information content of each character in sequence
				out.write(seq + "\t" + (-1 * log2(probability)) + nl);

				// DecimalFormat probFormat = new DecimalFormat("#0.00");
				// bitsToCode += (-1 * log2(probability));
				// System.out.println("\tBits: " + bitsToCode);

				// 3. Update counters:
				count++;
				totalCount++;
				charCount.put(seq, new Integer(count));
				charCount.put(totalSeq, new Integer(totalCount));
			}
			out.close();

		} catch (IOException e) {
		}
	}

	/**
	 * calcInfoContent: this function calculates the information content of a
	 * string, given an order.
	 * 
	 * @param s
	 *            String
	 * @param order
	 *            int
	 * @return double[] information content sequence
	 */

	public double[] calcInfoContent(char[] sequence, int order) {

		double[] infoCont = new double[sequence.length];
		int index = 0;

		HashMap<String, Integer> charCount = new HashMap<String, Integer>();
		String seq, totalSeq;
		int count = 1, totalCount = 1;

		double probability;

		String s = String.valueOf(sequence);

		// give probability of 1/charSet to all chars that don't have enough
		// <order> characters before them
		for (int i = 0; i < order; i++) {
			// out.write(s.charAt(i)+ "\t" + (-1*log2(1.0/myCharSet)) + nl);
			infoCont[index] = -1 * log2(1.0 / myCharSet);
			index++;
		}

		for (int i = 0; i < s.length() - order; i++) {
			// 1. Look at current subsequence in s
			seq = s.substring(i, i + order + 1);
			totalSeq = seq.substring(0, seq.length() - 1) + "-";

			// Get counters to predict probability:
			// counter for all the characters
			if (charCount.containsKey(totalSeq))
				totalCount = ((Integer) charCount.get(totalSeq)).intValue();
			else
				// assume all chars seen once at the start
				totalCount = myCharSet;

			// conunter for current character
			if (charCount.containsKey(seq))
				count = ((Integer) charCount.get(seq)).intValue();
			else
				// assume all chars seen once at start
				count = 1;

			// 2. Predict probability of current character according to
			// preceeding sequence (InfoContent)
			probability = (double) count / (double) totalCount;

			// output information content of each character in sequence
			infoCont[index] = (-1 * log2(probability));
			index++;

			// 3. Update counters:
			count++;
			totalCount++;
			charCount.put(seq, new Integer(count));
			charCount.put(totalSeq, new Integer(totalCount));
		}

		return infoCont;
	}

	/**
	 * Calculates log2 of argument given.
	 * 
	 * @param n
	 *            double
	 * @return double
	 */
	private double log2(double n) {
		return Math.log(n) / Math.log(2);
	}
}
