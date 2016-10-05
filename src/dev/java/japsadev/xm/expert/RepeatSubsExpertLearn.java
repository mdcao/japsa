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

package japsadev.xm.expert;

import japsa.seq.AbstractSequence;
import japsa.util.MyBitSet;
import japsa.xm.expert.Expert;
import japsa.xm.expert.RepeatExpert;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.StringTokenizer;

public class RepeatSubsExpertLearn extends RepeatExpert {

	public static double[][] countMatrix = { { 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0 } };

	public static int countSeq[] = { 0, 0, 0, 0 };
	public static int countBg[] = { 0, 0, 0, 0 };
	static int countSeqAll = 0, countBgAll = 0;

	// public static double[] dnaProb =
	// {.4043986778540554, 0.2972285786931096, 0.19946605644546148,
	// 0.09890668700737351};
	// {.25 ,.25 ,.25, .25};

	// #DEBUG_BEGIN
	public static int longestL = 0;
	public static double bestInfoGain = 0;
	// #DEBUG_END

	public static double[][] globalMatrix =
	// {{59.8817, 5.0548, 29.9753, 5.0882},{3.0111, 64.9124, 7.0311, 25.0454},
	// {29.9853, 4.0224, 59.9208, 6.0715}, {4.9642, 30.0038, 4.9152 , 60.116}};
	// {{40,10,10,10},{10,40,10,10},{10,10,40,10},{10,10,10,40}};
	{ { .50, .10, .30, .10 }, { .05, .60, .10, .25 }, { .20, .04, .70, .06 },
			{ .05, .05, .05, .80 } };

	// double[][] myMatrix;

	// double[] sum;

	double[][] myCount;

	// static double [] bgFreq = {0.1,0.2,.3,.4};
	// static double [] seqFreq = {0.25,0.25,.25,.25};

	public RepeatSubsExpertLearn(AbstractSequence seq, int start, MyBitSet b,
			int type) {

		super(seq, start, b, type);
		// Copy to my matrix
		// myMatrix = new
		// double[Expert.alphabet().size()][Expert.alphabet().size()];
		myCount = new double[Expert.alphabet().size()][Expert.alphabet().size()];

		// sum = new double[Expert.alphabet().size()];
		for (int i = 0; i < Expert.alphabet().size(); i++) {
			// sum[i] = 0.0;
			for (int j = 0; j < Expert.alphabet().size(); j++) {
				// myMatrix[i][j] = globalMatrix[i][j];
				// sum[i] += myMatrix[i][j];
				myCount[i][j] = 1.0;
			}
		}
	}

	public double probability(int character) {

		int match;
		if (expertType != PALIN_TYPE)
			match = character;
		else
			match = 3 - character;

		return globalMatrix[seq.symbolAt(currentPointer)][match];
		// myMatrix[japsa.seq[currentPointer]][match] /
		// sum[japsa.seq[currentPointer]];
	}

	public double update(int actual) {

		// if (encodeCount < countSize.length)
		// countSize[encodeCount] ++;

		// Move current pointer to the next char in knowledge
		if ((currentPointer - start) * expertType >= length) {
			return -1.0;
		}

		double prob = probability(actual);

		int match;
		if (expertType != PALIN_TYPE)
			match = actual;
		else
			match = 3 - actual;

		myCount[seq.symbolAt(currentPointer)][match] += 1.0;// weight /
		// Expert.repWeight;//
		// Expert.repWeight

		// myMatrix[japsa.seq[currentPointer]][match] += 1.0;
		// sum[japsa.seq[currentPointer]] += 1;

		currentPointer += expertType;
		updateCost(prob);
		return prob;

	}

	public static void readMatrix(BufferedReader in) throws IOException {
		int x = 0;
		String line = "";
		while ((line = in.readLine()) != null) {
			line = line.trim();
			if (line.startsWith("#"))
				continue;
			StringTokenizer tk = new StringTokenizer(line);
			for (int i = 0; i < 4; i++) {
				globalMatrix[x][i] = Double.parseDouble(tk.nextToken());
			}

			x++;
			if (x >= 4)
				break;
		}
		if (x < 4) {
			System.err.println("There are only " + x + " lines ");
		}
	}

	public static void printMatrix(PrintStream out) {
		for (int x = 0; x < RepeatSubsExpertLearn.globalMatrix.length; x++) {
			for (int y = 0; y < RepeatSubsExpertLearn.globalMatrix[x].length; y++) {
				out.printf("%6.4f  ", RepeatSubsExpertLearn.globalMatrix[x][y]);
			}
			out.println();
		}
	}

	public static void summary() {
		/** ****************************************************** */
		double countSum[] = { 0, 0, 0, 0 };
		double countSumSeq[] = { 0, 0, 0, 0 };
		for (int x = 0; x < RepeatSubsExpertLearn.countMatrix.length; x++) {
			for (int y = 0; y < RepeatSubsExpertLearn.countMatrix[x].length; y++) {
				System.out.printf("%6.2f  ",
						RepeatSubsExpertLearn.countMatrix[x][y]);
				countSum[x] += RepeatSubsExpertLearn.countMatrix[x][y];
				countSumSeq[y] += RepeatSubsExpertLearn.countMatrix[x][y];
			}
			System.out.println();
		}
		System.out
				.println("----------------------------------------------------------");

		// #DEBUG_BEGIN
		System.out.println(" Longest = " + longestL + "\n Best Infogain = "
				+ bestInfoGain);

		longestL = 0;
		bestInfoGain = 0;
		// #DEBUG_END

		System.out
				.println("----------------------------------------------------------");
		for (int y = 0; y < RepeatSubsExpertLearn.countMatrix.length; y++) {
			for (int x = 0; x < RepeatSubsExpertLearn.countMatrix[y].length; x++) {
				System.out.printf("%6.4f  ", 100
						* RepeatSubsExpertLearn.countMatrix[x][y]
						/ countSumSeq[y]);
			}
			System.out.println();
		}
		System.out
				.println("----------------------------------------------------------");

		for (int x = 0; x < RepeatSubsExpertLearn.countMatrix.length; x++) {
			for (int y = 0; y < RepeatSubsExpertLearn.countMatrix[x].length; y++) {
				// System.out.printf("%6.4f  ", 100
				// * RepeatSubsExpertLearn.countMatrix[x][y] / countSum[x]);

				RepeatSubsExpertLearn.globalMatrix[x][y] = (RepeatSubsExpertLearn.countMatrix[x][y] / countSum[x]);
				// * (bgFreq[x] / (countBg[x] * 1.0/countBgAll));
				countMatrix[x][y] = 0.0;
			}// normalise(RepeatSubsExpertLearn.globalMatrix[x]);
				// System.out.println();
		}
		printMatrix(System.out);
		// for (int x = 0; x < RepeatSubsExpertLearn.globalMatrix.length; x++) {
		// for (int y = 0; y < RepeatSubsExpertLearn.globalMatrix[x].length;
		// y++) {
		// System.out.printf("%6.4f  ",
		// RepeatSubsExpertLearn.globalMatrix[x][y]);
		// }
		// System.out.println();
		// }

		countSeqAll = countBgAll = 0;

		System.out
				.println("====================================================");
		/***********************************************************************
		 * for (int x = 0; x < RepeatExpert.countSize.length; x ++){ if
		 * (RepeatExpert.countSize[x] > 0) System.err.println(x+"\t" +
		 * JapsaMath.log2(RepeatExpert.countSize[x])); else
		 * System.err.println(x+"\t0.0"); }
		 * /*********************************************************
		 * 
		 * //double countProb = for (int x = 0; x < countProb.length; x ++){
		 * System.err.println((x * (COUNT_MAX-COUNT_MIN)/countProb.length)+"\t"
		 * + (countProb[x])); } /
		 **********************************************************************/
	}

	public static void normalise(double[] ar) {
		double sumA = 0;
		for (int i = 0; i < ar.length; i++) {
			sumA += ar[i];
		}
		for (int i = 0; i < ar.length; i++) {
			ar[i] /= sumA;
		}
	}

	public RepeatExpert duplicate(AbstractSequence seq_, int start_, MyBitSet b_) {
		return new RepeatSubsExpertLearn(seq_, start_, b_, expertType);
	}
}
