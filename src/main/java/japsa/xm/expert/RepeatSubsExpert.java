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

package japsa.xm.expert;

import japsa.seq.AbstractSequence;
import japsa.util.MyBitSet;

import java.io.*;




public class RepeatSubsExpert extends RepeatExpert {
	// #DEBUG_BEGIN
	// public static int goodCount = 0, badCount = 0;
	// public static int longestL = 0;
	// #DEBUG_END

	public static double[][] globalMatrix = { { 40, 10, 10, 10 },
			{ 10, 40, 10, 10 }, { 10, 10, 40, 10 }, { 10, 10, 10, 40 } };

	double[][] myMatrix;

	double[] sum;

	public RepeatSubsExpert(AbstractSequence seq, int start, MyBitSet b, int type) {

		super(seq, start, b, type);
		// Copy to my matrix
		myMatrix = new double[Expert.alphabet().size()][Expert.alphabet().size()];

		sum = new double[Expert.alphabet().size()];
		for (int i = 0; i < Expert.alphabet().size(); i++) {
			sum[i] = 0.0;
			for (int j = 0; j < Expert.alphabet().size(); j++) {
				myMatrix[i][j] = globalMatrix[i][j];
				sum[i] += myMatrix[i][j];
			}
		}
	}

	public double probability(int character) {

		int match;
		if (expertType != PALIN_TYPE)
			match = character;
		else
			match = 3 - character;

		return myMatrix[seq.symbolAt(currentPointer)][match] / sum[seq.symbolAt(currentPointer)];
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

		myMatrix[seq.symbolAt(currentPointer)][match] += 1.0;
		sum[seq.symbolAt(currentPointer)] += 1;

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
			String[] tks = line.split(" +");
			for (int i = 0; i < 4; i++) {
				globalMatrix[x][i] = 100 * Double.parseDouble(tks[i]);
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
		for (int x = 0; x < globalMatrix.length; x++) {
			for (int y = 0; y < globalMatrix[x].length; y++) {
				out.printf("%6.4f  ", globalMatrix[x][y]);
			}
			out.println();
		}
	}

	public RepeatExpert duplicate(AbstractSequence seq_, int start_, MyBitSet b_) {
		return new RepeatSubsExpert(seq_, start_, b_, expertType);
	}
}
