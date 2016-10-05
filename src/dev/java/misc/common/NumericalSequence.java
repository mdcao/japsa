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

package misc.common;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.BitSet;
import java.util.Vector;

/**
 * @author Minh Duc Cao this class represent a nummerical sequence and various
 *         operations on it.
 */
public class NumericalSequence {
	private double[] seqs;

	/**
	 * Create a numerical sequence that hold the data
	 * 
	 * @param aSeqs
	 */
	public NumericalSequence(double[] aSeqs) {
		seqs = new double[aSeqs.length];
		for (int index = 0; index < seqs.length; index++) {
			seqs[index] = aSeqs[index];//
		}
	}

	public double[] getSeqs() {
		return seqs;
	}

	public static double[] read(String fileName) throws IOException {

		BufferedReader bufRdr = new BufferedReader(new FileReader(fileName));
		String line = null;
		Vector<Double> v = new Vector<Double>();
		while ((line = bufRdr.readLine()) != null) {
			if (line.startsWith("#"))
				continue;// Comment lines
			String[] tokens = line.trim().split("\\s");// Broken up

			v.add(Double.parseDouble(tokens[tokens.length - 1]));// .. and get
																	// the last
																	// token
		}

		bufRdr.close();
		// Convert from vector to array
		double[] data = new double[v.size()];
		for (int i = 0; i < v.size(); i++)
			data[i] = v.get(i);

		return data;
	}

	public double[] smooth(int wSize) {
		double[] out = new double[seqs.length];
		double[] his = new double[wSize];
		int index = 0;
		double sum = 0.0;

		for (int i = 0; i < wSize; i++) {
			his[i] = 0.0;
		}

		for (int i = 0; i < seqs.length; i++) {
			index = i % wSize;
			sum = sum - his[index] + seqs[i];
			his[index] = seqs[i];

			if (i < wSize)
				out[i / 2] = sum / (i + 1);
			else
				out[i - wSize / 2] = sum / wSize;
		}

		for (int i = wSize - 1; i > 0; i--) {
			index = (index + 1) % wSize;
			sum -= his[index];
			out[out.length - i / 2 - 1] = sum / (i);
		}
		return out;
	}

	public double getSum(BitSet bs) {
		double sum = 0;
		for (int index = 0; index < seqs.length; index++) {
			if (bs.get(index))
				sum += seqs[index];
		}

		return sum;
	}

	public double getSum() {
		double sum = 0;
		for (int index = 0; index < seqs.length; index++) {
			sum += seqs[index];
		}
		return sum;
	}

	NumericalSequence difference(NumericalSequence aSeq) {
		// This may throw exception if the sizes are not match
		double[] anoSeq = aSeq.seqs;
		double[] newSeq = new double[seqs.length];

		for (int i = 0; i < newSeq.length; i++) {
			newSeq[i] = seqs[i] - anoSeq[i];
		}
		return new NumericalSequence(newSeq);
	}

	public boolean writeDataToFile(File file) {
		try {
			PrintWriter pw = new PrintWriter(new FileOutputStream(file));
			pw.println("# Double data written by DNAGraphTool");
			for (int i = 0; i < seqs.length; i++) {
				pw.println(i + "\t" + seqs[i]);
			}
			pw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return true;
	}

	static double[] difference(double[] seq1, double[] seq2) {
		double[] seqResult = new double[seq1.length];

		for (int i = 0; i < seqResult.length; i++) {
			seqResult[i] = seq1[i] - seq2[i];
		}

		return seqResult;
	}

	public static void main(String[] args) {
		try {
			double[] seq0 = read(args[0]), seq1 = read(args[1]);

			NumericalSequence output = new NumericalSequence(difference(seq0,
					seq1));
			output.writeDataToFile(new File(args[2]));

		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
