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

package japsa.bio.misc.dnaPlatform.function;

import java.io.*;
import java.util.Vector;

public class Smooth {
	String inFile, outFile;
	int wSize;

	public Smooth(int wSize, String inFile, String outFile) {
		this.inFile = inFile;
		this.outFile = outFile;
		this.wSize = wSize;
	}

	/**
	 * Smooth the data
	 * 
	 * @param in
	 * @param wSize
	 * @return
	 */

	public static double[] smooth(double[] in, int wSize) {
		double[] out = new double[in.length];
		double[] his = new double[wSize];
		int index = 0;
		double sum = 0.0;

		for (int i = 0; i < wSize; i++) {
			his[i] = 0.0;
		}

		for (int i = 0; i < in.length; i++) {
			index = i % wSize;
			sum = sum - his[index] + in[i];
			his[index] = in[i];

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

	// Read from stdin, smooth and write back to stdout
	public void smooth() {
		try {
			// This approach is rather not memory efficient
			BufferedReader in = new BufferedReader(new FileReader(inFile));
			Vector<Double> v = new Vector<Double>();
			String line;
			while ((line = in.readLine()) != null) {
				String arr[] = line.split(" |\t");
				double value = Double.parseDouble(arr[arr.length - 1]);
				v.add(value);
			}
			in.close();
			// Make an array of double
			double[] data = new double[v.size()];
			for (int i = 0; i < v.size(); i++)
				data[i] = v.get(i);

			double[] outData = smooth(data, wSize);

			PrintStream out = new PrintStream(new FileOutputStream(outFile));
			for (int i = 0; i < v.size(); i++) {
				out.println(i + "\t" + outData[i]);
			}

			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length < 3) {
			System.err.println("Smooth smooth inFile outFile");
			System.exit(1);
		}
		int s = Integer.parseInt(args[0]);
		Smooth sm = new Smooth(s, args[1], args[2]);
		sm.smooth();

		// double[] in = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
		// double[] out = smooth(in,1);
		// for (int i = 0;i < out.length; i++){
		// System.out.println(out[i]);
		// }
	}
}
