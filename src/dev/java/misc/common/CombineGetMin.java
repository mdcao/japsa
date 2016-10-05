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

import java.io.*;

public class CombineGetMin {
	String inFile1, inFile2;

	public CombineGetMin(String in1, String in2) {
		this.inFile1 = in1;
		this.inFile2 = in2;
	}

	String readLn(BufferedReader in) throws IOException {
		String ln = "#";
		while (ln.startsWith("#")) {
			ln = in.readLine();
			if (ln == null)
				return null;
		}

		return ln;
	}

	// Read from stdin, smooth and write back to stdout
	public void combineMin() {
		try {
			// This approach is rather not memory efficient
			BufferedReader in1 = new BufferedReader(new FileReader(inFile1));
			BufferedReader in2 = new BufferedReader(new FileReader(inFile2));

			String line1, line2;
			int count = 0;
			while (((line1 = readLn(in1)) != null)
					&& (line2 = readLn(in2)) != null) {

				String arr1[] = line1.split(" |\t");
				String arr2[] = line2.split(" |\t");

				double value1 = Double.parseDouble(arr1[arr1.length - 1]);
				double value2 = Double.parseDouble(arr2[arr2.length - 1]);
				System.out.println(count + "\t" + Math.min(value1, value2));
				count++;

			}
			in1.close();
			in2.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length < 2) {
			System.err.println("Smooth inFile1 inFile2");
			System.exit(1);
		}

		CombineGetMin sm = new CombineGetMin(args[0], args[1]);
		sm.combineMin();
	}
}
