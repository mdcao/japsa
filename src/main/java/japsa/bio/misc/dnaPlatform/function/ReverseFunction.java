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

//This class is written by Julie Bernal and subsequently modified and maintained
//by Minh Duc Cao

package japsa.bio.misc.dnaPlatform.function;

import japsa.bio.misc.dnaPlatform.OptionsHandle;
import japsa.bio.misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: ReverseFunction
 * </p>
 * 
 * <p>
 * Description: This class is used to reverse SequenceData.
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public class ReverseFunction implements Function {

	public ReverseFunction() {
	}

	public Class[] getTypeSequenceData() {
		Class[] typeData = { SequenceData.class };
		return typeData;
	}

	/**
	 * This method retusrns null as this function does not have any options.
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		return null;
	}

	/**
	 * Creates a new SequenceData from input SequenceData, reverses new
	 * SequenceData and adds this function to its history.
	 * 
	 * @param seqData
	 *            SequenceData
	 * @return SequenceData
	 */
	public SequenceData execute(OptionsHandle myOptions, SequenceData seqData) {

		boolean reversed = false;

		// 1. Create output SequenceData from given SequenceData
		SequenceData output = seqData.getNewSequenceData();

		// 2. Reverse new sequence and set output data to be reversed sequence
		if (output instanceof CharSequenceData) {
			((CharSequenceData) output)
					.setCharData(reverse(((CharSequenceData) output)
							.getCharData()));
			reversed = true;
		}

		else if (output instanceof DoubleSequenceData) {
			((DoubleSequenceData) output)
					.setDoubleData(reverse(((DoubleSequenceData) output)
							.getDoubleData()));
			reversed = true;
		}

		// reverse SequenceData as a sequence of objects if it isn't a
		// double or char sequence or if it wasn't reversed before
		else if (!reversed)
			output.setData(reverse(output.getData()));

		// 3. let output SequenceData know about this function
		output.addHistory(this);

		return output;

	}

	/**
	 * This functon reverses an array of Objects
	 * 
	 * @param data
	 *            Object[]
	 * @return Object[]
	 */
	private Object[] reverse(Object[] data) {
		Object temp;
		int last = data.length - 1;
		int first = 0;

		System.out.println("Reversing sequence ... ");

		while (first < last) {
			temp = data[first];
			data[first] = data[last];
			data[last] = temp;

			first++;
			last--;
		}

		return data;
	}

	/**
	 * This function reverses an array holding chars
	 * 
	 * @param data
	 *            Object[]
	 * @return Object[]
	 */
	private char[] reverse(char[] data) {
		char temp;
		int last = data.length - 1;
		int first = 0;

		System.out.println("Reversing sequence ... ");

		while (first < last) {
			temp = data[first];
			data[first] = data[last];
			data[last] = temp;

			first++;
			last--;
		}

		return data;
	}

	/**
	 * This function reverses an array holding doubles
	 * 
	 * @param data
	 *            Object[]
	 * @return Object[]
	 */
	private double[] reverse(double[] data) {
		double temp;
		int last = data.length - 1;
		int first = 0;

		System.out.println("Reversing sequence ... ");

		while (first < last) {
			temp = data[first];
			data[first] = data[last];
			data[last] = temp;

			first++;
			last--;
		}

		return data;
	}

	/* fixing toString() method */
	public String toString() {
		return "Reverse";
	}

}
