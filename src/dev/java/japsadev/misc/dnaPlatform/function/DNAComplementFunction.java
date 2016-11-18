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

package japsadev.misc.dnaPlatform.function;

import japsadev.misc.dnaPlatform.OptionsHandle;
import japsadev.misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: DNAComplementFunction
 * </p>
 * 
 * <p>
 * Description: This function transforms a DNA sequence stored in a
 * DNASequenceData object and creates a new DNASequenceData object containing
 * the complementary DNA sequence, in which A is matched with T, G is matched
 * with C and vice versa.
 * </p>
 * 
 * <p>
 * Copyright: Copyright (c) 2005
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class DNAComplementFunction implements Function {
	public DNAComplementFunction() {
	}

	@SuppressWarnings("rawtypes")
	public Class[] getTypeSequenceData() {
		Class[] typeData = { DNASequenceData.class };
		return typeData;
	}

	/**
	 * Returns null as the DNAComplement function has no options
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		return null;
	}

	/**
	 * Creates a new SequenceData from input SequenceData, smooths new
	 * SequenceData and adds this function to its history.
	 * 
	 * @param seqData
	 *            SequenceData
	 * @return SequenceData
	 */
	public SequenceData execute(OptionsHandle myOptions, SequenceData seqData)
			throws RuntimeException {

		if (!(seqData instanceof DNASequenceData)) {
			throw new RuntimeException(
					"Incorrect type of sequence, this function operates on DNA sequences");
		}

		// 1. Create output SequenceData from given SequenceData
		SequenceData output = seqData.getNewSequenceData();

		// 2. Smooth new sequence and set its data to be new smoothed data
		if (output instanceof DNASequenceData)
			((DNASequenceData) output)
					.setCharData(complement(((DNASequenceData) output)
							.getCharData()));

		// 3. let output SequenceData know about this function
		output.addHistory(this);

		return output;

	}

	/**
	 * This functon smooths array of doubles. It takes a window size as
	 * parameter, calculates the average of the information content of elements
	 * in array in window size and gives this value to the last element in the
	 * array.
	 * 
	 * @param data
	 *            Object[]
	 * @return Object[]
	 */
	private char[] complement(char[] data) throws RuntimeException {

		char[] newData = new char[data.length];

		for (int i = 0; i < data.length; i++) {
			switch (data[i]) {
			case 'a':
				newData[i] = 't';
				break;
			case 't':
				newData[i] = 'a';
				break;
			case 'g':
				newData[i] = 'c';
				break;
			case 'c':
				newData[i] = 'g';
				break;
			case 'A':
				newData[i] = 'T';
				break;
			case 'T':
				newData[i] = 'A';
				break;
			case 'G':
				newData[i] = 'C';
				break;
			case 'C':
				newData[i] = 'G';
				break;
			default:
				newData[i] = data[i];
				break;
			}
		}

		return newData;
	}

	/* fixing toString() method */
	public String toString() {
		return "DNA Complement";
	}

}
