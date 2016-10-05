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
 * Title: Difference
 * </p>
 * 
 * <p>
 * Description: This function takes two DoubleSequenceData sequences and
 * calculates the difference between their values.
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public class DifferenceFunction implements Function {
	public DifferenceFunction() {
	}

	public Class[] getTypeSequenceData() {
		Class[] typeData = { DoubleSequenceData.class };
		return typeData;
	}

	/**
	 * Returns an OptionsHandle object for a particular function which has
	 * parameters to run function.
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		OptionsHandle myOps = new OptionsHandle(this, 2);

		DoubleSequenceData temp = new DoubleSequenceData();
		temp.setDoubleData(new double[0]);

		myOps.addSequenceDataOption("Sequence", temp,
				"Substract this sequence from input");

		return myOps;
	}

	/**
	 * Creates a new DoubleSequenceData containing the difference between
	 * sequence given as parameter and sequence in OptionsHandle
	 * 
	 * 
	 * @param seqData
	 *            SequenceData
	 * @return SequenceData
	 */
	public SequenceData execute(OptionsHandle myOptions, SequenceData seqData)
			throws RuntimeException {

		SequenceData diffData = myOptions.getSequenceDataValue("Sequence");

		if (!(seqData instanceof DoubleSequenceData && diffData instanceof DoubleSequenceData)) {
			throw new RuntimeException(
					"Incorrect type of sequence, this function operates on sequences of doubles");
		}

		// 1. Create output SequenceData from given SequenceData
		DoubleSequenceData output = new DoubleSequenceData(seqData);

		// 2. Subtract sequence given as a parameter from input sequence and set
		// it to output
		output.setDoubleData(difference(
				((DoubleSequenceData) seqData).getDoubleData(),
				((DoubleSequenceData) diffData).getDoubleData()));

		// 3. let output SequenceData know about this function
		output.addHistory(this);
		output.addHistory(diffData.getHistory());

		return output;
	}

	/**
	 * This function cuts a given object array from firstIndex to lastIndex and
	 * returns created array of objects. An array is only cut when the indexes
	 * are bigger than 0 and less than the lenght of array.
	 * 
	 * @param firstIndex
	 *            int
	 * @param lastIndex
	 *            int
	 * @param data
	 *            Object[]
	 * @return Object[]
	 */
	private double[] difference(double[] data1, double[] data2) {

		System.out.println("Calculating difference of sequences ... ");

		// set length for output to the shortest sequence
		double[] temp = new double[Math.max(data1.length, data2.length)];

		int index = 0;
		for (int i = 0; i < temp.length; i++, index++)
			temp[index] = data1[i] - data2[i];

		return temp;
	}

	/* fixing toString() method */
	public String toString() {
		return "Difference";
	}

}
