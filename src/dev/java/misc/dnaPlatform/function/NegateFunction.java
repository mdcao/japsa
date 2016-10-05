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

package misc.dnaPlatform.function;

import misc.dnaPlatform.OptionsHandle;
import misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: NegateFunction
 * </p>
 * 
 * <p>
 * Description: This class implements a function that negates all numbers in a
 * double sequence.
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class NegateFunction implements Function {
	public NegateFunction() {
	}

	@SuppressWarnings("rawtypes")
	public Class[] getTypeSequenceData() {
		Class[] typeData = { DoubleSequenceData.class };
		return typeData;
	}

	/**
	 * Returns an OptionsHandle object for the smoothing function which has an
	 * option to change window size
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		return null;
	}

	/**
	 * Creates a new DoubleSequenceData from input SequenceData, negates new
	 * SequenceData and adds this function to its history.
	 * 
	 * @param seqData
	 *            SequenceData
	 * @return SequenceData
	 */
	public SequenceData execute(OptionsHandle myOptions, SequenceData seqData)
			throws RuntimeException {

		if (!(seqData instanceof DoubleSequenceData)) {
			throw new RuntimeException(
					"Incorrect type of sequence, this function operates on numeric sequences");
		}

		// 1. Create output SequenceData from given SequenceData
		SequenceData output = seqData.getNewSequenceData();

		// 2. Smooth new sequence and set its data to be new smoothed data
		if (output instanceof DoubleSequenceData)
			((DoubleSequenceData) output)
					.setDoubleData(negate(((DoubleSequenceData) output)
							.getDoubleData()));

		// 3. let output SequenceData know about this function
		output.addHistory(this);

		return output;

	}

	/**
	 * This functon negates an array of doubles
	 * 
	 * @param data
	 *            Object[]
	 * @return Object[]
	 */
	private double[] negate(double[] data) throws RuntimeException {

		double[] newData = new double[data.length];

		for (int i = 0; i < data.length; i++) {
			newData[i] = data[i] * -1;
		}

		return newData;
	}

	/* fixing toString() method */
	public String toString() {
		return "Negate";
	}

}
