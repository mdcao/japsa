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

package misc.dnaPlatform.function;

import misc.dnaPlatform.OptionsHandle;
import misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: SmoothingFunction
 * </p>
 * 
 * <p>
 * Description: This function takes a DoubleSequenceData object and a window
 * size and slides a window over data calculating average and setting this
 * average to element in the middle of the window. The first elements before the
 * first window is read are smoothed by setting their values as the average of
 * elements seen previously
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public class SmoothingFunction2 implements Function {
	public SmoothingFunction2() {
	}

	public Class[] getTypeSequenceData() {
		Class[] typeData = { DoubleSequenceData.class };
		return typeData;
	}

	/**
	 * Returns an OptionsHandle object for the smoothing2 function which has
	 * option to set window size.
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		OptionsHandle myOps = new OptionsHandle(this, 1);

		myOps.addIntOption("Window size", 1, "Sliding window size");
		return myOps;
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

		int winSize = myOptions.getIntValue("Window size");

		// 1. Create output SequenceData from given SequenceData
		SequenceData output = seqData.getNewSequenceData();

		// 2. Smooth new sequence and set its data to be new smoothed data
		if (output instanceof DoubleSequenceData)
			((DoubleSequenceData) output).setDoubleData(smooth(winSize,
					((DoubleSequenceData) output).getDoubleData()));

		// 3. let output SequenceData know about this function
		output.addHistory(myOptions);

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
	private double[] smooth(int winSize, double[] data) throws RuntimeException {

		double[] newData = new double[data.length];
		double tally = 0;

		if (winSize <= 0 && winSize >= data.length) {
			throw new RuntimeException("Bad window size " + winSize);
		}
		System.out.println("Smoothing double sequence data ... ");

		int index = 0;

		// coping first values from data into newData
		newData[0] = tally = data[0];
		for (index = 1; index < winSize; index++) {
			tally += data[index];
			newData[index] = tally / (index + 1);
		}

		// calculating rest of smoothed values and put into newData
		while (index < data.length) {
			tally = tally + data[index] - data[index - winSize];
			newData[index] = tally / winSize;
			index++;
		}

		return newData;
	}

	/* fixing toString() method */
	public String toString() {
		return "Smoothing II";
	}

}
