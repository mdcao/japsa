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

import japsa.bio.misc.dnaPlatform.OptionsHandle;
import japsa.bio.misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: CutFunction
 * </p>
 * 
 * <p>
 * Description: This class is used to cut SequenceData given a starting and
 * ending index in a SequenceData object
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public class SelectFunction implements Function {
	public SelectFunction() {
	}

	public Class[] getTypeSequenceData() {
		Class[] typeData = { SequenceData.class };
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

		myOps.addIntOption("Start position", 0,
				"Start position to cut (positions range from 1 to sequence length).");
		myOps.addIntOption("Last position", 0,
				"Last position to cut (positions range from 1 to sequence length).");

		return myOps;
	}

	/**
	 * Creates a new SequenceData from input SequenceData, cuts new SequenceData
	 * and adds itself to history of new SequenceData
	 * 
	 * @param seqData
	 *            SequenceData
	 * @return SequenceData
	 */
	public SequenceData execute(OptionsHandle myOptions, SequenceData seqData)
			throws RuntimeException {

		int firstIndex = myOptions.getIntValue("Start position") - 1;
		int lastIndex = myOptions.getIntValue("Last position") - 1;

		// 1. Create output SequenceData from given SequenceData
		SequenceData output = seqData.getNewSequenceData();

		// 2. Cut new sequence and set output data to be new sequence
		output.setData(cut(firstIndex, lastIndex, output.getData()));

		// 3. let output SequenceData know about this function
		output.addHistory(myOptions);

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
	private Object[] cut(int firstIndex, int lastIndex, Object[] data)
			throws RuntimeException {

		if (firstIndex < 0) {
			firstIndex = 0;
		}

		else if (firstIndex > data.length) {
			firstIndex = data.length;
		}

		if (lastIndex < 0) {
			lastIndex = 0;
		}

		else if (lastIndex > data.length) {
			lastIndex = data.length;
		}

		if (firstIndex > lastIndex) {
			throw new RuntimeException(
					"First position must be smaller than last position");
		}

		System.out.println("Selecting sequence ... [" + (firstIndex + 1) + "-"
				+ (lastIndex + 1) + "]");

		/*
		 * The new array ranges from firstIndex to lastIndex, including
		 * lastIndex
		 */
		Object[] temp = new Object[lastIndex - firstIndex + 1];

		for (int i = firstIndex, j = 0; i <= lastIndex; i++, j++)
			temp[j] = data[i];

		return temp;
	}

	/* fixing toString() method */
	public String toString() {
		return "Select";
	}

}
