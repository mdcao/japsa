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
 * Title: AppendFunction
 * </p>
 * 
 * <p>
 * Description: This function takes an input sequence and another sequnce as a
 * parameter and appends these sequences if they are of the same SequenceData
 * type
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public class AppendFunction implements Function {
	public AppendFunction() {
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

		CharSequenceData temp = new CharSequenceData();
		temp.setCharData(new char[0]);

		myOps.addSequenceDataOption("Sequence", temp, "Sequence to append");

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

		SequenceData appendSeq = myOptions.getSequenceDataValue("Sequence");

		int appendIndex = seqData.getData().length;

		// 1. Create output SequenceData from given SequenceData
		SequenceData output = seqData.getNewSequenceData();

		// 2. Cut new sequence and set output data to be new sequence
		output.setData(append(output.getData(), appendSeq.getData()));

		// 3. let output SequenceData know about this function
		output.addHistory(this + " [" + appendIndex + "] ");
		output.addHistory(appendSeq.getHistory());

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
	private Object[] append(Object[] data1, Object[] data2)
			throws RuntimeException {

		if (!data1.getClass().equals(data2.getClass())) {
			throw new RuntimeException(
					"Sequences to append must be of same type");
		}

		System.out.println("Appending sequences ... ");

		Object[] temp = new Object[data1.length + data2.length];

		int index = 0;
		for (int i = 0; i < data1.length; i++, index++)
			temp[index] = data1[i];

		for (int i = 0; i < data2.length; i++, index++)
			temp[index] = data2[i];

		return temp;
	}

	/* fixing toString() method */
	public String toString() {
		return "Append";
	}

}
