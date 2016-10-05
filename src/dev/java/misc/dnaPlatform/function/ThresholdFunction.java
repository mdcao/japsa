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

import japsa.seq.JapsaFeature;
import misc.dnaPlatform.OptionsHandle;
import misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: ThresholdFunction
 * </p>
 * 
 * <p>
 * Description: This function takes a DoubleSequenceData object and a window
 * size and slides a window over data calculating average and setting this
 * average to element in the middle of the window
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public class ThresholdFunction implements Function {
	public ThresholdFunction() {
	}

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
		OptionsHandle myOps = new OptionsHandle(this, 2);

		myOps.addDoubleOption("Threshold", 0.0, "Sliding window size");
		myOps.addBooleanOption("Up", true, "From threshold up?");

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

		if (!(seqData instanceof DoubleSequenceData)) {
			throw new RuntimeException(
					"Incorrect type of sequence, this function operates on sequences of doubles");
		}

		DoubleSequenceData doubleData = (DoubleSequenceData) seqData;
		double data[] = doubleData.getDoubleData();

		double thres = myOptions.getDoubleValue("Threshold");
		boolean upSide = myOptions.getBooleanValue("Up");

		// 1. Create output SequenceData from given SequenceData
		AnnotationSequenceData output = new AnnotationSequenceData(seqData);

		boolean in = false;
		int start = 0;
		double sum = 0.0;
		// Go thro the numerical sequence
		for (int i = 0; i < data.length; i++) {
			if (data[i] < thres == upSide && in) {// Outside
				in = false;
				JapsaFeature f = new JapsaFeature(start, i + 1 - start);
				f.setType("Threshold");
				f.addDesc("Length: " + f.getLength());
				f.addDesc("Significant:" + sum);
				f.setID("S" + start);

				output.addFeature(f);

			} else if (data[i] < thres != upSide) {
				if (!in) {
					// Start to be in the range
					in = true;
					start = i;
					sum = 0.0;// data[i];
				}
				sum += upSide ? (data[i] - thres) : (thres - data[i]);

			}
		}

		if (in) {
			JapsaFeature f = new JapsaFeature(start, data.length - start);
			f.setType("Threshold");
			f.addDesc("Length: " + f.getLength());
			f.addDesc("Significant:" + sum);
			f.setID("S" + start);

			output.addFeature(f);

		}

		System.out.println(output.size());

		// 3. let output SequenceData know about this function
		output.addHistory(myOptions);

		return output;

	}

	/* fixing toString() method */
	public String toString() {
		return "Thresholding";
	}

}
