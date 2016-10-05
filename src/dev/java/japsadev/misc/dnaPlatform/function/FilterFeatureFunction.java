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

import japsa.seq.JapsaFeature;
import japsadev.misc.dnaPlatform.OptionsHandle;
import japsadev.misc.dnaPlatform.sequence.*;

import java.util.Iterator;

/**
 * <p>
 * Title: FilterFeatureFunction
 * </p>
 * 
 * 
 * @author Minh Duc Cao
 * @version 1.0
 */
public class FilterFeatureFunction implements Function {
	public FilterFeatureFunction() {
	}

	@SuppressWarnings("rawtypes")
	public Class[] getTypeSequenceData() {
		Class[] typeData = { AnnotationSequenceData.class };
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

		myOps.addBooleanOption("Include", true, "Include or exclude");
		myOps.addStringOption("Features", "",
				"List of features separated by comma ','");

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

		if (!(seqData instanceof AnnotationSequenceData)) {
			throw new RuntimeException(
					"Incorrect type of sequence, this function operates on sequences of doubles");
		}

		AnnotationSequenceData annoData = (AnnotationSequenceData) seqData;
		Iterator<JapsaFeature> iter = annoData.iterator();

		String[] list = myOptions.getStringValue("Features").split(",");
		boolean include = myOptions.getBooleanValue("Include");

		for (int x = 0; x < list.length; x++) {
			list[x] = list[x].trim().toUpperCase();
		}

		// 1. Create output SequenceData from given SequenceData
		AnnotationSequenceData output = new AnnotationSequenceData(seqData);

		if (include) {// Include
			while (iter.hasNext()) {
				JapsaFeature feature = iter.next();
				for (int x = 0; x < list.length; x++) {
					if (list[x].equals(feature.getType().toUpperCase().trim())) {
						output.addFeature(feature.cloneFeature());
						break;
					}
				}

			}
		} else {// Exclude
			while (iter.hasNext()) {
				boolean added = true;
				JapsaFeature feature = iter.next();
				for (int x = 0; x < list.length; x++) {
					if (list[x].equals(feature.getType().toUpperCase().trim())) {
						added = false;
						break;
					}
					if (added)
						output.addFeature(feature.cloneFeature());
				}
			}

		}

		// 3. let output SequenceData know about this function
		output.addHistory(myOptions);
		return output;

	}

	/* fixing toString() method */
	public String toString() {
		return "Feature Filter";
	}

}
