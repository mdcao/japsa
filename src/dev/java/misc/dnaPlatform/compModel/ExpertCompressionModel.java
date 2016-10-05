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

package misc.dnaPlatform.compModel;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.xm.ExpertModel;
//import japsa.xm.ExpertModel;

import java.io.*;

import misc.dnaPlatform.OptionsHandle;
import misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: Fuzzy Model
 * </p>
 * 
 * <p>
 * Description: Class to use FyzzyDriver in DNAPlatform
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class ExpertCompressionModel implements CompressionModel {
	public ExpertCompressionModel() {
	}

	/**
	 * Returns an OptionsHandle with options to run expert model.
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		OptionsHandle myOptions = new OptionsHandle(this, 20);
		myOptions.addIntOption("hashSize", 11,
				"Hash size, should not exceed 13 as limitation of Java");
		myOptions.addIntOption("maxExpert", 100, "Expert Limit");
		myOptions.addIntOption("context", 20, "Context Length");

		myOptions.addDoubleOption("listenThreshold", .1, "Listen Threshold");

		myOptions.addIntOption("chances", 4,
				"Number of chances before expert  to be removed");

		myOptions.addStringOption("dna", "atgc",
				"Alphabet used by the sequence.");

		// myOptions.addDoubleOption()
		return myOptions;
	}

	/**
	 * Returns a Class array with one element, DNA_SEQUENCE indicating this
	 * model can only be applied to DNA sequences.
	 * 
	 * @return Class[]
	 */
	@SuppressWarnings("rawtypes")
	public Class[] getTypeSequenceData() {
		Class[] types = { DNASequenceData.class };
		return types;
	}

	/**
	 * Function to execute fuzzyLZ model using FuzzyDriver class.
	 * 
	 * @param options
	 *            OptionsHandle
	 * @param data
	 *            SequenceData
	 * @return SequenceData
	 * @throws IOException
	 * @throws RuntimeException
	 */
	public SequenceData execute(OptionsHandle options, SequenceData data) {
		// throws IOException, RuntimeException {

		if (options == null) {
			System.out.println("options is null!");
			return null;
		}

		else {
			System.out.println("Owner of options is " + options.getOwner());
		}

		DoubleSequenceData outputData = new DoubleSequenceData(data);
		if (data instanceof DNASequenceData) {
			int hashSize = options.getIntValue("hashSize");
			int context = options.getIntValue("context");
			int maxExpert = options.getIntValue("maxExpert");
			double listenThreshold = options.getDoubleValue("listenThreshold");
			int chances = options.getIntValue("chances");

			ExpertModel expertModel = new ExpertModel(hashSize,
					Alphabet.DNA4(), context, maxExpert, listenThreshold,
					chances, false);

			expertModel.printParams();
			DNASequenceData dnaData = (DNASequenceData) data;

			Sequence dna = new Sequence(Alphabet.DNA4(), dnaData.getCharData(),
					dnaData.getSequenceName());

			Sequence[] dnaArray = new Sequence[1];
			dnaArray[0] = dna;

			double[] costs = expertModel.encode(dnaArray);

			outputData.addHistory(this);
			outputData.addHistory(options);
			outputData.setDoubleData(costs);
		}
		return outputData;
	}

	/**
	 * Returns an array of Strings with names for all other files created by
	 * fuzzy.
	 * 
	 * @return String[]
	 */
	public String[] getOtherGraphs() {
		// other graphs created with fuzzyLZ are:
		return null;
	}

	public String toString() {
		return "Expert Model";
	}

}
