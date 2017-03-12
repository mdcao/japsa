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

package japsa.bio.misc.dnaPlatform.compModel;

import java.io.*;

import japsa.bio.misc.dnaPlatform.OptionsHandle;
import japsa.bio.misc.dnaPlatform.sequence.*;
import japsa.bio.misc.fuzzyLZ.FuzzyDriver;

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
public class FuzzyModel implements CompressionModel {

	String sequenceFile;
	String[] imageFiles;

	public FuzzyModel() {
	}

	/**
	 * Returns an OptionsHandle with options to run fuzzy model.
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		OptionsHandle myOptions = new OptionsHandle(this, 20);
		myOptions.addIntOption("maxIterations", 3,
				"Max number of iterations. (0 means until convergence)");
		myOptions.addStringOption("preFile", "",
				"Prepend 'preFile' to sequence.");
		myOptions
				.addStringOption(
						"seqModel",
						"markov(2)",
						"Base sequence model.  Use 'markov(n)' for a n-th order markov model, n=-1 for uniform");
		myOptions.addStringOption("dna", "atgc",
				"Alphabet used by the sequence.");

		myOptions.addStringOption("fwdMach", "3state",
				"Comma separated list of machines to use for forward matches.\n"
						+ "(Use an empty string '' for no forward machines.)\n"
						+ "Supported: 1state,3state");

		myOptions.addStringOption("revMach", "3state",
				"Comma separated list of machines to use for reverse matches.\n"
						+ "(Use an empty string '' for no reverse machines.)\n"
						+ "Supported: 1state,3state");

		myOptions.addBooleanOption("overwrite", true, "Overwrite msglen file.");
		myOptions.addIntOption("debug", 2,
				"Debug level (higher gives more verbose  )");
		myOptions.addIntOption("imageSize", 1024,
				"Maximum Image size in pixels");
		myOptions.addIntOption("imageFreq", 0,
				"Save an image every <n> seconds.  (0 - to disable)");
		myOptions.addIntOption("checkFreq", 0,
				"Save a checkpoint every <n> seconds.  (0 - to disable)");
		myOptions.addIntOption("statsFreq", 300,
				"Display some stats every <n> seconds.  (0 - to disable)");

		myOptions.addStringOption("msgFile", "",
				"Output file for encode length of each character.\n"
						+ "(The default is based on the input file name)");

		myOptions.addStringOption("outDir", "." + File.separatorChar,
				"Directory to save output files in.");

		// These options are for Matches_Sparse
		myOptions
				.addIntOption("hashSize", 10,
						"Window size to use for constructing hashtable (0 - for full N^2 algorithm)");
		myOptions.addIntOption("computeWin", 10,
				"Number of cells to activate on a hashtable hit");
		myOptions
				.addIntOption("cutML", 4,
						"When (cell_value - base_cell > cutML) then cell is killed. (in bits)");
		// myOptions.addBoolean("plotActive", false, "true: plot only active
		// cells, false: plot cell values");

		myOptions
				.addStringOption("paramFile", "",
						"Parameter file to read for various model parameters (see docs)");

		myOptions.addStringOption("resume", "",
				"Filename to resume from checkpoint");

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
	public SequenceData execute(OptionsHandle options, SequenceData data)
			throws IOException, RuntimeException {

		if (options == null) {
			System.out.println("options is null!");
			return null;
		}

		else {
			System.out.println("Owner of options is " + options.getOwner());
		}

		System.out.println(options.getStringValue("resume"));

		DoubleSequenceData outputData = new DoubleSequenceData(data);
		if (data instanceof DNASequenceData) {

			FuzzyDriver fuz = new FuzzyDriver();
			String infoContentFile = new String();

			// run fuzzyLZ
			if (data instanceof DNASequenceData) {
				infoContentFile = fuz.start(options, data.toString(),
						((DNASequenceData) data).getCharData());

				// add this model to the history of outputData
				outputData.addHistory(this);
				outputData.addHistory(options);

				// read double data obtained from Fuzzy
				outputData.readDataFromFile(infoContentFile);

				// set image files from fuzzy;
				imageFiles = fuz.getImageFileNames();
			}
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
		return imageFiles;
	}

	public String toString() {
		return "ARM";
	}

}
