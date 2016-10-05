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

import java.io.*;

import misc.dnaPlatform.OptionsHandle;
import misc.dnaPlatform.sequence.*;

//import org.biojavax.bio.seq.RichSequence.IOTools;
/**
 * 
 * @author hoangnguyen
 */
@SuppressWarnings("rawtypes")
public class SaveFileFunction implements Function {
	// public static int NUM_OF_CHARS_IN_LINE=60;
	private File file;

	/**
	 * Creates a new instance of saveFileFunction
	 */
	// non-argument constructor
	public SaveFileFunction() {
	}

	// constructor with argument
	public SaveFileFunction(File aFile) {
		file = aFile;
	}

	// override the method form super class
	public Class[] getTypeSequenceData() {
		Class[] typeData = { SequenceData.class };
		return typeData;
	}

	/**
	 * ReadFromFile does not have options therefore return null
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		return null;
	}

	/**
	 * OptionsHandle must be null as this function doesn't have options. seqData
	 * is sequence to read from a file. precondition: SequenceData is a
	 * CharSequenceData
	 * 
	 * @param myOptions
	 *            OptionsHandle
	 * @param seqData
	 *            SequenceData
	 * @return SequenceData
	 */
	public SequenceData execute(OptionsHandle myOptions, SequenceData seqData) {
		seqData.writeDataToFile(file);

		return null;
	}

	/**
	 * method returns the string representation of this function
	 * 
	 * @return String
	 */

	public String toString() {
		return "Save";
	}

	/**
	 * method sets the file object
	 * 
	 * @param fi
	 *            a File
	 * 
	 */
	public void setFile(File fi) {
		file = fi;
	}

}
