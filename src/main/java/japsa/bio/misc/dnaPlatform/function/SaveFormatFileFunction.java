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
import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

import java.io.*;

/**
 * <p>
 * Title: Save one or more sequecce to a file
 * </p>
 * 
 * <p>
 * Description: This function saves sequences to a file.
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public class SaveFormatFileFunction implements Function {
	private File file;
	Sequence charSeq = null;
	JapsaAnnotation annoSeq = null;

	public SaveFormatFileFunction(File f) {
		file = f;
	}

	public void setCharSequence(CharSequenceData charData) {
		if (charData instanceof DNASequenceData) {
			charSeq = new Sequence(Alphabet.DNA4(), charData.getCharData(),
					charData.getSequenceName());
		}

	}

	public void setAnnotationSequence(AnnotationSequenceData annoData) {
		if (annoData != null)
			this.annoSeq = annoData.getAnnotation();
	}

	public Class[] getTypeSequenceData() {
		Class[] typeData = { CharSequenceData.class,
				AnnotationSequenceData.class };
		return typeData;
	}

	/**
	 * ReadFromFile does not have options therefore return null
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle() {
		OptionsHandle myOps = new OptionsHandle(this, 2);
		DoubleSequenceData temp = new DoubleSequenceData();
		temp.setDoubleData(new double[0]);
		myOps.addSequenceDataOption("Annotation", temp,
				"The Annotation to save with");
		return myOps;

	}

	/**
	 * OptionsHandle must be null as this function doesn't have options. seqData
	 * is sequence to read from a file.
	 * 
	 * @param myOptions
	 *            OptionsHandle
	 * @param seqData
	 *            SequenceData
	 * @return SequenceData
	 */
	public SequenceData execute(OptionsHandle myOptions, SequenceData seqData)
			throws IOException {
		if (seqData instanceof DNASequenceData)
			charSeq = new Sequence(Alphabet.DNA4(),
					((DNASequenceData) seqData).getCharData(),
					seqData.getSequenceName());

		// Get from option
		SequenceData annoData = myOptions.getSequenceDataValue("Annotation");
		if (annoData instanceof AnnotationSequenceData) {
			setAnnotationSequence((AnnotationSequenceData) annoData);
		}

		SequenceOutputStream out = new SequenceOutputStream(
				new FileOutputStream(file));
		JapsaAnnotation.write(charSeq, annoSeq, out);
		out.close();

		return null;
	}

	public String toString() {
		return "Save japsa.seq with annotation";
	}

}
