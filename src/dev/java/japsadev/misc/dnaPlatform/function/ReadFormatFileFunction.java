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

import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFileFormat;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsadev.misc.dnaPlatform.OptionsHandle;
import japsadev.misc.dnaPlatform.sequence.*;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * <p>
 * Title: ReadFileFunction
 * </p>
 * 
 * <p>
 * Description: This function reads a sequence from a file.
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class ReadFormatFileFunction implements Function {

	String filename = "";

	public ReadFormatFileFunction() {
	}

	public ReadFormatFileFunction(String f) {
		filename = f;
	}

	@SuppressWarnings("rawtypes")
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
	 * is sequence to read from a file.
	 * 
	 * @param myOptions
	 *            OptionsHandle
	 * @param seqData
	 *            SequenceData
	 * @return SequenceData
	 */
	public SequenceData execute(OptionsHandle myOptions, SequenceData seqData) {
		// never get called
		// seqData = nummy
		// if (filename != "" && filename != null) {
		// return guessFormat(filename);
		// }
		return null;
	}

	public void setFile(String sequenceFile) {
		filename = sequenceFile;
	}

	public String toString() {
		return "Read";
	}

	public Iterator<SequenceData> guessFormat() {
		ArrayList<SequenceData> seqList = new ArrayList<SequenceData>();
		try {
			SequenceReader reader = SequenceReader.getReader(filename);
			if (reader instanceof JapsaFileFormat) {
				JapsaFileFormat bff = (JapsaFileFormat) reader;
				JapsaAnnotation annoRead;
				while ((annoRead = bff.readAnnotation()) != null) {
					AnnotationSequenceData anno = new AnnotationSequenceData(
							annoRead);
					if (anno.size() > 0) {
						seqList.add(anno);
						anno.addHistory(filename);
					}

					Sequence seq = annoRead.getSequence();
					if (seq != null && seq.length() > 0
							&& seq.alphabet() == Alphabet.DNA4()) {
						DNASequenceData dna = new DNASequenceData();
						dna.setCharData(seq.charSequence());

						seqList.add(dna);
						dna.addHistory(filename);
					}
				}
			} else {
				throw new RuntimeException("Unknown format");
			}

			/**************************************************
			 * BufferedReader in = SequenceReader.openFile(filename); if (in ==
			 * null) return null;
			 * 
			 * in.mark(10);
			 * 
			 * char[] buf = new char[10]; in.read(buf, 0, 10); in.reset();
			 * String format = new String(buf);
			 * 
			 * int mode = -1; if (format.startsWith(JapseFileFormat.HEADER)) {
			 * mode = 0; } else if (format.startsWith("LOCUS")) {// Genbank mode
			 * = 1;// Biojava System.out.println("Read as Genbank"); // seqsIt =
			 * SeqIOTools.readGenbank(in); // seqsIt =
			 * IOTools.readGenbankDNA(in,null); throw new RuntimeException(
			 * "Genbank format has not been supported"); } else if
			 * (format.startsWith(">")) { // Fasta mode = 1;// Biojava
			 * System.out.println("Read as Fasta"); // seqsIt =
			 * IOTools.readFastaDNA(in,null); throw new RuntimeException(
			 * "fasta format has not been supported"); } else if
			 * (format.startsWith("ID")) { // Fasta mode = 1;// Biojava
			 * System.out.println("Read as EMBL"); // seqsIt =
			 * IOTools.readEMBLDNA(in,null); throw new
			 * RuntimeException("EMBL format has not been supported"); } else
			 * throw new Exception("Unknown file format");
			 * 
			 * {// mode == 0 = my format BioCompFileFormat fileFormat = new
			 * BioCompFileFormat(in); try { Iterator<Sequence> iterSeq =
			 * fileFormat.getSequenceIterator(); while (iterSeq.hasNext()) {
			 * Sequence japsa.seq = iterSeq.next(); if (japsa.seq.length() > 0
			 * && japsa.seq.alphabet() == Alphabet.DNA4()) { DNASequenceData dna
			 * = new DNASequenceData();
			 * dna.setCharData(japsa.seq.charSequence());
			 * 
			 * seqList.add(dna); dna.addHistory(filename); } }
			 * 
			 * Iterator<JapsaAnnotation> iter =
			 * fileFormat.getAnnotationIterator(); while (iter.hasNext()) {
			 * AnnotationSequenceData anno = new AnnotationSequenceData(
			 * iter.next()); if (anno.size() > 0) { seqList.add(anno);
			 * anno.addHistory(filename); } }
			 * 
			 * } catch (Exception e) { e.printStackTrace(); } } /
			 **************************************************/
			return seqList.iterator();
		} catch (Exception e) {
			// System.err.println("Error reading '"+filename+"' "+e);
			e.printStackTrace();
		}
		return null;

	}

}
