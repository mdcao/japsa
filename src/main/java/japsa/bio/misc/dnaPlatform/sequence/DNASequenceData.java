/******************************************************************************
 * Copyright (C) 2006-2010 Minh Duc Cao                                        *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU General Public License as published by the Free  *
 * Software Foundation; either version 2 of the License, or (at your option)   *
 * any later version. This program is distributed in the hope that it will be  *
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General    *
 * Public License for more details.                                            *
 *                                                                             *
 * You should have received a copy of the GNU General Public License along with*
 * this program; if not, write to the Free Software  Foundation, Inc.,         *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                   *
 ******************************************************************************/

//This class is written by Julie Bernal and subsequently modified and maintained
//by Minh Duc Cao

package japsa.bio.misc.dnaPlatform.sequence;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * <p>
 * Title: DNASequenceData
 * </p>
 * 
 * <p>
 * Description: This class holds instances of DNA sequences
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class DNASequenceData extends CharSequenceData {

	public DNASequenceData() {
		super();

	}

	/**
	 * Creates an InfoContentSequenceData object from parent SequenceData.
	 * 
	 * @param parent
	 *            SequenceData
	 */

	public DNASequenceData(SequenceData parent) {
		super(parent);
	}

	public void readDataFromString(String str) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < str.length(); i++) {
			char c = Character.toLowerCase(str.charAt(i));
			if (c >= 'a' && c <= 'z')
				sb.append(c);
		}

		data = sb.toString().toCharArray();
	}

	public boolean writeDataToFile(File file) {
		try {
			Sequence aDNA = new Sequence(Alphabet.DNA4(), data, sequenceName);
			SequenceOutputStream out = new SequenceOutputStream(
					new FileOutputStream(file));

			aDNA.print(out);

			out.close();

			return true;
		} catch (Exception e) {
			e.printStackTrace();
		}

		return false;

	}

	/**
	 * Reads DNA sequence stored in file given as a parameter and stores that
	 * sequence in a character array.
	 * 
	 * @param sequenceFile
	 *            String
	 * @return int
	 */
	public int readDataFromFile(String sequenceFile) throws IOException {
		File file = new File(sequenceFile);
		addHistory(file);

		// set sequence of infoContent panel
		Sequence seq = SequenceReader.getReader(sequenceFile)
				.nextSequence(null);// (filename)IOTools.read(args[0]);
		if (seq == null) {
			System.err.println("Unable to read sequence file");
		}

		data = seq.charSequence();
		return data.length;
	}

	/**
	 * Function to return data stored in CharSequenceData as an array of
	 * Character objects, DNA sequences are composed of characters A, G, T, C
	 * 
	 * @return Character[]
	 */
	public Object[] getData() {

		Character[] tempData = new Character[data.length];
		for (int i = 0; i < data.length; i++) {
			if (data[i] == 'a' || data[i] == 'A')
				tempData[i] = 'a';
			else if (data[i] == 'g' || data[i] == 'G')
				tempData[i] = 'g';
			else if (data[i] == 't' || data[i] == 'T')
				tempData[i] = 't';
			else if (data[i] == 'c' || data[i] == 'C')
				tempData[i] = 'c';
			else
				tempData[i] = new Character(data[i]);
		}

		return tempData;
	}

	public String getProperty() {
		return "DNA sequence of " + data.length + " nucleotides";
	}

	/*
	 * need to modify toString() method so that every sequenceData object is
	 * unique
	 */
	public String toString() {
		if (sequenceName == null)
			return "DNA Sequence";
		return sequenceName + "(DNA)";
	}

}
