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

package misc.dnaPlatform.sequence;

import java.io.*;
import java.util.Vector;

/**
 * <p>
 * Title: CharSequenceData
 * </p>
 * 
 * <p>
 * Description: CharSequenceData is a DiscreteSequence data that stores actual
 * sequence data as characters
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class CharSequenceData extends SequenceData {
	/**
	 * data: the sequence of characters
	 */
	protected char data[];

	/**
	 * Initializes current CharSequenceData instance to hold 100,000 characters.
	 * Number of characters stored can be changed by reading a DNA file with
	 * method readDataFromDNAfile()
	 */
	public CharSequenceData() {
		super();
		data = new char[0];
	}

	/**
	 * Takes a parent SequenceData object to get information about how it was
	 * constructed.
	 * 
	 * @param parent
	 *            SequenceData
	 */

	public CharSequenceData(SequenceData parent) {
		super(parent);
		data = new char[0];
	}

	public String getProperty() {
		return "Char sequence of " + data.length + " charators";
	}

	/**
	 * Returns array of characters with data stored in object
	 * 
	 * @return char[]
	 */
	public char[] getCharData() {
		return (char[]) data.clone();
	}

	/**
	 * Sets data[] to be a copy of given charData array
	 * 
	 * @param charData
	 *            char[]
	 */
	public void setCharData(char[] charData) {
		data = charData;
	}

	/**
	 * method returns the string representation of data
	 * 
	 * @return String
	 * 
	 * 
	 */
	public String getStringData() {
		String s = "";
		for (int i = 0; i < data.length; i++)
			s = s + data[i];
		return s;
	}

	/**
	 * Function to return data stored in CharSequenceData as an array of
	 * Character objects
	 * 
	 * @return Character[]
	 */
	public Object[] getData() {
		Character[] tempData = new Character[data.length];
		for (int i = 0; i < data.length; i++)
			tempData[i] = new Character(data[i]);

		return tempData;
	}

	/**
	 * Function to set the data stored in a SequenceData object given an array
	 * of Character objects
	 * 
	 * @param newData
	 *            Character[]
	 */
	public void setData(Object[] newData) {
		data = new char[newData.length];
		for (int i = 0; i < newData.length; i++)
			data[i] = ((Character) newData[i]).charValue();
	}

	public int size() {
		return data.length;
	}

	/**
	 * Reads a file using a BufferedReader and stores characters found in file
	 * matching charRegex into char array
	 * 
	 * @param sequenceFile
	 *            name of the file contaning the sequence of characters to read
	 * @param charRegex
	 *            a regular expression describing the type of characters to be
	 *            read from the file
	 * @return the number of characters read from file
	 * 
	 *         public int readDataFromFile(String sequenceFile, String
	 *         charRegex) {
	 * 
	 *         int sequenceLength =
	 *         getSequenceLengthFromFile(sequenceFile,charRegex,""); data = new
	 *         char[sequenceLength];
	 * 
	 *         int totalLen = 0; if (sequenceFile == null) { return 0; }
	 * 
	 *         try { // Open a file of the given name. File file = new
	 *         File(sequenceFile);
	 * 
	 *         //store file in constructor vector addHistory(file);
	 * 
	 *         FileReader fr = new FileReader(file); BufferedReader bufRdr = new
	 *         BufferedReader(fr);
	 * 
	 *         String line = null; while ( (line = bufRdr.readLine()) != null) {
	 *         String[] words = line.split(""); for (int j = 0; j <
	 *         words.length; j++) { if (words[j].matches(charRegex)) { // if
	 *         there is an element in data equal to new element // point to that
	 *         same element in data array data[totalLen] = words[j].charAt(0);
	 *         totalLen++; //if totalLength = sequenceLength, stop looping if
	 *         (totalLen == sequenceLength) { break; } } } } //while
	 * 
	 *         bufRdr.close(); } catch (IOException ioex) {
	 *         System.err.println(ioex); }
	 * 
	 *         System.out.println("graph created!");
	 * 
	 *         sequenceLength = totalLen; return totalLen; }
	 */

	/**
	 * read a file specify which type the file is if the file is FASTA or
	 * GENBANK, use the appropriate from biojava otherwise use the normal
	 * method(read by BufferedReader)
	 * 
	 */
	public int readDataFromFile(String sequenceFile, String charRegex) {

		if (sequenceFile == null) {
			return 0;
		}
		Vector<Character> v = new Vector<Character>();

		try {
			// Open a file of the given name.
			File file = new File(sequenceFile);

			// store file in constructor vector
			addHistory(file);

			FileReader fr = new FileReader(file);
			BufferedReader bufRdr = new BufferedReader(fr);

			String line = null;
			while ((line = bufRdr.readLine()) != null) {
				String[] words = line.split("");
				for (int j = 0; j < words.length; j++) {
					if (words[j].matches(charRegex)) {
						v.add(words[j].charAt(0));
					}
				}
			}
			data = new char[v.size()];
			for (int i = 0; i < v.size(); i++)
				data[i] = v.get(i);

			bufRdr.close();
		} catch (IOException ioex) {
			ioex.printStackTrace();
		}
		System.out.println("graph created!");
		return data.length;

	}

	public boolean writeDataToFile(File file) {
		return true;
	}

	/**
	 * This function reads alphanumeric characters from file
	 * 
	 * @param sequenceFile
	 *            String
	 * @return int
	 * @throws Exception
	 */
	public int readDataFromFile(String sequenceFile) throws IOException {
		return readDataFromFile(sequenceFile, "[acgtACGT]");
	}

	/*
	 * need to modify toString() method so that every sequenceData object is
	 * unique
	 */
	public String toString() {
		if (sequenceName == null)
			return "Character sequence";
		return sequenceName + "(Character sequence)";
	}

}
