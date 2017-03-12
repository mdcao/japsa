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

import java.io.*;

import japsa.bio.misc.common.NumericalSequence;

/**
 * <p>
 * Title: DoubleSequenceData
 * </p>
 * 
 * <p>
 * Description: This class holds values assigned to bases throughout a sequence.
 * Values are represented as doubles and are stored in an array of fixed length.
 * </p>
 * 
 * <p>
 * Copyright: Copyright (c) 2005
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class DoubleSequenceData extends SequenceData {
	private double data[];

	/**
	 * Creates a DoubleSequenceData to hols a maximum of given doubles.
	 * 
	 * @param sequenceLen
	 *            int
	 */
	public DoubleSequenceData() {
		super();
		data = new double[0];
	}

	/**
	 * Creates a DoubleSequenceData object from parent SequenceData.
	 * DoubleSequenceData is set to hold the same number of objects as parent.
	 * 
	 * @param parent
	 *            SequenceData
	 */

	public DoubleSequenceData(SequenceData parent) {
		super(parent);
		data = new double[0];
	}

	public String getProperty() {
		return "Numerical sequence of " + data.length + " numbers";
	}

	/**
	 * Returns an array with data stored in DoubleSequenceData object.
	 * 
	 * @return double[]
	 */
	public double[] getDoubleData() {
		return (double[]) data.clone();
	}

	/**
	 * Sets data[] to be a copy of given doubleData array
	 * 
	 * @param doubleData
	 *            double[]
	 */
	public void setDoubleData(double[] doubleData) {
		data = doubleData;
	}

	/**
	 * Function to return data stored in DoubleSequenceData as an array of
	 * Double objects
	 * 
	 * @return Double[]
	 */
	public Object[] getData() {
		Double[] tempData = new Double[data.length];
		for (int i = 0; i < data.length; i++)
			tempData[i] = new Double(data[i]);

		return tempData;
	}

	/**
	 * Function to set the data stored in a DoubleSequenceData object given an
	 * array of Double objects
	 * 
	 * @param newData
	 *            Double[]
	 */
	public void setData(Object[] newData) {

		data = new double[newData.length];
		for (int i = 0; i < newData.length; i++)
			data[i] = ((Double) newData[i]).doubleValue();
	}

	/**
	 * readDataFromFile: This function reads output in a file from a compression
	 * model on a sequence. Points are stored in an array of fixed length,
	 * infoContent, to access elements of the array an integer is used.
	 */
	public int readDataFromFile(String filename) {

		try {
			// Open a file of the given name.
			File file = new File(filename);
			addHistory(file);
			data = NumericalSequence.read(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("graph created!" + data.length);
		return data.length;
	}

	public boolean writeDataToFile(File file) {
		return (new NumericalSequence(data)).writeDataToFile(file);
		/*******************************************
		 * try{ PrintWriter pw = new PrintWriter(new FileOutputStream(file));
		 * pw.println("# Double data written by DNAGraphTool"); for (int i = 0;
		 * i < data.length; i++){ pw.println(i + "\t" + data[i]); }
		 * 
		 * pw.close(); }catch (Exception e){ e.printStackTrace(); }
		 * 
		 * return true; /
		 *******************************************/

	}

	public int size() {
		return data.length;
	}

	/*
	 * need to modify toString() method so that every sequenceData object is
	 * unique
	 */
	public String toString() {
		if (sequenceName == null)
			return "Numerical sequence";
		return sequenceName + "(Numerical sequence)";
	}

}
