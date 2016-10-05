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

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceOutputStream;

import java.io.*;
import java.util.Iterator;

import misc.dnaPlatform.function.ReadFormatFileFunction;

/**
 * <p>
 * Title: FeatureSequenceData
 * </p>
 * 
 * <p>
 * Description: FeatureSequenceData is a DiscreteSequence data that stores
 * actual features
 * </p>
 * The features are in increasing order of starting points.
 * 
 * @author Minh Duc
 * @version 1.0
 */
public class AnnotationSequenceData extends SequenceData {
	JapsaAnnotation annotation;

	/**
	 * Initializes current CharSequenceData instance to hold 100,000 characters.
	 * Number of characters stored can be changed by reading a DNA file with
	 * method readDataFromDNAfile()
	 */
	public AnnotationSequenceData() {
		super();
		annotation = new JapsaAnnotation();
	}

	public AnnotationSequenceData(JapsaAnnotation x) {
		super();
		annotation = x;
	}

	/**
	 * Takes a parent SequenceData object to get information about how it was
	 * constructed.
	 * 
	 * @param parent
	 *            SequenceData
	 */

	public AnnotationSequenceData(SequenceData parent) {
		super(parent);
		annotation = new JapsaAnnotation();
	}

	/**
	 * Function to return data stored in FeatureSequenceData as an array of
	 * Character objects
	 * 
	 * @return Character[]
	 */
	public Object[] getData() {
		return annotation.getFeatureList().toArray();
	}

	public int size() {
		return annotation.numFeatures();
	}

	/**
	 * Function to set the data stored in a FeatureSequenceData object given an
	 * array of Features objects
	 * 
	 * @param newData
	 *            Character[]
	 */
	public void setData(Object[] newData) {
		annotation = new JapsaAnnotation();

		for (int i = 0; i < newData.length; i++)
			annotation.add((JapsaFeature) newData[i]);
	}

	public Iterator<JapsaFeature> iterator() {
		return annotation.iterator();
	}

	public boolean writeDataToFile(File file) {
		try {
			SequenceOutputStream ps = new SequenceOutputStream(
					new FileOutputStream(file));
			annotation.writeAnnotation(ps);
			ps.close();

		} catch (Exception e) {
		}
		return true;
	}

	public String getProperty() {
		return "Annotation containing " + size() + " features\n"
				+ annotation.getDescription();
	}

	public JapsaFeature getFeature(int idx) {
		return annotation.getFeature(idx);
	}

	public void addFeature(JapsaFeature f) {
		annotation.add(f);
	}

	public void addDescription(String desc) {
		annotation.addDescription(desc);
	}

	/**
	 * This function reads alphanumeric characters from file
	 * 
	 * @param sequenceFile
	 *            String
	 * @return int
	 */
	public int readDataFromFile(String sequenceFile) throws IOException {

		annotation = null;
		Iterator<SequenceData> iterSeq = (new ReadFormatFileFunction(
				sequenceFile)).guessFormat();
		while (iterSeq.hasNext()) {
			SequenceData seq = iterSeq.next();
			if (seq instanceof AnnotationSequenceData) {
				annotation = ((AnnotationSequenceData) seq).annotation;
				return annotation.numFeatures();
			}
		}
		throw new RuntimeException("No annotation found ");

	}

	public JapsaAnnotation getAnnotation() {
		return annotation;
	}

	public void setAnnotation(JapsaAnnotation annotation) {
		this.annotation = annotation;
	}

	/*
	 * need to modify toString() method so that every sequenceData object is
	 * unique
	 */
	public String toString() {
		if (sequenceName == null)
			return "Annotation sequence";
		return sequenceName + "(Annotation Sequence)";
	}

}
