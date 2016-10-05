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

import java.util.*;
import java.io.*;

/**
 * <p>
 * Title: SequenceData
 * </p>
 * 
 * <p>
 * Description: This is an abstract class to represent all types of sequences
 * within the DNAPlatform. Sequences hold actual sequence data and a vector
 * containing references to objects used to create them.
 * </p>
 * 
 * <p>
 * Copyright: Copyright (c) 2005
 * </p>
 * 
 * Walked through my Minh Duc Cao 14/11/2007
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public abstract class SequenceData implements Cloneable {

	// Vector to hold references to objects used to create current sequence data
	// Objects used to create data include sequenceFiles, modelHandles and
	// Functions
	protected Vector<Object> constructor;
	protected String sequenceName;
	protected SequenceData myParent;

	/**
	 * Creates an instance of SequenceData.
	 * 
	 */
	public SequenceData() {
		constructor = new Vector<Object>();
		myParent = null;
	}

	/**
	 * Creates an instance of SequenceData from parent SequenceData, parent
	 * history is added to new instance.
	 * 
	 * @param parent
	 *            SequenceData
	 */

	public SequenceData(SequenceData parent) {
		// constructor = new Vector();
		constructor = parent.getHistory();
		myParent = parent;
	}

	/**
	 * method returns the sequenceName of the SequenceData
	 * 
	 * @return String sequence name
	 * 
	 */
	public String getSequenceName() {
		return sequenceName;
	}

	/**
	 * method sets the sequenceName
	 * 
	 * @param String
	 *            sequence name
	 * 
	 */
	public void setSequenceName(String name) {
		sequenceName = name;
	}

	public SequenceData getParentSequence() {
		return myParent;
	}

	/*
	 * The name of sequence file is a String in Vector constructor
	 * 
	 * @return String
	 */
	public File getSequenceFile() {
		Iterator i = constructor.iterator();
		while (i.hasNext()) {
			Object o = i.next();
			if (o instanceof File) {
				return (File) o;
			}
		}
		return null;
	}

	/**
	 * This method adds a new object to the history of sequence
	 * 
	 * @param o
	 *            object
	 */
	public void addHistory(Object o) {
		constructor.add(o);
	}

	/**
	 * This method returns a vector holding references to objects used to create
	 * the sequence.
	 * 
	 * @return Vector
	 */
	public Vector<Object> getHistory() {
		Vector<Object> history = new Vector<Object>();
		history.add(this);
		history.add(constructor);
		return history;
	}

	/**
	 * Returns clone of current sequence data object.
	 * 
	 * @return SequenceData
	 */
	// Do i need to clone???
	public SequenceData getNewSequenceData() {
		SequenceData newSequenceData = null;
		try {
			newSequenceData = (SequenceData) this.clone();
			newSequenceData.constructor = this.getHistory();
		} catch (Exception ex) {
			System.err.println(ex);
			ex.printStackTrace();
		}
		return newSequenceData;
	}

	/**
	 * Function to return data stored in SequenceData object.
	 * 
	 * @return Object[]
	 */
	public abstract Object[] getData();

	/**
	 * Function to set the data stored in a SequenceData object.
	 * 
	 * @param newData
	 *            Object[]
	 */
	public abstract void setData(Object[] newData);

	/*
	 * need to modify toString() method so that every sequenceData object is
	 * unique
	 */
	public String toString() {
		if (sequenceName == null)
			return "Sequence";
		return sequenceName + "(Sequence)";
	}

	/*
	 * All sequence data objects should know how to read sequences from files
	 */
	public abstract int readDataFromFile(String file) throws IOException;

	public abstract boolean writeDataToFile(File file);

	public abstract String getProperty();

	public abstract int size();

}
