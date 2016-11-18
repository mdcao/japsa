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

package japsadev.misc.dnaPlatform;

import japsadev.misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: OptionsHandle
 * </p>
 * 
 * <p>
 * Description: A FunctionHandle holds options for functions. Options for
 * functions include integers, boolean, doubles, strings and SequenceData. An
 * OptionsHandle can also have other OptionsHandle as option
 * </p>
 * 
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class OptionsHandle {

	// this variable holds a reference to the owner of options
	private Object owner;

	private int MAXOPTIONS = 100;
	Option[] myOptions;
	private int numberOptions;

	private String newline;

	public OptionsHandle(Object optionsOwner) {
		owner = optionsOwner;
		myOptions = new Option[MAXOPTIONS];
		numberOptions = 0;

		newline = System.getProperty("line.separator");
	}

	public OptionsHandle(Object optionsOwner, int options) {
		owner = optionsOwner;

		MAXOPTIONS = options;
		myOptions = new Option[MAXOPTIONS];
		numberOptions = 0;

		newline = System.getProperty("line.separator");
	}

	/**
	 * This function returns the class of the object that created OptionsHandle
	 * 
	 * @return Class
	 */
	public Object getOwner() {
		return owner;
	}

	/* getNumberOfOptions: this function returns the number of options */
	public int getNumberOptions() {
		return numberOptions;
	}

	/* getOptionAt: returns the option name at given index */
	public String getOptionAt(int index) {
		if (index < 0 || index >= numberOptions)
			return null;
		else
			return myOptions[index].getOption();

	}

	/*
	 * optionSet: this function returns true if the value of a given option has
	 * been changed
	 */
	public boolean optionSet(String option) {
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)) {
				return myOptions[i].valueChanged();
			}
		}
		return false;
	}

	/* getHelp: given an option it returns a string with help */
	public String getHelp(String option) {
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)) {
				return myOptions[i].getHelp();
			}
		}
		return null;
	}

	/*** Functions to add options ***/

	/**
	 * Function to add an option of any type.
	 * 
	 * @param option
	 *            String
	 * @param defaultVal
	 *            Object
	 * @param help
	 *            String
	 */
	public void addOption(String option, Object defaultVal, String help) {
		if (numberOptions < MAXOPTIONS) {
			myOptions[numberOptions] = new Option(option, defaultVal, help);
			numberOptions++;
		}
	}

	/**
	 * Function to add option holding booleans
	 * 
	 * @param option
	 *            String
	 * @param defaultVal
	 *            boolean
	 * @param help
	 *            String
	 */
	public void addBooleanOption(String option, boolean defaultVal, String help) {
		if (numberOptions < MAXOPTIONS) {
			myOptions[numberOptions] = new BooleanOption(option, defaultVal,
					help);
			numberOptions++;
		}
	}

	/**
	 * Function to add option holding integers
	 * 
	 * @param option
	 *            String
	 * @param defaultVal
	 *            int
	 * @param help
	 *            String
	 */
	public void addIntOption(String option, int defaultVal, String help) {
		if (numberOptions < MAXOPTIONS) {
			myOptions[numberOptions] = new IntOption(option, defaultVal, help);
			numberOptions++;
		}
	}

	/**
	 * Function to add option holding doubles
	 * 
	 * @param option
	 *            String
	 * @param defaultVal
	 *            double
	 * @param help
	 *            String
	 */
	public void addDoubleOption(String option, double defaultVal, String help) {
		if (numberOptions < MAXOPTIONS) {
			myOptions[numberOptions] = new DoubleOption(option, defaultVal,
					help);
			numberOptions++;
		}
	}

	/**
	 * Function to add options to hold strings
	 * 
	 * @param option
	 *            String
	 * @param defaultVal
	 *            String
	 * @param help
	 *            String
	 */
	public void addStringOption(String option, String defaultVal, String help) {
		if (numberOptions < MAXOPTIONS) {
			myOptions[numberOptions] = new StringOption(option, defaultVal,
					help);
			numberOptions++;
		}
	}

	/**
	 * Function to add options holding SequenceData
	 * 
	 * @param option
	 *            String
	 * @param defaultVal
	 *            SequenceData
	 * @param help
	 *            String
	 */
	public void addSequenceDataOption(String option, SequenceData defaultVal,
			String help) {
		if (numberOptions < MAXOPTIONS) {
			myOptions[numberOptions] = new SequenceDataOption(option,
					defaultVal, help);
			numberOptions++;
		}
	}

	/**
	 * Function to add options holding OptionsHandle
	 * 
	 * @param option
	 *            String
	 * @param defaultVal
	 *            SequenceData
	 * @param help
	 *            String
	 */
	public void addOptionsHandleOption(String option, OptionsHandle defaultVal,
			String help) {
		if (numberOptions < MAXOPTIONS) {
			myOptions[numberOptions] = new OptionsHandleOption(option,
					defaultVal, help);
			numberOptions++;
		}
	}

	/*** function to get values of options ***/
	public Object getOptionValue(String option) {
		Option o = null;
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option))
				o = myOptions[i];
		}
		return o.getValue();
	}

	public boolean getBooleanValue(String option) {
		Option o = null;
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)
					&& myOptions[i].getValue() instanceof Boolean)
				o = myOptions[i];
		}
		return ((Boolean) o.getValue()).booleanValue();
	}

	public int getIntValue(String option) {
		Option o = null;
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)
					&& myOptions[i].getValue() instanceof Integer)
				o = myOptions[i];

		}
		return ((Integer) o.getValue()).intValue();
	}

	public double getDoubleValue(String option) {
		Option o = null;
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)
					&& myOptions[i].getValue() instanceof Double)
				o = myOptions[i];

		}
		return ((Double) o.getValue()).doubleValue();
	}

	public String getStringValue(String option) {
		Option o = null;
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)
					&& myOptions[i].getValue() instanceof String)
				o = myOptions[i];

		}
		return (String) o.getValue();
	}

	@SuppressWarnings({ "unused", "null" })
	public SequenceData getSequenceDataValue(String option) {
		Option o = null;
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)
					&& myOptions[i].getValue() instanceof SequenceData)
				o = myOptions[i];
			return (SequenceData) o.getValue();

		}
		return null;
	}

	public OptionsHandle getOptionsHandleValue(String option) {
		Option o = null;
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)
					&& myOptions[i].getValue() instanceof OptionsHandle)
				o = myOptions[i];

		}
		return (OptionsHandle) o.getValue();
	}

	/*** functions to set values of options ***/
	public void setOptionValue(String option, Object val) {
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)
					&& myOptions[i].getValue().getClass()
							.equals(val.getClass()))
				myOptions[i].setValue(val);
		}
	}

	public void setSequenceDataValue(String option, SequenceData val) {
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)
					&& myOptions[i] instanceof SequenceDataOption)
				myOptions[i].setValue(val);
		}
	}

	public void setOptionsHandleValue(String option, OptionsHandle val) {
		for (int i = 0; i < numberOptions; i++) {
			if (myOptions[i].getOption().equals(option)
					&& myOptions[i] instanceof OptionsHandleOption)
				myOptions[i].setValue(val);
		}
	}

	/**
	 * Function to create string representing an OptionsHandle object given as a
	 * parameter.
	 * 
	 * @param ops
	 *            OptionsHandle
	 * @return String
	 */
	private String optionsToString(OptionsHandle ops) {
		String optionString = "";
		for (int i = 0; i < ops.getNumberOptions(); i++) {

			String o = ops.getOptionAt(i);
			Object val = ops.getOptionValue(o);

			if (val instanceof Boolean)
				optionString += "\n" + o + " (boolean) = "
						+ ((Boolean) val).booleanValue();
			else if (val instanceof Integer)
				optionString += "\n" + o + " (integer) = "
						+ ((Integer) val).intValue();
			else if (val instanceof Double)
				optionString += "\n" + o + " (double) = "
						+ ((Double) val).doubleValue();
			else if (val instanceof String)
				optionString += "\n" + o + " (string) = " + val;
			// else if(val instanceof SequenceData)
			// optionString += "\n"+ o +" = " + val;
			else
				optionString += "\n" + o + " = " + val;

		}

		return optionString;
	}

	public String toString() {
		String options = owner + "{";
		options += optionsToString(this);
		options += "\n}";
		return options;
	}

	/**
	 * <p>
	 * Title: Option
	 * </p>
	 * 
	 * <p>
	 * Description: This class holds options and ModelHandle
	 * </p>
	 * 
	 */
	private class Option {

		protected String option; // holds name of option
		protected String help;
		protected Object defaultVal;
		protected Object value;

		public Option(String option, Object defaultVal, String help) {
			this.option = option;
			this.help = help;
			this.defaultVal = value = defaultVal;
		}

		public String getOption() {
			return option;
		}

		@SuppressWarnings("unused")
		public Object getDefaultVal() {
			return defaultVal;
		}

		public Object getValue() {
			return value;
		}

		public void setValue(Object val) {
			value = val;
		}

		// if value != defaultVal it means value has been changed
		public boolean valueChanged() {
			return value != defaultVal;
		}

		public String getHelp() {
			return help + newline + "(default = " + defaultVal + ")";
		}

	}

	/**
	 * <p>
	 * Title: BooleanOption
	 * </p>
	 * 
	 * <p>
	 * Description: This class holds boolean options
	 * </p>
	 * 
	 */
	private class BooleanOption extends Option {

		public BooleanOption(String option, boolean defaultVal, String help) {
			super(option, new Boolean(defaultVal), help);
		}

	}

	/**
	 * <p>
	 * Title: IntOption
	 * </p>
	 * 
	 * <p>
	 * Description: This class holds integer options
	 * </p>
	 * 
	 */
	private class IntOption extends Option {
		public IntOption(String option, int defaultVal, String help) {
			super(option, new Integer(defaultVal), help);
		}

	}

	/**
	 * <p>
	 * Title: DoubleOption
	 * </p>
	 * 
	 * <p>
	 * Description: This class holds double options
	 * </p>
	 * 
	 */
	private class DoubleOption extends Option {
		public DoubleOption(String option, double defaultVal, String help) {
			super(option, new Double(defaultVal), help);
		}
	}

	/**
	 * <p>
	 * Title: StringOption
	 * </p>
	 * 
	 * <p>
	 * Description: This class holds string options and ModelHandle
	 * </p>
	 * 
	 */
	private class StringOption extends Option {
		public StringOption(String option, String defaultVal, String help) {
			super(option, defaultVal, help);
		}

	}

	/**
	 * <p>
	 * Title: SequenceDataOption
	 * </p>
	 * 
	 * <p>
	 * Description: This class holds SequenceData options
	 * </p>
	 * 
	 */
	private class SequenceDataOption extends Option {

		public SequenceDataOption(String option, SequenceData defaultVal,
				String help) {
			super(option, defaultVal, help);
		}

	}

	/**
	 * <p>
	 * Title: OptionsHandleOption
	 * </p>
	 * 
	 * <p>
	 * Description: This class holds OptionHandle options
	 * </p>
	 * 
	 */
	private class OptionsHandleOption extends Option {

		public OptionsHandleOption(String option, OptionsHandle defaultVal,
				String help) {
			super(option, defaultVal, help);
		}

	}

	/*
	 * Main function for testing
	 */
	public static void main(String args[]) {

		String owner = "Master";
		OptionsHandle func = new OptionsHandle(owner);

		System.out.println("Testing function handle");
		System.out.println();

		// adding options into ModelHandle
		System.out.println("Creating options ... ");
		func.addIntOption("foo", 7, "Some integer");
		func.addBooleanOption("fom", false, "A bool");
		func.addDoubleOption("doub", 8.42, "A double");
		func.addStringOption("str", "cthulu", "A string");
		CharSequenceData s = new CharSequenceData();
		s.setSequenceName("char_1");
		func.addSequenceDataOption("aSequence", s,
				"This is a character sequence");

		OptionsHandle h = new OptionsHandle("Slave", 2);
		h.addBooleanOption("poor", true, "very unfortunate");
		h.addDoubleOption("money", 0.0, "All his money");

		func.addOptionsHandleOption("slave", h, "his worker");

		// printing option values
		System.out.println("These are the options: ");
		System.out.println("foo  = " + func.getOptionValue("foo"));
		System.out.println("fom  = " + func.getOptionValue("fom"));
		System.out.println("doub = " + func.getOptionValue("doub"));
		System.out.println("str  = " + func.getOptionValue("str"));

		// setting option values
		/*
		 * System.out.println("Setting option values ...");
		 * func.setOptionValue("foo",new Integer(4));
		 * func.setOptionValue("fom",new Boolean(true));
		 * func.setOptionValue("doub",new Double(4.5));
		 * func.setOptionValue("str","hello hello");
		 * 
		 * //printing new option values System.out.println("New values:");
		 * System.out.println("foo  = " + func.getOptionValue("foo"));
		 * System.out.println("fom  = " + func.getOptionValue("fom"));
		 * System.out.println("doub = " + func.getOptionValue("doub"));
		 * System.out.println("str  = " + func.getOptionValue("str"));
		 * 
		 * //attempting to set wrong value for options
		 * System.out.println("Setting options with wrong value type ...");
		 * func.setOptionValue("foo",new Boolean(false));
		 * func.setOptionValue("fom",new Integer(666));
		 * func.setOptionValue("doub","this is a string!");
		 * func.setOptionValue("str",new Double(0.0));
		 * 
		 * //printing new option values
		 * System.out.println("Values should be same as above:");
		 * System.out.println("foo  = " + func.getOptionValue("foo"));
		 * System.out.println("fom  = " + func.getOptionValue("fom"));
		 * System.out.println("doub = " + func.getOptionValue("doub"));
		 * System.out.println("str  = " + func.getOptionValue("str"));
		 * 
		 * 
		 * //getting types of options
		 * System.out.println("These are the types of the options");
		 * System.out.println("foo type  = " +
		 * func.getOptionValue("foo").getClass());
		 * System.out.println("fom type  = " +
		 * func.getOptionValue("fom").getClass());
		 * System.out.println("doub type = " +
		 * func.getOptionValue("doub").getClass());
		 * System.out.println("str type  = " +
		 * func.getOptionValue("str").getClass());
		 * 
		 * 
		 * System.out.println("Checking for types of options");
		 * if(func.getOptionValue("foo") instanceof Integer)
		 * System.out.println("foo is an integer!");
		 * if(func.getOptionValue("fom") instanceof Boolean)
		 * System.out.println("fom is a boolean!");
		 * if(func.getOptionValue("doub") instanceof Double)
		 * System.out.println("doub is a double!");
		 */

		// printing toString method of options handle dialog:
		System.out.println("toString(): ");
		System.out.println(func);

		System.out.println();
		System.out.println("Attempting to get all children out of string");
		String description = func.toString();

		System.out.println("\t" + description + "\n");

		// if toString() returns a string with curly braces
		// add all lines in between as child nodes
		if (description.indexOf('{') > -1) {
			System.out.println("There are children!");

			// add anything before brakets into tree
			String desStart = description.substring(0,
					description.indexOf('{') + 1);
			if (desStart.matches("\\S+"))
				System.out.println("\t" + desStart);

			// add children to a vector and call function to add children to
			// tree
			description = description.substring((description.indexOf('{') + 1),
					(description.lastIndexOf('}')));

			String[] children = description.split("[\\n]");

			for (int i = 0; i < children.length; i++)
				System.out.println("\t\t\t" + children[i]);

			// add anything after children
			System.out.println("The end bit is:");
			String end = description.substring(description.lastIndexOf('}'),
					description.length() - 1);
			if (end.matches("\\S+"))
				System.out.println("\t" + end);
		}

	}

}
