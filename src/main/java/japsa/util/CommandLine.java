/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
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
 * 3. Neither the names of the institutions nor the names of the contributors*
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

/*                           Revision History                                
 * 08/01/2012 - Minh Duc Cao: Revised                                        
 * 30/12/2012 - Minh Duc Cao: Add option to make a parameter compulsory
 ****************************************************************************/

package japsa.util;

import java.util.ArrayList;

import japsa.util.deploy.Deploy;

/**
 * An implementation of commandLine utilities. This class was written based
 * heavily on code from David Powell
 * 
 * @author Minh Duc Cao (rewrite from David Powell code)
 */

public class CommandLine {
	/**
	 * String describe the usage of the program (program -i input -o output f1 f2 ...)
	 */
	private String usage = "";
	/**
	 * A concise description of the program
	 */
	private String desc = "";
	//List of options
	ArrayList<Option> options;

	//Error message during parse
	String errors = null;

	//Keep the commandline
	String fullCmd = "";

	public CommandLine() {
		options = new ArrayList<Option>();
	}

	public CommandLine(String usageMsg) {
		this();
		usage = usageMsg;
	}

	public CommandLine(String usageMsg, String desc) {
		this(usageMsg);
		this.desc = desc;
	}

	public void setUsage(String usg){
		usage = usg;
	}

	public void setDesc(String desc){
		this.desc = desc;
	}

	public String errors(){
		return errors;
	}

	public String fullCmd(){
		return fullCmd;
	}

	private void addError(String errorMsg){
		if (errors == null)
			errors = errorMsg + "\n";
		else
			errors = errors + errorMsg + "\n";
	}

	private static String spaces(int num) {
		if (num <=1) return " ";
		char[] s = new char[num];
		for (int j = 0; j < num; j++)
			s[j] = ' ';
		return new String(s);
	}

	private static String indentLines(String s, int indent) {
		String[] lines = s.split("\n");
		StringBuffer res = new StringBuffer();
		for (int i = 0; i < lines.length; i++) {
			if (i > 0)
				res.append("\n" + spaces(indent));
			res.append(lines[i]);
		}
		return res.toString();
	}


	public String usage(){
		return usage;
	}

	/**
	 * Get the string describing all options of the program
	 * @return
	 */

	public String options(){
		int indent = 18;
		StringBuffer res = new StringBuffer();
		for (Option option:options) {
			res.append("  --" + option.optName);

			int len = option.optName.length();
			switch (option.optType) {
			case 'b':
				break;
			case 'i':
				res.append("=i");
				len += 2;
				break;
			case 'f':
				res.append("=d");
				len += 2;
				break;
			case 's':
				res.append("=s");
				len += 2;
				break;
			}

			res.append(spaces(indent - len - 4));			

			res.append(indentLines(option.optHelp, indent));
			res.append("\n" + spaces(indent) + (option.required?"(REQUIRED)":"(default='" + option.defaultValue + "')"));
			res.append("\n");
		}
		return res.toString();
	}

	/**
	 * Deprecated: use usageString instead
	 * @return
	 */
	@Deprecated
	public String usageMessage() {
		int indent = 18;
		StringBuffer res = new StringBuffer(usage + "\nOptions:\n");
		for (Option option:options) {
			res.append("  --" + option.optName);

			int len = option.optName.length();
			switch (option.optType) {
			case 'b':
				break;
			case 'i':
				res.append("=i");
				len += 2;
				break;
			case 'f':
				res.append("=d");
				len += 2;
				break;
			case 's':
				res.append("=s");
				len += 2;
				break;
			}

			res.append(spaces(indent - len - 4));			

			res.append(indentLines(option.optHelp, indent));
			res.append("\n" + spaces(indent) + (option.required?"(REQUIRED)":"(default='" + option.defaultValue + "')"));
			res.append("\n");
		}

		return res.toString();
	}

	private Option addOption(String opt, char type, Object def, String help, boolean req) {
		Option option = new Option(opt, type, def, help, req); 
		options.add(option);
		return option;
	}

	public Option addBoolean(String opt, boolean def, String help, boolean req) {
		return addOption(opt, 'b', new Boolean(def), help, req);
	}

	public Option addInt(String opt, int def, String help, boolean req) {
		return addOption(opt, 'i', new Integer(def), help, req);
	}

	public Option addDouble(String opt, double def, String help, boolean req) {
		return addOption(opt, 'f', new Double(def), help, req);
	}

	public Option addString(String opt, String def, String help, boolean req) {
		return addOption(opt, 's', def, help, req);
	}	

	public Option addBoolean(String opt, boolean def, String help) {
		return addOption(opt, 'b', new Boolean(def), help, false);
	}

	public Option addInt(String opt, int def, String help) {
		return addOption(opt, 'i', new Integer(def), help, false);
	}

	public Option addDouble(String opt, double def, String help) {
		return addOption(opt, 'f', new Double(def), help, false);
	}

	public Option addString(String opt, String def, String help) {
		return addOption(opt, 's', def, help, false);
	}

	public boolean optionSet(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return false;
		}
		return options.get(o).optionSet;
	}

	Object getVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return null;
		}
		return options.get(o).value;
	}

	public int getIntVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return 0;
		}

		if (options.get(o).optType != 'i') {
			System.err.println("ERROR: Option '" + opt
					+ "' is not an int option in getIntVal");
			return 0;
		}

		return ((Integer) options.get(o).value);
	}

	public double getDoubleVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return 0;
		}

		if (options.get(o).optType != 'f') {
			System.err.println("ERROR: Option '" + opt
					+ "' is not a double option in getDoubleVal");
			return 0;
		}

		return ((Double) options.get(o).value);
	}

	public String getStringVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return null;
		}

		if (options.get(o).optType != 's') {
			System.err.println("ERROR: Option '" + opt
					+ "' is not a string option in getStringVal");
			return null;
		}

		return (String) options.get(o).value;
	}

	public boolean getBooleanVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return false;
		}

		if (options.get(o).optType  != 'b') {
			System.err.println("ERROR: Option '" + opt
					+ "' is not a boolean option in getBooleanVal");
			return false;
		}

		return ((Boolean) options.get(o).value).booleanValue();
	}

	int isOption(String opt, int noDashOk) {
		if (opt.startsWith("--"))
			opt = opt.substring(2);
		else if (opt.startsWith("-"))
			opt = opt.substring(1);
		else if (noDashOk == 0)
			return -1;

		if (opt.indexOf("=") >= 0) {
			opt = opt.substring(0, opt.indexOf("="));
		}

		int match = -2;
		for (int i = 0; i < options.size(); i++) {
			if (opt.length() > options.get(i).optName.length())
				continue;

			String optStr = options.get(i).optName.substring(0, opt.length());
			if (opt.compareToIgnoreCase(optStr) == 0) {
				if (match >= 0) {
					addError("ERROR: Ambiguous option '" + opt
							+ "' could be '" + options.get(i).optName + "' or '" + options.get(match).optName
							+ "'");
					return match;
				}
				match = i;
			}
		}
		return match;
	}

	public String errorString(){
		return (errors == null?"":errors) + "\n" + "Usage: " + usage() + "\nOptions:\n" + options();		
	}
	public String usageString(){
		return desc + "\n\n" + "Usage: " + usage() + "\nOptions:\n" + options();
	}	

	public String[] stdParseLine(String[] args) {				
		/**********************************************************************/
		String[] ret = parseLine(args);
		//System.out.println(optionValues());
		if (isOption("help", 1) >=0 && getBooleanVal("help")){
			System.out.println(usageString());			
			System.exit(0);
		}

		if (errors != null) {
			System.out.println(errorString());			
			System.exit(-1);
		}	
		/**********************************************************************/
		return ret;
	}

	public String optionValues(){
		String ret = "";
		for (Option option:options){
			ret = ret + String.format("%20s = ",option.optName) + ((option.value == null)?"(null)":option.value) + "\n";			
		}		

		return ret;
	}
	public String[] parseLine(String[] args) {
		//Keep the original command line
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < args.length;i++){
			sb.append(' ');
			sb.append(args[i]);			
		}
		fullCmd = sb.toString();		


		int i = 0;
		int j = 0;
		String[] res = new String[args.length];

		// First setup defaults
		for (Option option:options)
			option.value = option.defaultValue;

		while (i < args.length) {
			int o = isOption(args[i], 0);
			if (o == -2) {
				addError("ERROR: Unknown option '"+args[i]+"'");
				return null;
			}

			if (o >= 0) {
				Option option = options.get(o);
				option.optionSet = true;

				try {
					switch (option.optType) {
					case 'b':
						if (args[i].indexOf("=") >= 0) {
							String s = args[i]
									.substring(args[i].indexOf("=") + 1);
							if (s.equalsIgnoreCase("true"))
								option.value = Boolean.TRUE;
							else if (s.equalsIgnoreCase("yes"))
								option.value = Boolean.TRUE;
							else if (s.equalsIgnoreCase("on"))
								option.value = Boolean.TRUE;
							else if (s.equalsIgnoreCase("1"))
								option.value = Boolean.TRUE;
							else if (s.equalsIgnoreCase("false"))
								option.value = Boolean.FALSE;
							else if (s.equalsIgnoreCase("no"))
								option.value = Boolean.FALSE;
							else if (s.equalsIgnoreCase("off"))
								option.value = Boolean.FALSE;
							else if (s.equalsIgnoreCase("0"))
								option.value = Boolean.FALSE;
							else {
								System.err
								.println("ERROR: Unknown boolean option parameter '"
										+ s + "'");
								return null;
							}
						}else
							option.value = Boolean.TRUE;
						break;
					case 'i':
						if (args[i].indexOf("=") >= 0)
							option.value = new Integer(args[i].substring(args[i].indexOf("=") + 1));
						else {
							option.value = new Integer(args[i + 1]);
							i++;
						}
						break;
					case 'f':
						if (args[i].indexOf("=") >= 0)
							option.value = new Double(args[i].substring(args[i].indexOf("=") + 1));
						else {
							option.value = new Double(args[i + 1]);
							i++;
						}
						break;
					case 's':
						if (args[i].indexOf("=") >= 0)
							option.value = args[i]
									.substring(args[i].indexOf("=") + 1);
						else {
							option.value = args[i + 1];
							i++;
						}
						break;
					}
				} catch (Exception e) {
					addError("ERROR: Bad option '" + args[i] + "'");
					return null;
				}
			} else {
				res[j++] = args[i];
			}

			i++;
		}

		//check if any required params not set

		boolean pass = true;
		for (Option option:options)
			if (option.required & (!option.optionSet)){
				addError("ERROR: The required param '"+option.optName+"' is not specified.");
				pass = false;
			}
		if (!pass)
			return null;

		String[] r = new String[j];
		for (int k = 0; k < j; k++)
			r[k] = res[k];
		return r;
	}

	///Some standard options	
	public void addStdInputFile(){
		addString("input", null, "Name of the input file, - for standard input", true);
	}

	public void addStdOutputFile(){
		addString("output", null, "Name of the output file, - for standard output", true);
	}

	public void addStdAlphabet(){
		addString("alphabet", "DNA", "Alphabet of the input file. Options: DNA (DNA=DNA16), DNA4\n(ACGT), DNA5(ACGTN), DNA16 and Protein");
	}

	public void addStdHelp(){
		addBoolean("help", false, "Display this usage and exit");
	}

	boolean galaxy = false;
	String cmd = null;	
	
	public void setGalaxy(String cmd){
		galaxy = true;
		this.cmd = cmd;
	}

	public boolean galaxy(){
		return galaxy;
	}
	
	

	public String generateGalaxyWrapper(){

		String  inputs = " <inputs>\n";
		String outputs = " <outputs>\n";

		String wrapper =
				"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" +		      
						"<tool id=\"" + cmd.replace("." , "_") +"\" name=\"" + cmd + "\" version=\"" + Deploy.VERSION + "\">\n" +
						" <description>" + desc + "</description>\n" +
						" <stdio>\n" +
						"  <exit_code range =\"1:\" />\n" +
						"  <exit_code range =\":-1\" />\n" +
						" </stdio>\n";

		wrapper = wrapper +
				" <command>\n"  + cmd;	

		for (Option option:options){
			GalaxySetting galaxySetting = option.galaxySetting;
			if (galaxySetting != null){
				wrapper += " --" + option.optName + " $" + option.optName;
				String params = galaxySetting.isOutput?
						"  <data name=\""+ option.optName + "\""
						:
							"  <param name=\""+ option.optName + "\"";
				if (galaxySetting.label != null){
					params +=  " label=\"" + galaxySetting.label + "\"";
				}else
					params +=  " label=\"" + option.optHelp + "\"";

				if (galaxySetting.format !=null){
					params += " format=\"" + galaxySetting.format + "\"";					
				}				

				if (galaxySetting.help != null){
					params +=  " help=\"" + galaxySetting.help + "\"";
				}

				if (!galaxySetting.isOutput){
					if (option.optType == 'b'){
						params += " type=\"boolean\" checked=\"" + option.defaultValue + "\" truevalue=\"true\" falsevalue=\"false\"";	
					}else if (option.optType == 'i'){
						params += " type=\"integer\" value=\"" + option.defaultValue +"\"";	
					}				
					else if (galaxySetting.type != null){
						params += " type=\""+ galaxySetting.type  + "\"";
					}
				}

				params += " />\n";


				if (galaxySetting.isOutput){					
					outputs = outputs + params;						 

				}else{
					inputs = inputs + params;
				}
			}
		}
		wrapper = wrapper + " </command>\n";
		wrapper = wrapper + inputs + " </inputs>\n" + outputs + " </outputs>\n";
		wrapper = wrapper + 
				"  <help>\n" +
				"**What it does**\n" +
				desc + "\n" +
				"  </help>\n";
		wrapper = wrapper + 
				"</tool>";

		/************************************************************************		



        <inputs>

                <param name="title_col" type="data_column" data_ref="input" multiple="True" numerical="False" label="Title column(s)" help="Multi-select list - hold the appropriate key while clicking to select multiple columns"/>
                <param name="seq_col" type="data_column" data_ref="input" numerical="False" label="Sequence column" />
        </inputs>
        <outputs>
                <data name="output" format="fasta"/>
        </outputs>
        <tests>
                <test>
                        <param name="input" value="solexa.tabular" />
                        <param name="title_col" value="1,2,3,4" />
                        <param name="seq_col" value="5" />
                        <output name="output" file="tabular_to_fasta_out1.fasta" />
                </test>
        </tests>
        <help>

		 **What it does**

Converts tab delimited data into FASTA formatted sequences.

-----------

		 **Example**

Suppose this is a sequence file produced by Illumina (Solexa) sequencer::

        5       300     902     419     GACTCATGATTTCTTACCTATTAGTGGTTGAACATC
        5       300     880     431     GTGATATGTATGTTGACGGCCATAAGGCTGCTTCTT

Selecting **c3** and **c4** as the **Title column(s)** and **c5** as the **Sequence column** will result in::

        &gt;902_419
        GACTCATGATTTCTTACCTATTAGTGGTTGAACATC
        &gt;880_431
        GTGATATGTATGTTGACGGCCATAAGGCTGCTTCTT

        </help>
</tool>

		/************************************************************************/


		return wrapper;
	}


	/**
	 * Representing option for galaxy
	 * 
	 * @author minhduc
	 *
	 */
	public static class GalaxySetting{
		//For use in galaxy wraper:		
		String  format = null;//bam, text, etc
		boolean isOutput = false;
		String  type = null;
		String  label = null;
		String  help = null;

		public GalaxySetting(String type, String format,  boolean isOutput){
			this.format = format;
			this.type = type;
			this.isOutput = isOutput;
		}


		public GalaxySetting(String type, String format){
			this (type,format,false);
		}

		/**
		 * @return the galaxyFormat
		 */
		public String getGalaxyFormat() {
			return format;
		}
		/**
		 * @param galaxyFormat the galaxyFormat to set
		 */
		public void setGalaxyFormat(String format) {
			this.format = format;
		}	

		public void setLabel(String label){
			this.label = label;
		}

		public String getLabel(){
			return label;
		}


	}
	/**
	 * Represent an option from the command line
	 * @author minhduc
	 *
	 */
	public static class Option{		
		String optName; //The name of the option
		char optType;   //The type of the option
		Object defaultValue;
		Object value = null;
		String optHelp;
		boolean required;
		boolean optionSet = false;

		private GalaxySetting galaxySetting = null;

		public Option(String opt, char type, Object def, String help, boolean req) {
			optName = opt;			
			optType = type;
			defaultValue = def;
			optHelp = help;			
			required = req;
		}
		public Option(String opt, char type, Object def, String help) {
			this(opt, type, def,help,true);
		}

		public GalaxySetting getGalaxySetting(){
			return galaxySetting;
		}

		public void setGalaxySetting(GalaxySetting galaxy){
			galaxySetting = galaxy;
		}
		//

	}	
}
