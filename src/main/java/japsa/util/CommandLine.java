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

/**
 * An implementation of commandLine utilities. This class was written based
 * heavily on code from David Powell
 * 
 * @author Minh Duc Cao (rewrite from David Powell code)
 */

public class CommandLine {
	
	private String usg = "";
	private String desc = "";
		
	private int maxOptions = 100;
	private String[] opts;
	private char[] types;
	private Object[] defaults;
	private String[] helps;
	private boolean[] option_set;
	private boolean[] requires;
	private int numOptions;

	Object[] values;	
	String errors = null;

	public CommandLine(int maxO) {
		maxOptions = maxO;
		opts = new String[maxOptions];
		types = new char[maxOptions];
		defaults = new Object[maxOptions];
		helps = new String[maxOptions];

		option_set = new boolean[maxOptions];
		requires = new boolean[maxOptions];
		values = new Object[maxOptions];
		numOptions = 0;
	}
	
	public CommandLine() {
		this(100);
	}

	public CommandLine(String usageMsg) {
		this(100);
		usg = usageMsg;
	}
	
	public CommandLine(String usageMsg, String desc) {
		this(usageMsg);
		this.desc = desc;
	}

	public Object[] getDefaults() {
		return defaults;
	}

	public void setDefaults(Object[] defaults) {
		this.defaults = defaults;
	}
	
	public String errors(){
		return errors;
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

	public String usage() {
		int indent = 18;
		StringBuffer res = new StringBuffer(usg + "\nOptions:  \n");
		for (int i = 0; i < numOptions; i++) {
			res.append("  --" + opts[i]);

			int len = opts[i].length();
			switch (types[i]) {
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

			res.append(indentLines(helps[i], indent));
			res.append("\n" + spaces(indent) + (requires[i]?"(REQUIRED)":"(default='" + defaults[i] + "')"));
			res.append("\n");
		}

		return res.toString();
	}

	private void addOption(String opt, char type, Object def, String help, boolean req) {
		if (numOptions >= maxOptions) {
			System.err.println("ERROR: Too many options to CommandLine class");
			System.exit(1);
		}
		opts[numOptions] = opt;
		types[numOptions] = type;
		defaults[numOptions] = def;
		helps[numOptions] = help;
		option_set[numOptions] = false;
		requires[numOptions] = req;
		
		numOptions++;
	}
	
	public void addBoolean(String opt, boolean def, String help, boolean req) {
		addOption(opt, 'b', new Boolean(def), help, req);
	}

	public void addInt(String opt, int def, String help, boolean req) {
		addOption(opt, 'i', new Integer(def), help, req);
	}

	public void addDouble(String opt, double def, String help, boolean req) {
		addOption(opt, 'f', new Double(def), help, req);
	}

	public void addString(String opt, String def, String help, boolean req) {
		addOption(opt, 's', def, help, req);
	}	

	public void addBoolean(String opt, boolean def, String help) {
		addOption(opt, 'b', new Boolean(def), help, false);
	}

	public void addInt(String opt, int def, String help) {
		addOption(opt, 'i', new Integer(def), help, false);
	}

	public void addDouble(String opt, double def, String help) {
		addOption(opt, 'f', new Double(def), help, false);
	}

	public void addString(String opt, String def, String help) {
		addOption(opt, 's', def, help, false);
	}

	public boolean optionSet(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return false;
		}
		return option_set[o];
	}

	Object getVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return null;
		}
		return values[o];
	}

	public int getIntVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return 0;
		}

		if (types[o] != 'i') {
			System.err.println("ERROR: Option '" + opt
					+ "' is not an int option in getIntVal");
			return 0;
		}

		return ((Integer) values[o]).intValue();
	}

	public double getDoubleVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return 0;
		}

		if (types[o] != 'f') {
			System.err.println("ERROR: Option '" + opt
					+ "' is not a double option in getDoubleVal");
			return 0;
		}

		return ((Double) values[o]).doubleValue();
	}

	public String getStringVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return null;
		}

		if (types[o] != 's') {
			System.err.println("ERROR: Option '" + opt
					+ "' is not a string option in getStringVal");
			return null;
		}

		return (String) values[o];
	}

	public boolean getBooleanVal(String opt) {
		int o = isOption(opt, 1);
		if (o < 0) {
			System.err.println("ERROR: Attempt to lookup non-defined option '" + opt
					+ "'");
			return false;
		}

		if (types[o] != 'b') {
			System.err.println("ERROR: Option '" + opt
					+ "' is not a boolean option in getBooleanVal");
			return false;
		}

		return ((Boolean) values[o]).booleanValue();
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
		for (int i = 0; i < numOptions; i++) {
			if (opt.length() > opts[i].length())
				continue;

			String optStr = opts[i].substring(0, opt.length());
			if (opt.compareToIgnoreCase(optStr) == 0) {
				if (match >= 0) {
					addError("ERROR: Ambiguous option '" + opt
							+ "' could be '" + opts[i] + "' or '" + opts[match]
							+ "'");
					return match;
				}
				match = i;
			}
		}
		return match;
	}
	
	
		
	
	public String[] stdParseLine(String[] args) {
		addStdHelp();		
		/**********************************************************************/
		String[] ret = parseLine(args);
		if (getBooleanVal("help")){
			System.out.println(desc + "\n" + usage());			
			System.exit(0);
		}
		
		if (errors != null) {
			System.err.println(errors + usage());
			System.exit(-1);
		}	
		/**********************************************************************/
		return ret;
	}

	public String[] parseLine(String[] args) {
		int i = 0;
		int j = 0;
		String[] res = new String[args.length];

		// First setup defaults
		for (int k = 0; k < numOptions; k++)
			values[k] = defaults[k];

		while (i < args.length) {
			int o = isOption(args[i], 0);
			if (o == -2) {
				 addError("ERROR: Unknown option '"+args[i]+"'");
				return null;
			}

			if (o >= 0) {
				option_set[o] = true;

				try {
					switch (types[o]) {
					case 'b':
						if (args[i].indexOf("=") >= 0) {
							String s = args[i]
									.substring(args[i].indexOf("=") + 1);
							if (s.equalsIgnoreCase("true"))
								values[o] = Boolean.TRUE;
							else if (s.equalsIgnoreCase("yes"))
								values[o] = Boolean.TRUE;
							else if (s.equalsIgnoreCase("on"))
								values[o] = Boolean.TRUE;
							else if (s.equalsIgnoreCase("1"))
								values[o] = Boolean.TRUE;
							else if (s.equalsIgnoreCase("false"))
								values[o] = Boolean.FALSE;
							else if (s.equalsIgnoreCase("no"))
								values[o] = Boolean.FALSE;
							else if (s.equalsIgnoreCase("off"))
								values[o] = Boolean.FALSE;
							else if (s.equalsIgnoreCase("0"))
								values[o] = Boolean.FALSE;
							else {
								System.err
										.println("ERROR: Unknown boolean option parameter '"
												+ s + "'");
								return null;
							}
						} else
							values[o] = Boolean.TRUE;
						break;
					case 'i':
						if (args[i].indexOf("=") >= 0)
							values[o] = new Integer(args[i].substring(args[i].indexOf("=") + 1));
						else {
							values[o] = new Integer(args[i + 1]);
							i++;
						}
						break;
					case 'f':
						if (args[i].indexOf("=") >= 0)
							values[o] = new Double(args[i].substring(args[i].indexOf("=") + 1));
						else {
							values[o] = new Double(args[i + 1]);
							i++;
						}
						break;
					case 's':
						if (args[i].indexOf("=") >= 0)
							values[o] = args[i]
									.substring(args[i].indexOf("=") + 1);
						else {
							values[o] = args[i + 1];
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
		for (int x = 0 ; x < numOptions; x++)
			if (requires[x] & (!option_set[x])){
				addError("ERROR: The required param '"+opts[x]+"' is not specified.");
				pass = false;
			}
		if (!pass)
			return null;
		
			

		String[] r = new String[j];
		for (int k = 0; k < j; k++)
			r[k] = res[k];
		return r;
	}

	public void printOptions() {
		System.out
				.println("======================================================");
		for (int i = 0; i < numOptions; i++) {
			System.out.printf("%15s = " + values[i] + "\n", opts[i]);
		}
		System.out
				.println("======================================================");
	}

	///Some standard options	
	public void addStdInputFile(){
		addString("input", null, "Name of the input file, - for standard input", true);
	}
	
	public void addStdOutputFile(){
		addString("output", null, "Name of the output file, - for standard output", true);
	}
	
	public void addStdAlphabet(){
		addString("dna", "DNA", "Alphabet of the input file. Options: DNA (extended DNA), DNA4\n(ACGT), DNA5(ACGTN), DNA16 and Protein");
	}
	
	public void addStdHelp(){
		addBoolean("help", false, "Display this usage and exit");
	}
}
