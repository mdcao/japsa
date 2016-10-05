/*
 * Copyright (c) David Powell <david@drp.id.au>
 *
 *
 * This file is part of FuzzyLZ
 *
 * FuzzyLZ is a program orginally intended for the compression of DNA sequeces.
 * It can be viewed as a compression model like Lempel-Ziv 77, but instead of
 * exact matches, allowing matches that contain inserts/deletes/mismatches.
 *
 */

package misc.fuzzyLZ;

import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import misc.common.MarkovN;
import misc.common.Misc;
import misc.common.Params;
import misc.common.Seq_Model;
import misc.dnaPlatform.OptionsHandle;

public class FuzzyDriver implements Serializable {
	private static final long serialVersionUID = 1L;

	static final String VERSION = "1.3";

	private static final double CONVERGE_CUTOFF = 0.0001;

	// If any object variables are added here, be use to think about
	// serialisation of them...
	// if the variable is needed after a 'resume', add the variable to the
	// writeObject() and readObject() functions.

	int maxIterations;

	char alphabet[];

	int DEBUG;

	String fname;

	int imageFreq;

	int checkpointFreq;

	int statsFreq;

	String msgFname;

	String outDir;

	String fprefix;

	// FileWriter msgFile;
	PrintStream msgOut;

	Params params;

	char[] str;

	char[] preStr;

	char[] joinedStr; // Simple concatenation of preStr + str

	String seqModelStr; // Sequence model for base state.

	boolean overwrite; // Overwrite files

	static CommandLine defaultCmdLine = new CommandLine();
	static {
		defaultCmdLine.addInt("maxIterations", 1,
				"Max number of iterations. (0 means until convergence)");
		defaultCmdLine.addString("preFile", "",
				"Prepend 'preFile' to sequence.");
		defaultCmdLine
				.addString(
						"seqModel",
						"markov(2)",
						"Base sequence model.  Use 'markov(n)' for a n-th order markov model, n=-1 for uniform");
		defaultCmdLine.addString("dna", "atgc",
				"Alphabet used by the sequence.");

		defaultCmdLine.addString("fwdMach", "3state",
				"Comma separated list of machines to use for forward matches.\n"
						+ "(Use an empty string '' for no forward machines.)\n"
						+ "Supported: 1state,3state");

		defaultCmdLine.addString("revMach", "3state",
				"Comma separated list of machines to use for reverse matches.\n"
						+ "(Use an empty string '' for no reverse machines.)\n"
						+ "Supported: 1state,3state");

		defaultCmdLine.addBoolean("overwrite", false, "Overwrite msglen file.");
		defaultCmdLine.addInt("debug", 2,
				"Debug level (higher gives more verbose  )");
		defaultCmdLine
				.addInt("imageSize", 1024, "Maximum Image size in pixels");
		defaultCmdLine.addInt("imageFreq", 0,
				"Save an image every <n> seconds.  (0 - to disable)");
		defaultCmdLine.addInt("checkFreq", 0,
				"Save a checkpoint every <n> seconds.  (0 - to disable)");
		defaultCmdLine.addInt("statsFreq", 300,
				"Display some stats every <n> seconds.  (0 - to disable)");

		defaultCmdLine.addString("msgFile", "",
				"Output file for encode length of each character.\n"
						+ "(The default is based on the input file name)");

		defaultCmdLine.addString("outDir", "." + File.separatorChar,
				"Directory to save output files in.");

		// These options are for Matches_Sparse
		defaultCmdLine
				.addInt("hashSize", 10,
						"Window size to use for constructing hashtable (0 - for full N^2 algorithm)");
		defaultCmdLine.addInt("computeWin", 10,
				"Number of cells to activate on a hashtable hit");
		defaultCmdLine
				.addInt("cutML", 4,
						"When (cell_value - base_cell > cutML) then cell is killed. (in bits)");
		// defaultCmdLine.addBoolean("plotActive", false, "true: plot only
		// active
		// cells, false: plot cell values");

		defaultCmdLine
				.addString("paramFile", "",
						"Parameter file to read for various model parameters (see docs)");

		defaultCmdLine.addBoolean("resume", false, "Resume from a checkpoint");
	}

	@SuppressWarnings("unchecked")
	public static void main(String args[]) throws Exception {
		CommandLine cmdLine = FuzzyDriver.defaultCmdLine;

		args = cmdLine.parseLine(args);

		System.out
				.println("Approximate  Repeat Model for DNA sequence compression");

		System.out.println("   L. Allison, T. Edgoose, T. I. Dix.");
		System.out
				.println("   Compression of Strings with Approximate Repeats");
		System.out.println("   Intell. Sys. in Molec. Biol., pp.8-16, 1998\n");

		if (args == null || args.length != 1) {
			System.err
					.println("Usage: java " + FuzzyDriver.class.getName()
							+ " [options] <seqFile|checkpointFile>\n"
							+ cmdLine.usageMessage());
			System.exit(1);
		}

		boolean resume = cmdLine.getBooleanVal("resume");

		FuzzyDriver me = null;

		// Start compression
		if (resume) {
			// Resuming from a checkpoint
			try {
				System.out.println("Reloading from checkpoint...");
				File f = new File(args[0]);
				ObjectInputStream is = new ObjectInputStream(
						new FileInputStream(f));
				me = (FuzzyDriver) is.readObject();
				System.out.println("Successfully loaded checkpoint...");
				is.close();
			} catch (Exception e) {
				System.err.println("Unable to resume from checkpoint: " + e);
				System.exit(1);
			}

			// use defaults.
			if (cmdLine.optionSet("maxIterations"))
				me.maxIterations = cmdLine.getIntVal("maxIterations");
			if (cmdLine.optionSet("imageFreq"))
				me.imageFreq = cmdLine.getIntVal("imageFreq");
			if (cmdLine.optionSet("checkFreq"))
				me.checkpointFreq = cmdLine.getIntVal("checkFreq");
			if (cmdLine.optionSet("statsFreq"))
				me.statsFreq = cmdLine.getIntVal("statsFreq");
			if (cmdLine.optionSet("imageSize"))
				FuzzyLZ.img_width = FuzzyLZ.img_height = cmdLine
						.getIntVal("imageSize");
			if (cmdLine.optionSet("debug")) {
				me.DEBUG = cmdLine.getIntVal("debug");
				FuzzyLZ.DEBUG = me.DEBUG;
			}

		}// end if resume
		else {// start from scratch
				// Starting anew (no checkpoint)
			me = new FuzzyDriver();
			me.params = new Params();

			me.maxIterations = cmdLine.getIntVal("maxIterations");
			me.alphabet = cmdLine.getStringVal("dna").toCharArray();
			me.DEBUG = cmdLine.getIntVal("debug");
			FuzzyLZ.DEBUG = me.DEBUG;
			FuzzyLZ.img_width = FuzzyLZ.img_height = cmdLine
					.getIntVal("imageSize");
			me.imageFreq = cmdLine.getIntVal("imageFreq");
			me.checkpointFreq = cmdLine.getIntVal("checkFreq");
			me.statsFreq = cmdLine.getIntVal("statsFreq");
			me.msgFname = cmdLine.getStringVal("msgFile");
			me.outDir = cmdLine.getStringVal("outDir");
			me.overwrite = cmdLine.getBooleanVal("overwrite");
			Matches_Sparse.def_winSize = cmdLine.getIntVal("hashSize");
			Matches_Sparse.def_computeWin = cmdLine.getIntVal("computeWin");
			Matches_Sparse.def_cutML = cmdLine.getIntVal("cutML");

			// Matches_Sparse.def_plotActive =
			// cmdLine.getBooleanVal("plotActive");

			String paramFile = cmdLine.getStringVal("paramFile");
			if (paramFile.length() > 0) {
				try {
					BufferedReader rdr = new BufferedReader(new FileReader(
							paramFile));
					String line;
					while ((line = rdr.readLine()) != null) {
						int i = line.indexOf('=');
						if (i >= 0) {
							String key = line.substring(0, i);
							String valStr = line.substring(i + 1);
							double val = Double.parseDouble(valStr);
							me.params.put(key, val);
						} else {
							if (line.length() > 0 && !line.startsWith("#"))
								System.err
										.println("WARNING: Ignoring param line '"
												+ line + "'");
						}
					}
					rdr.close();
				} catch (FileNotFoundException e) {
					System.err.println("ERROR: Unable to read file '"
							+ paramFile + "'.  Ignoring...");
				} catch (NumberFormatException e) {
					System.err
							.println("ERROR: Converting string from paramFile to number :"
									+ e);
					me.params = new Params();
				} catch (IOException e) {
					System.err.println("ERROR: Reading file '" + paramFile
							+ "'.  Ignoring...");
					me.params = new Params();
				}
			}

			me.seqModelStr = cmdLine.getStringVal("seqModel");

			@SuppressWarnings("rawtypes")
			Vector machs = new Vector();
			try {
				FuzzyLZ.def_numFwd = parseMachineNames(
						cmdLine.getStringVal("fwdMach"), machs);
				FuzzyLZ.def_numRev = parseMachineNames(
						cmdLine.getStringVal("revMach"), machs);
			} catch (Exception e) {
				System.err.println(e);
				System.exit(1);
			}

			FuzzyLZ.MutationModels = 
					(String[]) machs.toArray(new String[machs.size()]);

			String preFile = cmdLine.getStringVal("preFile");
			if (preFile == "") {
				me.preStr = new char[0];
			} else {
				Sequence seq = SequenceReader.getReader(preFile).nextSequence(
						null);// (filename)IOTools.read(args[0]);
				if (seq == null) {
					System.err.println("Unable to read prefix file: '"
							+ preFile + "'");
					System.exit(1);
				}
				me.preStr = seq.charSequence();

				for (int i = 0; i < me.preStr.length; i++)
					me.preStr[i] = Character.toLowerCase(me.preStr[i]);

				Misc.printf("Pre-sequence length = %d\n", me.preStr.length);
			}

			me.fname = args[0];
			Sequence seq = SequenceReader.getReader(me.fname)
					.nextSequence(null);// (filename)IOTools.read(args[0]);
			if (seq == null) {
				System.err.println("Unable to read sequence file");
				System.exit(1);
			}
			me.str = seq.charSequence();
			for (int i = 0; i < me.str.length; i++)
				me.str[i] = Character.toLowerCase(me.str[i]);

			Misc.printf("Sequence length = %d\n", me.str.length);

			// Display the machines we are going to use
			for (int i = 0; i < FuzzyLZ.def_numFwd; i++)
				System.err.println("FwdMachine[" + i + "]: "
						+ FuzzyLZ.MutationModels[i]);
			for (int i = 0; i < FuzzyLZ.def_numRev; i++)
				System.err.println("RevMachine[" + i + "]: "
						+ FuzzyLZ.MutationModels[i + FuzzyLZ.def_numFwd]);

			// Create joinedStr
			me.joinedStr = new char[me.preStr.length + me.str.length];
			System.arraycopy(me.preStr, 0, me.joinedStr, 0, me.preStr.length);
			System.arraycopy(me.str, 0, me.joinedStr, me.preStr.length,
					me.str.length);

			// Get the file prefix to use for output filenames
			me.fprefix = me.outDir;
			if (me.fprefix.length() > 0)
				me.fprefix += File.separatorChar;
			me.fprefix += (new File(me.fname)).getName();

			// Open the file for msglen output
			// try {
			// File f;
			if (me.msgFname.equals(""))
				me.msgFname = me.fprefix + "-msglen.txt";
			// f = new File(me.msgFname);
			// if (me.overwrite)
			// f.delete();

			// if (f.exists()) {
			// System.err.println("Output file '" + me.msgFname
			// + "' already exists.");
			// System.exit(1);
			// }
			// if (!f.createNewFile()) {
			// System.err.println("Unable to create Output file '"
			// + me.msgFname);
			// System.exit(1);
			// }
			//
			// me.msgFile = new FileWriter(f);
			// Create new file output stream
			// me.msgOut = new PrintStream(f);
			// }
			// catch (IOException e) {
			// System.err.println("Error creating msglen output file: " + e);
			// System.exit(1);
			// }
		}
		try {
			me.go(resume);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * This Function is called by FuzzyModel to run FuzzyLZ given a char[]
	 * instead of reading sequence from a file.
	 * 
	 * @param myOptions
	 *            OptionsHandle
	 * @param myFile
	 *            String
	 * @param sequence
	 *            char[]
	 * @throws RuntimeException
	 * @throws IOException
	 * @return filename of Information content sequence created
	 */
	@SuppressWarnings("unchecked")
	public String start(OptionsHandle myOptions, String myFile, char[] sequence)
			throws RuntimeException, IOException {

		String resume = myOptions.getStringValue("resume");

		FuzzyDriver me = null;

		if (resume.length() > 0) {
			// Resuming from a checkpoint
			try {
				System.out.println("Reloading from checkpoint...");
				File f = new File(myOptions.getStringValue("resume"));
				ObjectInputStream is = new ObjectInputStream(
						new FileInputStream(f));
				me = (FuzzyDriver) is.readObject();
				System.out.println("Successfully loaded checkpoint...");
				is.close();
			} catch (Exception e) {
				throw new IOException("Unable to resume from checkpoint: " + e);
			}
			// We can re-set some of the parameters here.
			// Only reset ones that are defined on this commandline, ie. don't
			// use defaults.
			if (myOptions.optionSet("maxIterations"))
				me.maxIterations = myOptions.getIntValue("maxIterations");
			if (myOptions.optionSet("imageFreq"))
				me.imageFreq = myOptions.getIntValue("imageFreq");
			if (myOptions.optionSet("checkFreq"))
				me.checkpointFreq = myOptions.getIntValue("checkFreq");
			if (myOptions.optionSet("statsFreq"))
				me.statsFreq = myOptions.getIntValue("statsFreq");
			if (myOptions.optionSet("imageSize"))
				FuzzyLZ.img_width = FuzzyLZ.img_height = myOptions
						.getIntValue("imageSize");
			if (myOptions.optionSet("debug")) {
				me.DEBUG = myOptions.getIntValue("debug");
				FuzzyLZ.DEBUG = me.DEBUG;
			}

		} else {
			// Starting anew (no checkpoint)
			me = new FuzzyDriver();
			me.params = new Params();

			me.maxIterations = myOptions.getIntValue("maxIterations");
			me.alphabet = myOptions.getStringValue("dna").toCharArray();
			me.DEBUG = myOptions.getIntValue("debug");
			FuzzyLZ.DEBUG = me.DEBUG;
			FuzzyLZ.img_width = FuzzyLZ.img_height = myOptions
					.getIntValue("imageSize");
			me.imageFreq = myOptions.getIntValue("imageFreq");
			me.checkpointFreq = myOptions.getIntValue("checkFreq");
			me.statsFreq = myOptions.getIntValue("statsFreq");
			me.msgFname = myOptions.getStringValue("msgFile");
			me.outDir = myOptions.getStringValue("outDir");
			me.overwrite = myOptions.getBooleanValue("overwrite");
			Matches_Sparse.def_winSize = myOptions.getIntValue("hashSize");
			Matches_Sparse.def_computeWin = myOptions.getIntValue("computeWin");
			Matches_Sparse.def_cutML = myOptions.getIntValue("cutML");
			// Matches_Sparse.def_plotActive =
			// myOptions.getBooleanValue("plotActive");

			String paramFile = myOptions.getStringValue("paramFile");
			if (paramFile.length() > 0) {
				try {
					BufferedReader rdr = new BufferedReader(new FileReader(
							paramFile));
					String line;
					while ((line = rdr.readLine()) != null) {
						int i = line.indexOf('=');
						if (i >= 0) {
							String key = line.substring(0, i);
							String valStr = line.substring(i + 1);
							double val = Double.parseDouble(valStr);
							me.params.put(key, val);
						} else {
							if (line.length() > 0 && !line.startsWith("#"))
								System.err
										.println("WARNING: Ignoring param line '"
												+ line + "'");
						}
					}
					rdr.close();
				} catch (FileNotFoundException e) {
					System.err.println("ERROR: Unable to read file '"
							+ paramFile + "'.  Ignoring...");
				} catch (NumberFormatException e) {
					System.err
							.println("ERROR: Converting string from paramFile to number :"
									+ e);
					me.params = new Params();
				} catch (IOException e) {
					System.err.println("ERROR: Reading file '" + paramFile
							+ "'.  Ignoring...");
					me.params = new Params();
				}
			}

			me.seqModelStr = myOptions.getStringValue("seqModel");

			@SuppressWarnings("rawtypes")
			Vector machs = new Vector();
			FuzzyLZ.def_numFwd = parseMachineNames(
					myOptions.getStringValue("fwdMach"), machs);
			FuzzyLZ.def_numRev = parseMachineNames(
					myOptions.getStringValue("revMach"), machs);
			FuzzyLZ.MutationModels = 
					(String[]) machs.toArray(new String[machs.size()]);

			String preFile = myOptions.getStringValue("preFile");
			if (preFile == "") {
				me.preStr = new char[0];
			} else {
				Sequence seq = SequenceReader.getReader(preFile).nextSequence(
						null);// (filename)IOTools.read(args[0]);
				if (seq == null) {
					throw new IOException("Unable to read prefix file: '"
							+ preFile + "' ");
				}
				me.preStr = seq.charSequence();

				for (int i = 0; i < me.preStr.length; i++)
					me.preStr[i] = Character.toLowerCase(me.preStr[i]);

				Misc.printf("Pre-sequence length = %d\n", me.preStr.length);
			}

			me.fname = myFile;
			/*
			 * DNA japsa.seq = DNA.guess_format(me.fname); if (japsa.seq ==
			 * null) { System.err.println("Unable to read sequence file");
			 * System.exit(1); }
			 */
			me.str = sequence;
			Misc.printf("Sequence length = %d\n", me.str.length);

			// Display the machines we are going to use
			for (int i = 0; i < FuzzyLZ.def_numFwd; i++)
				System.err.println("FwdMachine[" + i + "]: "
						+ FuzzyLZ.MutationModels[i]);
			for (int i = 0; i < FuzzyLZ.def_numRev; i++)
				System.err.println("RevMachine[" + i + "]: "
						+ FuzzyLZ.MutationModels[i + FuzzyLZ.def_numFwd]);

			// Create joinedStr
			me.joinedStr = new char[me.preStr.length + me.str.length];
			System.arraycopy(me.preStr, 0, me.joinedStr, 0, me.preStr.length);
			System.arraycopy(me.str, 0, me.joinedStr, me.preStr.length,
					me.str.length);

			// Get the file prefix to use for output filenames
			me.fprefix = me.outDir;
			if (me.fprefix.length() > 0)
				me.fprefix += File.separatorChar;
			me.fprefix += (new File(me.fname)).getName();

			// Open the file for msglen output
			try {
				File f;
				if (me.msgFname.equals(""))
					me.msgFname = me.fprefix + "-msglen.txt";
				f = new File(me.msgFname);
				if (me.overwrite)
					f.delete();
				if (f.exists()) {
					throw new IOException("Output file '" + me.msgFname
							+ "' already exists.");
				}
				if (!f.createNewFile()) {
					throw new IOException("Unable to create Output file '"
							+ me.msgFname);
				}
				// me.msgFile = new FileWriter(f);
				me.msgOut = new PrintStream(f);
			} catch (IOException e) {
				throw new IOException("Error creating msglen output file: " + e);
			}
		}
		me.go(resume.length() > 0);

		// return filename where information content was written
		return me.msgFname;

	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	static int parseMachineNames(String l, Vector res) throws IOException {
		String s[] = l.split(",");
		int num = 0;
		for (int i = 0; i < s.length; i++) {
			if (s[i].equals(""))
				continue;
			if (s[i].compareToIgnoreCase("1state") == 0) {
				res.add("misc.common.Mutation_1State");
			} else if (s[i].compareToIgnoreCase("3state") == 0) {
				res.add("misc.common.Mutation_3State");
			} else {
				throw new IOException("Unknown machine: '" + s[i] + "'");
			}
			num++;
		}
		return num;
	}

	// Note: these object variables must also be saved/loaded in the
	// writeObject/readObject functions.

	FuzzyLZ mdl;

	double tot_msglen;
	double last_msglen = -1;

	int iteration;

	int inner_i;

	long iterationStartTime; // Start time

	void init_iteration() throws IOException {
		if (DEBUG >= 3) {
			// TODO check this pls
			// S.printf("Starting (unnormalized) costs:\n%s\n", params);
		}

		// Ensure Params 'p' has no funky numbers
		int numParams = params.get_num();
		for (int i = 0; i < numParams; i++) {
			String name = params.get_name_by_id(i);
			double v = params.get(name);
			if (Double.isInfinite(v) || Double.isNaN(v)) {
				System.err.println("Potential parameter problem: '" + name
						+ "' has bad value=" + v);
			}
		}

		Seq_Model seqModel = parseSeqModel(seqModelStr);
		if (seqModel == null) {
			throw new IOException("ERROR: Unable to parse seqModel : '"
					+ seqModelStr + "'");
		}

		// I don't think the Buffered_Seq_Model actually
		// saves any computation!
		// seqModel = new Buffered_Seq_Model(seqModel); // reduces calls to
		// encodeLen

		// Train the seqModel up on preStr if it has been specified
		for (int i = 0; i < preStr.length; i++) {
			seqModel.update(preStr[i], i);
		}

		mdl = new FuzzyLZ(params, seqModel, joinedStr, alphabet.length,
				preStr.length);

		tot_msglen = 0;
		iterationStartTime = System.currentTimeMillis();
	}

	/**
	 * Parse a string representing a sequence model, and create the model.
	 * 
	 * @param tandemRepeat
	 * @return The sequence model. null on error
	 */
	private Seq_Model parseSeqModel(String str_) {
		Pattern p = Pattern.compile("markov\\((-?\\d+)\\)",
				Pattern.CASE_INSENSITIVE);
		Matcher m = p.matcher(str_);
		if (m.matches()) {
			int n = Integer.parseInt(m.group(1));
			return new MarkovN(n, alphabet);
		}

		return null;
	}

	void inner_loop() throws IOException {
		char c = str[inner_i];
		double d = mdl.update(c, preStr.length + inner_i);

		tot_msglen += d;
		Misc.my_assert(d > 0, "Bugger! -ve bits to encode char");

		if (DEBUG >= 3)
			System.out.printf("%c %03d: m=%.2f tot(m)=%.2f\n", c, inner_i, d,
					tot_msglen);

		// try {
		msgOut.printf("%c\t%7d\t%f\n", c, inner_i, d);

		// msgFile.write(VNTRReadDepth.sprintf("%s %03d %f\n", new
		// VNTRReadDepth.VarArgs(c).add(inner_i).add(d)));
		// msgFile.flush();
		// }
		// catch (IOException e) {
		// throw new IOException("Error writing to msglen output file: " + e);
		// }
	}

	/**
	 * The main function to compress
	 * 
	 * @param resume
	 */
	void go(boolean resume) throws IOException {
		if (resume) {
			System.out.println("Resuming from iteration=" + iteration
					+ " at character=" + inner_i + " of " + str.length);
		} else {
			iteration = 0;
		}

		long last_checkpoint = System.currentTimeMillis();
		long last_image = System.currentTimeMillis();
		long last_stats = System.currentTimeMillis();

		for (; maxIterations == 0 || iteration < maxIterations; iteration++) {

			// Re initilise the new file to write
			msgOut = new PrintStream(new FileOutputStream(this.msgFname
					+ iteration));
			msgOut.println("# Information content generated by Approximate Repeat Model");
			msgOut.println("#  Compression of Strings with Approximate Repeats "
					+ "\n# L. Allison, T. Edgoose, T. I. Dix., Intell. Sys. in Molec. Biol., pp.8-16, 1998.");

			if (DEBUG >= 1)
				Misc.printf("\n\nIteration %d\nParams:\n%s\n", new Object[] {
						new Integer(iteration), params });

			if (!resume) {
				try {
					init_iteration();
					inner_i = 0;
				} catch (Exception e) {
					System.err.println(e);
					return;
				}
			}

			resume = false;

			for (; inner_i < str.length; inner_i++) {
				// Display some stats?
				if (statsFreq > 0
						&& (System.currentTimeMillis() - last_stats) / 1000 >= statsFreq) {
					last_stats = System.currentTimeMillis();
					System.out.println("Stats as at " + new java.util.Date());
					System.out.println("Current compression: "
							+ (tot_msglen / inner_i) + " bits/char.");

					System.out.println("At " + inner_i + " of " + str.length
							+ " (" + (100.0 * inner_i / str.length) + "%)");
					double eta = 1.0
							* (System.currentTimeMillis() - iterationStartTime)
							* str.length * str.length
							/ (1.0 * inner_i * inner_i);
					int secs = (int) (eta / 1000);
					int mins = secs / 60;
					int hours = mins / 60;
					System.out.println("iteration ETA: " + hours + " hrs "
							+ (mins % 60) + " mins " + (secs % 60) + " secs.");
					mdl.display_stats();
				}

				// Save a checkpoint?
				if (checkpointFreq > 0
						&& (System.currentTimeMillis() - last_checkpoint) / 1000 >= checkpointFreq) {
					try {
						String f = fprefix + "-checkpoint-" + iteration + "-"
								+ inner_i + ".obj";
						if (DEBUG >= 1)
							System.out.println("Saving checkpoint : " + f);
						ObjectOutputStream o = new ObjectOutputStream(
								new FileOutputStream(new File(f)));
						o.writeObject(this);
						o.close();
					} catch (IOException e) {
						System.err.println("Failed to save checkpoint: " + e);
					}
					last_checkpoint = System.currentTimeMillis();
				}

				// Save an image?
				if (imageFreq > 0
						&& (System.currentTimeMillis() - last_image) / 1000 >= imageFreq) {
					mdl.plot.save(Misc.sprintf(fprefix + "-tmp-%02d-%07d.ppm",
							new Misc.VarArgs(iteration).add(inner_i)), "Seq: "
							+ fname);
					mdl.plotActive.save(Misc.sprintf(fprefix
							+ "-tmpActive-%02d-%07d.ppm", new Misc.VarArgs(
							iteration).add(inner_i)), "Seq: " + fname);
					mdl.plotHits.save(Misc.sprintf(fprefix
							+ "-tmpHits-%02d-%07d.ppm", new Misc.VarArgs(
							iteration).add(inner_i)), "Seq: " + fname);
					last_image = System.currentTimeMillis();
				}
				// Do the work here
				try {
					inner_loop();
				} catch (Exception e) {
					System.err.println(e);
					return;
				}
			}

			if (DEBUG >= 0)
				Misc.printf("Iteration " + iteration
						+ " : Total for mdl = %.4f  (%.4f bits/char)\n",
						tot_msglen, tot_msglen / str.length);

			if (DEBUG >= 1 || statsFreq > 0) {
				System.out.println("Stats as at " + new java.util.Date());
				mdl.display_stats();
			}

			if (DEBUG >= 2)
				mdl.display();

			params = mdl.counts_to_params();

			mdl.plot.save(Misc.sprintf(fprefix + "-final-iter%02d.ppm",
					new Misc.VarArgs(iteration)), "Seq: " + fname);
			mdl.plotActive.save(
					Misc.sprintf(fprefix + "-finalActive-iter%02d.ppm",
							new Misc.VarArgs(iteration)), "Seq: " + fname);
			mdl.plotHits.save(Misc.sprintf(fprefix + "-finalHits-iter%02d.ppm",
					new Misc.VarArgs(iteration)), "Seq: " + fname);

			if (last_msglen > 0) {
				if (last_msglen < tot_msglen)
					System.err.println("WARNING: Problem with convergence");
				else if (last_msglen - tot_msglen < CONVERGE_CUTOFF) {
					System.err
							.println("Converged to within " + CONVERGE_CUTOFF);
					break;
				}
			}
			last_msglen = tot_msglen;
			System.out.println("On average " + (tot_msglen * 1.0 / str.length));
			this.msgOut.close();
			File f = new File(this.msgFname + iteration), f2 = new File(
					this.msgFname);
			copy(f, f2);

			System.out.println("============="
					+ f.renameTo(new File(this.msgFname)));
			// Copy this into the file
		}
	}

	// Copies src file to dst file.
	// If the dst file does not exist, it is created
	static void copy(File src, File dst) throws IOException {
		InputStream in = new FileInputStream(src);
		OutputStream out = new FileOutputStream(dst);

		// Transfer bytes from in to out
		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0) {
			out.write(buf, 0, len);
		}
		in.close();
		out.close();
	}

	/**
	 * Used by FuzzyModel to get the image file names
	 * 
	 * @return String[]
	 */
	public String[] getImageFileNames() {
		String[] graphs = new String[3];
		graphs[0] = fprefix + "-final-iter%02d.ppm";
		graphs[1] = fprefix + "-finalActive-iter%02d.ppm";
		graphs[2] = fprefix + "-finalHits_iter%02d.ppm";

		return graphs;
	}

	// Write our own serization handler.
	// Just store everything. Must re-open the msglen output file
	private void writeObject(java.io.ObjectOutputStream out) throws IOException {
		out.writeInt(maxIterations);
		out.writeObject(alphabet);
		out.writeInt(DEBUG);
		out.writeObject(fname);
		out.writeInt(imageFreq);
		out.writeInt(checkpointFreq);
		out.writeInt(statsFreq);
		out.writeObject(msgFname);
		out.writeObject(outDir);

		out.writeObject(fprefix);
		out.writeObject(params);
		out.writeObject(str);
		out.writeObject(preStr);
		out.writeObject(seqModelStr);

		out.writeObject(mdl);
		out.writeDouble(tot_msglen);
		out.writeInt(iteration);
		out.writeInt(inner_i);
	}

	private void readObject(java.io.ObjectInputStream in) throws IOException,
			ClassNotFoundException {
		maxIterations = in.readInt();
		alphabet = (char[]) in.readObject();
		DEBUG = in.readInt();
		fname = (String) in.readObject();
		imageFreq = in.readInt();
		checkpointFreq = in.readInt();
		statsFreq = in.readInt();
		msgFname = (String) in.readObject();
		outDir = (String) in.readObject();

		fprefix = (String) in.readObject();
		params = (Params) in.readObject();
		str = (char[]) in.readObject();
		preStr = (char[]) in.readObject();
		seqModelStr = (String) in.readObject();

		mdl = (FuzzyLZ) in.readObject();
		tot_msglen = in.readDouble();
		iteration = in.readInt();
		inner_i = in.readInt();

		// re-Create joinedStr
		joinedStr = new char[preStr.length + str.length];
		System.arraycopy(preStr, 0, joinedStr, 0, preStr.length);
		System.arraycopy(str, 0, joinedStr, preStr.length, str.length);

		// Open the file for msglen output
		try {
			File f;
			f = new File(msgFname);
			if (!f.exists()) {
				throw new IOException("Output file '" + msgFname
						+ "' does not exist. Unable to resume");
			}
			// msgFile = new FileWriter(f, true);
			msgOut = new PrintStream(f);
			msgOut.println("# Resuming from checkpoint...\n");
			// msgFile.write("# Resuming from checkpoint...\n");
		} catch (IOException e) {
			throw new IOException("Error creating msglen output file: " + e);
		}
	}
}
