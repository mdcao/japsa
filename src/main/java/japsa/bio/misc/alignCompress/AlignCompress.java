/*  
 *  Copyright (c) David Powell <david@drp.id.au>
 *
 * This file is part of AlignCompress.
 *
 * AlignCompress aligns two sequences using a modified DPA
 * that allows the use of a _model_ for each sequence.
 *

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

 */

package japsa.bio.misc.alignCompress;

import japsa.bio.misc.common.*;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;

class BlendModel implements Seq_Model {
	private static final long serialVersionUID = 1L;

	Seq_Model models[];
	double m_lens[];

	public BlendModel(char[] alpha) {
		MarkovN_exact m0 = new MarkovN_exact(-1, alpha);
		MarkovN_exact m1 = new MarkovN_exact(1, alpha);

		m1.setProb("aa", 0.28);
		m1.setProb("at", 0.7);
		m1.setProb("ag", 0.01);
		m1.setProb("ac", 0.01);

		m1.setProb("ta", 0.7);
		m1.setProb("tt", 0.28);
		m1.setProb("tg", 0.01);
		m1.setProb("tc", 0.01);

		m1.setProb("ga", 0.49);
		m1.setProb("gt", 0.49);
		m1.setProb("gg", 0.01);
		m1.setProb("gc", 0.01);

		m1.setProb("ca", 0.49);
		m1.setProb("ct", 0.49);
		m1.setProb("cg", 0.01);
		m1.setProb("cc", 0.01);

		models = new Seq_Model[2];
		m_lens = new double[2];
		models[0] = m0;
		models[1] = m1;
		m_lens[0] = 1;
		m_lens[1] = 1;
	}

	public double encodeLen(char a, int i) {
		double sum = Double.POSITIVE_INFINITY;
		double weights[] = new double[models.length];

		for (int j = 0; j < models.length; j++)
			sum = MyMath.logplus(sum, m_lens[j]);

		for (int j = 0; j < models.length; j++)
			weights[j] = MyMath.exp2(sum - m_lens[j]);

		double p = 0;
		for (int j = 0; j < models.length; j++)
			p += weights[j] * MyMath.exp2(-models[j].encodeLen(a, i));

		return -MyMath.log2(p);
	}

	public double update(char a, int i) {
		double res = encodeLen(a, i);

		for (int j = 0; j < models.length; j++) {
			m_lens[j] += models[j].update(a, i);
			m_lens[j] /= 1.5;
		}

		return res;
	}
}

/**
 * BufferModel takes a model and a sequence and the sequence dna. It
 * precomputes the encoding length of every dna character for every
 * position in the sequence. The cumulative encoding length is also computed
 * along the sequence.
 **/
class BufferModel implements Seq_Model {
	private static final long serialVersionUID = 1L;

	String str;
	char[] alphabet;
	double[][] enc;
	double[] encCumulative;

	BufferModel(Seq_Model model, String str, char[] alpha) {
		this.str = str;
		this.alphabet = alpha;

		enc = new double[str.length()][alphabet.length];
		encCumulative = new double[str.length() + 1];
		encCumulative[0] = 0;

		for (int i = 0; i < str.length(); i++) {
			for (int j = 0; j < alphabet.length; j++) {
				enc[i][j] = model.encodeLen(alphabet[j], i);
			}
			model.update(str.charAt(i), i);
			encCumulative[i + 1] = encCumulative[i]
					+ enc[i][char2int(str.charAt(i))];
		}

		// Test 'model' is efficient and doesn't cheat.
		for (int i = 0; i < str.length(); i++) {
			double s = Double.POSITIVE_INFINITY;
			for (int j = 0; j < alphabet.length; j++) {
				s = MyMath.logplus(s, enc[i][j]);
			}
			Misc.my_assert(Math.abs(s) < 1E-10, "logplus() != 0.  s=" + s);
		}
	}

	public String toString() {
		StringBuffer s = new StringBuffer();
		for (int i = 0; i < str.length(); i++) {
			s.append(str.charAt(i) + "    ");
			for (int j = 0; j < alphabet.length; j++) {
				double n = Math.rint(enc[i][j] * 100) / 100;
				s.append(alphabet[j] + ": " + n + " ");
			}
			s.append("\n");
		}
		return s.toString();
	}

	int char2int(char a) {
		for (int i = 0; i < alphabet.length; i++)
			if (a == alphabet[i])
				return i;
		Misc.my_assert(false, "Charater '" + a + "' not in defined dna");
		return -1;
	}

	public double encodeLen(char a, int i) {
		return enc[i][char2int(a)];
	}

	public double update(char a, int i) {
		return enc[i][char2int(a)];
	}

	public double encodeCumulative(int i) {
		return encCumulative[i];
	}
}

class AlignCompress {
	char[] alphabet;
	String seqA;
	String seqB;

	String paramString;
	int markovOrder;
	int maxIterations;
	int verbose;
	boolean useBlendModel;
	boolean linearCosts;
	boolean localAlign;
	boolean sumAlignments;

	boolean doTraceBack;

	// Encode length l: 0..infinity
	static private double encode_length(double l) {
		Misc.my_assert(l >= 0, "Bad length to encode:" + l);
		// return MyMath.logstar_continuous(l+1);
		// return -(l*MyMath.log2(PcharEncode) + MyMath.log2(1-PcharEncode)); //
		// Geometric distribution
		// return 0.104808 * l; // This is to match SW costs. (TESTING!)
		return 0;
	}

	public AlignCompress() {
	}

	// Return the appropriate cell of D[][]. Use i%2 if only keeping 2 rows
	private Mutation_FSM cell(Mutation_FSM D[][], int i, int j) {
		if (doTraceBack) {
			return D[i][j];
		} else {
			return D[i % 2][j];
		}
	}

	private String getAlignment(Mutation_FSM D[][],
			Mutation_FSM.TraceBack_Info fcell) {
		Mutation_FSM.TraceBack_Data d = fcell.get_tbdata();
		int i = d.i;
		int j = d.j;

		StringBuffer resA = new StringBuffer();
		StringBuffer resB = new StringBuffer();

		if (i < 0 && j < 0) {
			d = fcell.get_from(d);
			i = d.i;
			j = d.j;
		}

		int end_i = i, end_j = j;

		String endA = seqA.substring(i);
		String endB = seqB.substring(j);

		while (true) {
			// System.err.println("\n\n(i,j)="+i+","+j+"\nA="+resA+"\nB="+resB+"\n");

			d = ((Mutation_FSM.TraceBack_Info) D[i][j]).get_from(d);

			if (d == null)
				break;

			resA.append((d.i == i) ? '-' : seqA.charAt(d.i));
			resB.append((d.j == j) ? '-' : seqB.charAt(d.j));

			i = d.i;
			j = d.j;
		}

		System.out.println("ALIGNMENT CUTS: A[" + i + ".." + end_i + "] B[" + j
				+ ".." + end_j + "]");

		resA.reverse();
		resB.reverse();

		StringBuffer res = new StringBuffer();

		// Pretty up the alignment. For local alignments, put the front and end
		// bits on.
		for (int x = 0; x < j - (i < j ? i : j); x++)
			res.append(" ");
		res.append(seqA.substring(0, i));
		if (i > 0 || j > 0)
			res.append("   ");
		res.append(resA.toString().toUpperCase());
		res.append("   ");
		res.append(endA);

		res.append("\n");

		for (int x = 0; x < i - (i < j ? i : j); x++)
			res.append(" ");
		res.append(seqB.substring(0, j));
		if (i > 0 || j > 0)
			res.append("   ");
		res.append(resB.toString().toUpperCase());
		res.append("   ");
		res.append(endB);

		return res.toString();
	}

	public double doAlign() {
		System.out.println("# SeqA = "
				+ seqA
				+ "\n# SeqB = "
				+ seqB
				+ "\n# japsa.seq model = "
				+ (useBlendModel ? "blend model" : "markov order "
						+ markovOrder) + "\n# max iterations=" + maxIterations
				+ "\n# linearCosts=" + linearCosts + "\n# sumAlignments="
				+ sumAlignments + "\n# local Alignment=" + localAlign
				+ "\n# verbosity=" + verbose + "\n# params=" + paramString
				+ "\n");

		Params p = new Params();
		p.fromString(paramString);

		BufferModel modelA, modelB;
		{
			Seq_Model mdlA, mdlB;
			if (useBlendModel) {
				mdlA = new BlendModel(alphabet);
				mdlB = new BlendModel(alphabet);
			} else {
				mdlA = new MarkovN_fitted(markovOrder, alphabet, seqA);
				mdlB = new MarkovN_fitted(markovOrder, alphabet, seqB);
			}

			modelA = new BufferModel(mdlA, seqA, alphabet);
			modelB = new BufferModel(mdlB, seqB, alphabet);
		}

		// System.err.println("modelA\n"+modelA+"\n");
		// System.err.println("modelB\n"+modelB+"\n");

		int mdlCounts = Model_SeqAB.required_counts();

		int fsmCounts = -1;
		if (!linearCosts && !sumAlignments)
			fsmCounts = Mutation_1State.One.required_counts();
		if (!linearCosts && sumAlignments)
			fsmCounts = Mutation_1State.All.required_counts();
		if (linearCosts && !sumAlignments)
			fsmCounts = Mutation_3State.One.required_counts();
		if (linearCosts && sumAlignments)
			fsmCounts = Mutation_3State.All.required_counts();

		Misc.my_assert(fsmCounts >= 0, "Bad number of fsmCounts");

		int totCounts = mdlCounts + fsmCounts;

		double bestDiff = Double.NEGATIVE_INFINITY;
		double lastAlignment = 0;

		int iter = 0;
		while (maxIterations < 0 || iter < maxIterations) {
			int countPos = 0;
			Two_Seq_Model_Counts model = new Model_SeqAB(p, modelA, modelB,
					countPos);
			countPos += mdlCounts;
			Mutation_FSM fsmType = null;

			if (!linearCosts && !sumAlignments)
				fsmType = new Mutation_1State.One(model, p, totCounts, countPos);
			if (!linearCosts && sumAlignments)
				fsmType = new Mutation_1State.All(model, p, totCounts, countPos);
			if (linearCosts && !sumAlignments)
				fsmType = new Mutation_3State.One(model, p, totCounts, countPos);
			if (linearCosts && sumAlignments)
				fsmType = new Mutation_3State.All(model, p, totCounts, countPos);

			Misc.my_assert(fsmType != null, "Unable to construct fsmType");

			countPos += fsmCounts;
			Misc.my_assert(countPos == totCounts,
					"Internal error: countPos!=totCounts");

			Counts initialCounts = new Counts(totCounts);
			for (int i = 0; i < totCounts; i++)
				initialCounts.inc(i, 0.5); // Initialise all counts.

			// Determine if fsmType supports traceBack. ie. Does it implement
			// TraceBack_Info
			doTraceBack = false;
			if (fsmType instanceof Mutation_FSM.TraceBack_Info)
				doTraceBack = true;

			// Setup the DPA matrix
			Mutation_FSM D[][];
			Mutation_FSM final_cell;
			if (!doTraceBack) {
				// No traceback info, only keep 2 rows in D[][]
				D = new Mutation_FSM[2][seqB.length() + 1];
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < seqB.length() + 1; j++) {
						D[i][j] = (Mutation_FSM) fsmType.clone();
					}
				}
				final_cell = (Mutation_FSM) fsmType.clone();
			} else {
				// Keep all traceback info, so keep all rows in D[][]
				D = new Mutation_FSM[seqA.length() + 1][seqB.length() + 1];
				for (int i = 0; i < seqA.length() + 1; i++) {
					for (int j = 0; j < seqB.length() + 1; j++) {
						D[i][j] = (Mutation_FSM) fsmType.clone();
						Mutation_FSM.TraceBack_Info c = (Mutation_FSM.TraceBack_Info) D[i][j];
						c.set_tbdata(new Mutation_FSM.TraceBack_Data(i, j));
					}
				}
				final_cell = (Mutation_FSM) fsmType.clone();
				((Mutation_FSM.TraceBack_Info) final_cell)
						.set_tbdata(new Mutation_FSM.TraceBack_Data(-1, -1));
			}

			// Initialise the first cell of the DPA matrix
			cell(D, 0, 0).init_counts(initialCounts); // Initialise counts
			if (localAlign)
				cell(D, 0, 0).init_val(Double.POSITIVE_INFINITY);
			else
				cell(D, 0, 0).init_val(0);

			final_cell.init_val(Double.POSITIVE_INFINITY);

			if (verbose >= 1) {
				System.out.println("\n\nIteration: " + iter);
				System.out.println(model);
				System.out.println(cell(D, 0, 0).paramsToString());
			}

			// Do the DPA!
			// NB. The Mutation_FSM calc function calculates the value of this
			// cell
			// on the three neighbouring cell. This is the reverse of the way
			// the
			// DPA is usually expressed.
			for (int i = 0; i < seqA.length() + 1; i++) {

				// Reset the next row of the DPA matrix
				if (!doTraceBack)
					for (int j = 0; j < seqB.length() + 1; j++) {
						cell(D, i + 1, j).reset();
					}

				for (int j = 0; j < seqB.length() + 1; j++) {
					Mutation_FSM v = (i == seqA.length() ? null : cell(D,
							i + 1, j));
					Mutation_FSM h = (j == seqB.length() ? null : cell(D, i,
							j + 1));
					Mutation_FSM d = (i == seqA.length() || j == seqB.length() ? null
							: cell(D, i + 1, j + 1));
					char aChar = (i == seqA.length() ? '-' : seqA.charAt(i));
					char bChar = (j == seqB.length() ? '-' : seqB.charAt(j));

					if (localAlign) {
						double val;

						// Compute contribution of a local alignment that starts
						// at (i,j)
						val = modelA.encodeCumulative(i)
								+ modelB.encodeCumulative(j);
						val += encode_length(i);
						val += encode_length(j);
						cell(D, i, j).or(val, initialCounts);

						// Compute contribution of a local alignment that ends
						// at (i,j)
						val = cell(D, i, j).get_val()
								+ (modelA.encodeCumulative(seqA.length()) - modelA
										.encodeCumulative(i))
								+ (modelB.encodeCumulative(seqB.length()) - modelB
										.encodeCumulative(j));

						val += encode_length(seqA.length() - i);
						val += encode_length(seqB.length() - j);

						if (doTraceBack)
							((Mutation_FSM.TraceBack_Info) final_cell).or(val,
									cell(D, i, j).get_counts(), cell(D, i, j));
						else
							final_cell.or(val, cell(D, i, j).get_counts());
					}

					cell(D, i, j).calc(h, v, d, aChar, bChar, i, j);

					if (verbose >= 3) {
						System.err.println("Calc outputs from D[" + i + "]["
								+ j + "]");
						System.err.println(cell(D, i, j));
					}

				}
			}

			// If global alignment, then the final cell is really the bottom
			// right of D[][]
			if (!localAlign)
				final_cell = cell(D, seqA.length(), seqB.length());

			double encAlignModel = final_cell.encode_params();
			double encAlignment = encAlignModel + final_cell.get_val();
			if (localAlign) {
				// Add a cost for the length of the alignment.
				// Assume we know the length of the sequences. Need to encode
				// the start and end
				// of the alignment. Assume uniform over all positions. Have 4
				// cut-points to
				// encode, but only need to encode 3 (kind of).
				// So encode 2 from the shorter sequence, and 1 from the longer.
				double l1 = (seqA.length() < seqB.length() ? seqA.length()
						: seqB.length());
				double l2 = (seqA.length() > seqB.length() ? seqA.length()
						: seqB.length());
				encAlignment += MyMath.log2(l1) * 2 - 1;
				encAlignment += MyMath.log2(l2);
			}

			double encA = modelA.encodeCumulative(seqA.length());
			double encB = modelB.encodeCumulative(seqB.length());
			double encNull = encA + encB;

			if (localAlign) {
				// Encode lengths for the null theory
				// Note that global alignments ignore lengths, ignore for null
				// when doing global.
				encNull += encode_length(seqA.length());
				encNull += encode_length(seqB.length());
			}

			if (doTraceBack) {
				String s = getAlignment(D,
						(Mutation_FSM.TraceBack_Info) final_cell);
				System.out.println("ALIGNMENT:\n" + s);
			}

			// Get new parameters for next iteration
			p = fsmType.counts_to_params(final_cell.get_counts());
			/*
			 * if (localAlign) { double n = final_cell.get_counts().get(0);
			 * p.put("PalignChar", (n-1)/n); }
			 */

			if (verbose >= 2) {
				System.out.println("encA=" + encA + " encB=" + encB
						+ " encNull=" + (encNull));

				System.out.print("Mutual Encoding = " + encAlignment + " bits");
				System.out.println(" (model=" + encAlignModel + " data="
						+ (encAlignment - encAlignModel) + ")");
				System.out.println("\nCounts:\n" + final_cell.get_counts());
				System.out.println("\nEstimated Parameters:\n" + p);
			}
			System.out.println((encAlignment < encNull ? "related"
					: "unrelated")
					+ " ("
					+ (encNull - encAlignment)
					+ ")"
					+ "  log odds ratio = "
					+ (encNull - encAlignment)
					+ " bits");

			if (bestDiff < encNull - encAlignment)
				bestDiff = encNull - encAlignment;

			if (iter > 0 && verbose >= 1 && encAlignment > lastAlignment) {
				System.err.println("NON-CONVERGENCE: this=" + encAlignment
						+ " last=" + lastAlignment);
			}

			// Done we change little in alignment length?
			if (iter > 0 && encAlignment <= lastAlignment
					&& lastAlignment - encAlignment < 0.1)
				break;

			lastAlignment = encAlignment;
			iter++;
		} // End iterations

		return bestDiff;
	}

	public static void printLicense() {
		System.err.println("");
		System.err
				.println("  AlignCompress, Copyright (C) 2004 David Powell <david@drp.id.au>");
		System.err
				.println("  AlignCompress comes with ABSOLUTELY NO WARRANTY; and is provided");
		System.err
				.println("  under the GNU Public License v2, for details see file COPYRIGHT");
		System.err.println("");
		System.err.println("Please cite:");
		System.err.println("  L. Allison, D. R. Powell and T. I. Dix");
		System.err.println("  \"Compression and Approximate Matching\"");
		System.err.println("  The Computer Journal, 1999 42:1 pp1-10");
		System.err.println("");
	}

	public static void main(String args[]) throws Exception {
		printLicense();

		CommandLine cmdLine = new CommandLine();
		cmdLine.addInt("markov", 0,
				"Order of Markov Model to use for sequence models.");
		cmdLine.addInt("maxIterations", -1, "Maximum number of iterations.");
		cmdLine.addBoolean("blendModel", false,
				"Use blend model (a fixed model).");
		cmdLine.addBoolean("linearCosts", true, "Use linear gap costs.");
		cmdLine.addBoolean("sumAlignments", false, "Sum over all alignments.");
		cmdLine.addBoolean("local", true, "Compute using local alignments.");
		cmdLine.addInt("verbose", 0,
				"Display verbose output (larger num means more verbosity).");
		cmdLine.addBoolean("exSeq", false,
				"Command line options are explicit sequence, not filenames");
		cmdLine.addBoolean("protein", false, "Sequences are protein data");
		cmdLine.addString("params", "",
				"Params to pass to all classes (comma separated)");

		args = cmdLine.parseLine(args);

		AlignCompress a = new AlignCompress();
		boolean readFile = !cmdLine.getBooleanVal("exSeq");

		if (args == null || args.length != 2) {
			System.err
					.println("Usage: java AlignCompress [options] <seqA> <seqB>\n"
							+ cmdLine.usageMessage());
			System.exit(1);
		} else if (readFile && args.length == 2) {
			// Two filenames on commandline
			Sequence d = SequenceReader.getReader(args[0]).nextSequence(null);// (filename)IOTools.read(args[0]);
			a.seqA = d.toString();
			d = SequenceReader.getReader(args[1]).nextSequence(null);// (filename)IOTools.read(args[0]);
			a.seqB = d.toString();
		} else if (args.length == 2) {
			// Two strings on commandline, assume they are sequences
			a.seqA = args[0];
			a.seqB = args[1];
		}

		if (a.seqA == null || a.seqB == null) {
			System.err.println("Unable to read both sequences");
			System.exit(1);
		}

		a.seqA = a.seqA.toLowerCase();
		a.seqB = a.seqB.toLowerCase();

		a.markovOrder = cmdLine.getIntVal("markov");
		a.maxIterations = cmdLine.getIntVal("maxIterations");
		a.useBlendModel = cmdLine.getBooleanVal("blendModel");
		a.linearCosts = cmdLine.getBooleanVal("linearCosts");
		a.sumAlignments = cmdLine.getBooleanVal("sumAlignments");
		a.localAlign = cmdLine.getBooleanVal("local");
		a.verbose = cmdLine.getIntVal("verbose");
		a.paramString = cmdLine.getStringVal("params");

		if (cmdLine.getBooleanVal("protein")) {
			a.alphabet = new char[] { 'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g',
					'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v',
					'b', 'z', 'x', '*' };
		} else {
			a.alphabet = new char[] { 'a', 't', 'g', 'c' };
			// a.alphabet = new char[] {'a', 't', 'g', 'c',
			// 'y', 'r', 'w', 'n'};
		}

		System.out.println("-log odds ratio = " + a.doAlign() + " bits");

	} // End main()
}
