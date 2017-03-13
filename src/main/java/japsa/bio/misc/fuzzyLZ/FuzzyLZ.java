/*
 *  Copyright (c) David Powell <david@drp.id.au>
 *
 *
 * This file is part of FuzzyLZ
 *
 * FuzzyLZ is a program orginally intended for the
 * compression of DNA sequeces.  It can be viewed as a
 * compression model like Lempel-Ziv 77, but instead of
 * exact matches, allowing matches that contain
 * inserts/deletes/mismatches.
 *
 */

package japsa.bio.misc.fuzzyLZ;


import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;

import japsa.bio.misc.common.Counts;
import japsa.bio.misc.common.Misc;
import japsa.bio.misc.common.Mutation_FSM;
import japsa.bio.misc.common.MyMath;
import japsa.bio.misc.common.Params;
import japsa.bio.misc.common.Seq_Model;
import japsa.bio.misc.common.Two_Seq_Model_Counts;



/**
 * A model of 2 DNA sequences (for alignments) with counts. Has a seperate cost
 * for every pair of DNA characters. Thus 16 counts and parameters.
 */
class Model_SeqA_DNA implements Two_Seq_Model_Counts {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	static final int alphaSize = 4;

	char chars[] = { 'a', 't', 'g', 'c' };

	double costs[];

	int countIndex;

	public Model_SeqA_DNA(Params p, int alphaSize, int countIndex) {
		this.countIndex = countIndex;

		Misc.my_assert(chars.length == Model_SeqA_DNA.alphaSize,
				"Internal class error: char array and alphaSize don't match");
		Misc.my_assert(Model_SeqA_DNA.alphaSize == alphaSize,
				"Model_SeqA_DNA only works for an dna of "
						+ Model_SeqA_DNA.alphaSize);

		costs = new double[alphaSize * alphaSize];

		if (!p.exists("match_" + chars[0] + chars[0])) {
			set_default_costs();
		} else {
			for (int i = 0; i < alphaSize; i++)
				for (int j = 0; j < alphaSize; j++) {
					String s = "match_" + chars[i] + "|" + chars[j];
					costs[pair_to_int(chars[i], chars[j])] = p.get(s);
				}
		}

		normalize_costs();
	}

	private int pair_to_int(char a, char b) {
		int i, r = 0;
		for (i = 0; i < alphaSize; i++)
			if (a == chars[i]) {
				r = i;
				break;
			}
		Misc.my_assert(i < alphaSize, "Unknown character in DNA sequence");
		r *= alphaSize;

		for (i = 0; i < alphaSize; i++)
			if (b == chars[i]) {
				r += i;
				break;
			}
		Misc.my_assert(i < alphaSize, "Unknown character in DNA sequence");

		return r;
	}

	void set_default_costs() {
		for (int i = 0; i < alphaSize; i++) {
			for (int j = 0; j < alphaSize; j++) {
				costs[pair_to_int(chars[i], chars[j])] = (i == j ? 1 : 4);
			}
		}
	}

	void normalize_costs() {
		for (int j = 0; j < alphaSize; j++) {
			double sum = 0;
			for (int i = 0; i < alphaSize; i++)
				sum += MyMath.exp2(-costs[pair_to_int(chars[i], chars[j])]);
			for (int i = 0; i < alphaSize; i++)
				costs[pair_to_int(chars[i], chars[j])] += MyMath.log2(sum);
		}
	}

	public double encA(char a, int i) {
		return MyMath.log2(alphaSize);
	}

	public double encB(char a, int i) {
		return 0;
	}

	public double encBoth(char a, char b, int i, int j) {
		// if (a==b) return -MyMath.log2(0.8);
		// else return -MyMath.log2(0.2) + MyMath.log2(alphaSize-1);
		return costs[pair_to_int(a, b)];
	}

	public static int required_counts() {
		return alphaSize * alphaSize;
	}

	public Params counts_to_params(Counts counts) {
		Params par = new Params();

		for (int j = 0; j < alphaSize; j++) {
			double sum = 0;
			for (int i = 0; i < alphaSize; i++)
				sum += counts.counts[countIndex
						+ pair_to_int(chars[i], chars[j])];
			for (int i = 0; i < alphaSize; i++)
				par.put("match_" + chars[i] + "|" + chars[j],
						-MyMath.log2(counts.counts[countIndex
								+ pair_to_int(chars[i], chars[j])]
								/ sum));
		}

		return par;
	}

	public void update_count_encA(Counts c, double w, char a, int i) {
	};

	public void update_count_encB(Counts c, double w, char a, int i) {
	};

	public void update_count_encBoth(Counts c, double w, char a, char b, int i,
			int j) {
		c.inc(countIndex + pair_to_int(a, b), w);
	}

	public double encode_params(double N) {
		return 0;
	}
}

/**
 * Implements the FuzzyLZ algorithm of "Compression of Strings with Approximate
 * Repeats" by L. Allison, T. Edgoose, T.I. Dix in Intell. Sys. in Mol. Biol.
 * '98
 */
@SuppressWarnings("rawtypes")
public class FuzzyLZ implements Seq_Model {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	static int DEBUG = 5;

	/**
	 * Parent class for forward/reverse matches. Has parameters to continue a
	 * match (cont_cost) and to end a match (end_cost)
	 */
	static abstract class Matches implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		Mutation_FSM fsmType;

		Plot plot;

		Plot plotActive;

		Plot plotHits;

		int seqLen;

		char[] sequence;

		double encEnd, encContinue;

		int countIndex;

		protected final int contIndex = 0, endIndex = 1;

		Matches(Mutation_FSM fsm, Params p, int countIndex, char[] seq) {
			fsmType = fsm;
			this.countIndex = countIndex;
			sequence = seq;
			seqLen = sequence.length;

			if (!p.exists("end_cost")) {
				set_default_costs();
			} else {
				encEnd = p.get("end_cost");
				encContinue = p.get("cont_cost");
			}
			normalize_costs();

			plot = null;
		}

		// These three really oughtn't be here.
		void constructHash() {
		};

		void setHash(ExactMatches h) {
		};

		ExactMatches getHash() {
			return null;
		}

		void setPlot(Plot p) {
			plot = p;
		};

		void setPlotActive(Plot p) {
			plotActive = p;
		};

		void setPlotHits(Plot p) {
			plotHits = p;
		};

		public void plotVals(int i, double base) {
			return;
		};

		void set_default_costs() {

			// encEnd = -MyMath.log2(0.1);
			// encContinue = -MyMath.log2(0.9);
			// mdc changed to the values after a few iteration
			encEnd = 8.052274663541743;
			encContinue = 0.0054452480179498315;
		}

		void normalize_costs() {
			double sum = MyMath.exp2(-encEnd) + MyMath.exp2(-encContinue);
			encEnd = encEnd + MyMath.log2(sum);
			encContinue = encContinue + MyMath.log2(sum);
			if (DEBUG >= 2)
				System.out.printf("end_cost=%f cont_cost=%f\n", encEnd,
						encContinue);
			System.out.println("==end_cost=" + encEnd + " cont_cost="
					+ encContinue);
		}

		public abstract void beginLinks(int i, double startLen,
				Counts startCounts);

		public abstract double update(char a, int i, Counts retCounts);

		public abstract double msgLen();

		public abstract double normalise(double base);

		public static int required_counts() {
			return 2;
		};

		public Params counts_to_params(Counts c) {
			Params p = fsmType.counts_to_params(c);
			double sum = c.get(countIndex + endIndex)
					+ c.get(countIndex + contIndex);
			p.put("end_cost", -MyMath.log2(c.get(countIndex + endIndex) / sum));
			p.put("cont_cost",
					-MyMath.log2(c.get(countIndex + contIndex) / sum));
			return p;
		}
	}

	static int img_width = 800, img_height = 800;

	static String[] MutationModels = { "common.Mutation_3State$All" };

	static int def_numFwd = 1, def_numRev = 1;

	Seq_Model seqModel;

	int seqLen, alphaSize, totCounts;

	char[] sequence;

	int numFwd, numRev;

	double encNoStart;

	double[] encStartMachines;

	Matches[] machines;

	double base[];

	Counts baseCounts[];

	double last_msgLen;

	private final int noStartIndex = 0;

	private int myCounts = 1;

	public Plot plot;

	public Plot plotActive;

	public Plot plotHits;

	@SuppressWarnings("unchecked")
	FuzzyLZ(Params params, Seq_Model seqModel, char[] sequence, int alphaSize,
			int plotStart) throws RuntimeException {
		this.seqModel = seqModel;
		this.sequence = sequence;
		this.alphaSize = alphaSize;
		seqLen = sequence.length;
		numFwd = def_numFwd;
		numRev = def_numRev;

		Class[] M_class = new Class[numFwd + numRev];
		Method[] M_counts = new Method[numFwd + numRev];
		Constructor[] M_new = new Constructor[numFwd + numRev];
		for (int i = 0; i < M_class.length; i++) {
			try {
				M_class[i] = Class.forName(MutationModels[i
						% MutationModels.length]);
				M_counts[i] = M_class[i].getMethod("required_counts", null);
				M_new[i] = M_class[i].getConstructor(new Class[] {
						Class.forName("japsa.bio.misc.common.Two_Seq_Model_Counts"),
						Class.forName("japsa.bio.misc.common.Params"), Integer.TYPE,
						Integer.TYPE });
			} catch (Exception e) {
				throw new RuntimeException(
						"Error finding class/method for class '"
								+ MutationModels[i % MutationModels.length]
								+ "' : " + e);
			}
		}

		encStartMachines = new double[numFwd + numRev];
		machines = new Matches[numFwd + numRev];

		myCounts += numFwd + numRev;
		int countPos = myCounts;

		int mdlCounts = Model_SeqA.required_counts();
		int matCounts = Matches.required_counts();
		totCounts = myCounts;

		int fsmCounts[] = new int[numFwd + numRev];

		for (int i = 0; i < numFwd + numRev; i++) {
			try {
				fsmCounts[i] = ((Integer) M_counts[i].invoke(null, null))
						.intValue();
			} catch (Exception e) {
				System.err.println("required_counts Exception : " + e);
			}
			totCounts += (mdlCounts + fsmCounts[i] + matCounts);
		}

		if (!params.exists("nostart_cost")) {
			set_default_costs();
		} else {
			encNoStart = params.get("nostart_cost");
			for (int i = 0; i < numFwd + numRev; i++)
				encStartMachines[i] = params.get("start" + i + "_cost");
		}
		normalize_costs();

		plot = new Plot(seqLen + 1, seqLen + 1, img_width, img_height,
				plotStart);
		plotActive = new Plot(seqLen + 1, seqLen + 1, img_width, img_height,
				plotStart);
		plotHits = new Plot(seqLen + 1, seqLen + 1, img_width, img_height,
				plotStart);

		// Now initialise the various mutation machines
		for (int m = 0; m < numFwd + numRev; m++) {
			Params p = new Params();
			for (int i = 0; i < params.get_num(); i++) {
				String n = params.get_name_by_id(i);
				if (n.startsWith("MACHINE" + m + "_"))
					p.put(n.substring(9), params.get(n));
			}
			Two_Seq_Model_Counts model = new Model_SeqA(p, alphaSize, countPos);
			// Two_Seq_Model_Counts model = new Model_SeqA_SubModel(seqModel, p,
			// alphaSize, countPos);

			countPos += mdlCounts;
			Mutation_FSM fsmType = null;
			try {
				fsmType = (Mutation_FSM) M_new[m]
						.newInstance(new Object[] { model, p,
								new Integer(totCounts), new Integer(countPos) });
			} catch (Exception e) {
				System.err.println("newInstance Exception : " + e);
			}
			countPos += fsmCounts[m];
			machines[m] = new Matches_Sparse((m < numFwd), fsmType, p,
					countPos, sequence);
			countPos += matCounts;

			if (m == 0)
				machines[m].constructHash();
			else
				machines[m].setHash(machines[0].getHash());

			machines[m].setPlot(plot);
			machines[m].setPlotActive(plotActive);
			machines[m].setPlotHits(plotHits);

			// VNTRReadDepth.printf("Starting costs:\n%s:\n", fsmType.paramsToString());
		}
		Misc.my_assert(countPos == totCounts,
				"Internal error: countPos!=totCounts");

		// Setup the base states.
		base = new double[2];
		base[0] = 0;

		baseCounts = new Counts[2];
		baseCounts[0] = new Counts(totCounts);
		baseCounts[1] = new Counts(totCounts);

		last_msgLen = 0;
	}

	void set_default_costs() {
		// for (int i = 0; i < numFwd + numRev; i++)
		// encStartMachines[i] = -MyMath.log2(0.1 / (numFwd + numRev));
		// encNoStart = -MyMath.log2(0.9);

		// mdc changed to an good param
		encNoStart = 0.001602767095846256;// -MyMath.log2(0.9)

		for (int i = 0; i < numFwd + numRev; i++) {
			if (i == 0)
				encStartMachines[i] = 10.422286251724662;
			// else if (i == 1) encStartMachines[i] = 11.355628428955441;
			else
				encStartMachines[i] = 11.355705924947152;// -MyMath.log2(0.1 /
															// (numFwd +
															// numRev));
		}

	}

	void normalize_costs() {
		double sum = MyMath.exp2(-encNoStart);
		for (int i = 0; i < numFwd + numRev; i++)
			sum += MyMath.exp2(-encStartMachines[i]);

		encNoStart += MyMath.log2(sum);
		for (int i = 0; i < numFwd + numRev; i++)
			encStartMachines[i] += MyMath.log2(sum);

		if (DEBUG >= 2) {
			System.out.printf("nostart_cost=%f ", encNoStart);
			for (int i = 0; i < numFwd + numRev; i++) {
				System.out.printf("start" + i + "_cost=%f ",
						encStartMachines[i]);
				System.out.print(" == start" + i + "_cost="
						+ encStartMachines[i]);
			}

			System.out.println();
		}
	}

	public double encodeLen(char a, int i) throws RuntimeException {
		throw new RuntimeException("encodeLen not implemented here");
	}

	public double update(char aChar, int i) {
		// Start base[i+1] as base[i] with noStart and the char (also the
		// counts)
		base[1] = base[0] + seqModel.encodeLen(aChar, i)
				+ (i == 0 ? 0 : encNoStart);
		baseCounts[1].duplicate(baseCounts[0]);
		baseCounts[1].inc(noStartIndex, 1);

		for (int m = 0; m < numFwd + numRev; m++) {
			Counts tmpCounts = new Counts(totCounts);
			tmpCounts.duplicate(baseCounts[0]);
			tmpCounts.inc(noStartIndex + 1 + m, 1);

			machines[m].beginLinks(i, base[0] + encStartMachines[m]
					+ encStartPos(i), tmpCounts);
		}

		// update matches and calc their contribution to base
		double msgLen = base[1];
		for (int m = 0; m < numFwd + numRev; m++) {
			Counts retCounts = new Counts(totCounts);
			double len = machines[m].update(aChar, i, retCounts);

			baseCounts[1].combine_with_lens(base[1], retCounts, len);
			base[1] = MyMath.logplus(base[1], len);

			msgLen = MyMath.logplus(msgLen, machines[m].msgLen());
		}

		// Update the sequence model with the new character.
		seqModel.update(aChar, i);

		double p = MyMath.exp2(msgLen - base[1]);
		p = Math.pow(p, 0.3);
		plot.putMax(i + 1, i + 1, p, p, p);
		plotActive.putMax(i + 1, i + 1, p, p, p);
		plotHits.putMax(i + 1, i + 1, p, p, p);
		for (int m = 0; m < numFwd + numRev; m++) {
			machines[m].plotVals(i + 1, msgLen);
			if (DEBUG >= 4)
				((Matches_Sparse) machines[m]).display_stats();
		}

		if (msgLen > 100) { // Renormalise values so we don't lose accuracy
			double factor = last_msgLen;
			base[1] -= factor;

			msgLen = base[1];
			for (int m = 0; m < numFwd + numRev; m++) {
				msgLen = MyMath.logplus(msgLen, machines[m].normalise(factor));
			}
			last_msgLen -= factor;
			// System.err.println("Normalising. factor="+factor+"
			// last_msgLen="+last_msgLen);
		}

		base[0] = base[1];
		Counts t = baseCounts[0];
		baseCounts[0] = baseCounts[1];
		baseCounts[1] = t;

		double tm = msgLen;
		msgLen -= last_msgLen;
		last_msgLen = tm;

		return msgLen;
	}

	/**
	 * The start position of a copy is simply encoded from a uniform prior over
	 * possible start positions from [0..n]
	 */
	double encStartPos(int n) {
		return MyMath.log2(n + 1);
	}

	public Params counts_to_params() {
		Counts counts = baseCounts[0];

		Params res = new Params();

		double sum = counts.get(noStartIndex);
		// Prepend the parameters for the matches with 'MACHINEx_'
		for (int m = 0; m < numFwd + numRev; m++) {
			Params p = machines[m].counts_to_params(counts);
			for (int i = 0; i < p.get_num(); i++) {
				String n = p.get_name_by_id(i);
				res.put("MACHINE" + m + "_" + n, p.get(n));
			}

			sum += counts.get(noStartIndex + 1 + m);
		}

		// Finally parameters from us...
		res.put("nostart_cost", -MyMath.log2(counts.get(noStartIndex) / sum));
		for (int m = 0; m < numFwd + numRev; m++) {
			res.put("start" + m + "_cost",
					-MyMath.log2(counts.get(noStartIndex + 1 + m) / sum));
		}

		return res;
	}

	public void display() {
		// VNTRReadDepth.printf("Final counts:\n%s\n", baseCounts[seqLen]);
		System.out.printf("Final costs:\n%s\n", counts_to_params());
		// TODO: check this pls
	}

	public void display_stats() {
		for (int m = 0; m < numFwd + numRev; m++) {
			((Matches_Sparse) machines[m]).display_stats();
		}
		System.out.println("");
	}

	public void saveState(String fname) {
		try {
			File f = new File(fname);
			ObjectOutputStream oos = new ObjectOutputStream(
					new FileOutputStream(f));
			oos.writeObject(this);
			oos.close();
		} catch (Exception e) {
			System.err.println("Error writing file: " + e);
		}
	}

	public static FuzzyLZ readState(String fname) {
		FuzzyLZ r = null;
		Object o = null;
		try {
			FileInputStream istream = new FileInputStream(fname);
			ObjectInputStream p = new ObjectInputStream(istream);
			o = p.readObject();
			istream.close();
		} catch (Exception e) {
			System.err.println("Error reading file: " + e);
		}

		r = (FuzzyLZ) o;
		return r;
	}

}
