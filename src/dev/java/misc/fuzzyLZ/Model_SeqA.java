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

package misc.fuzzyLZ;

import misc.common.*;

/**
 * A model of 2 sequences (for alignments) with counts. Uses two parameters:
 * match_cost, change_cost
 * 
 * This is really a pseudo 2-sequence model. It assumes that sequence B is
 * already known to the receiver. So, it encodes a character from sequence A
 * using a cost 'match_cost' if it is the same as the character from B. If it is
 * different it is has a cost 'change_cost' plus a cost for the character using
 * a uniform model over the rest of the dna
 */

class Model_SeqA implements Two_Seq_Model_Counts {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	int alphaSize;

	double match_cost, change_cost;

	int countIndex;

	final private int matchIndex = 0, changeIndex = 1;

	public Model_SeqA(Params p, int alphaSize, int countIndex) {
		this.alphaSize = alphaSize;
		this.countIndex = countIndex;

		if (!p.exists("match_cost")) {
			set_default_costs();
		} else {
			match_cost = p.get("match_cost");
			change_cost = p.get("change_cost");
		}

		normalize_costs();
	}

	void set_default_costs() {
		// match_cost = -MyMath.log2(9.0);
		// change_cost = -MyMath.log2(1.0);
		// mdc changed to a good set of variable

		match_cost = 0.20554840023089077;//
		change_cost = 2.912770498714653;// -MyMath.log2(1.0);
	}

	void normalize_costs() {
		double sum = MyMath.exp2(-match_cost) + MyMath.exp2(-change_cost);
		match_cost = match_cost + MyMath.log2(sum);
		change_cost = change_cost + MyMath.log2(sum);

		if (FuzzyLZ.DEBUG >= 2) {
			System.out.printf("match_cost=%f change_cost=%f\n", match_cost,
					change_cost);
			System.out.println(" ==match_cost=" + match_cost + " change_cost="
					+ change_cost);
		}

	}

	public double encA(char a, int i) {
		return MyMath.log2(alphaSize);
	}

	public double encB(char a, int i) {
		return 0;
	}

	public double encBoth(char a, char b, int i, int j) {
		return ((a == b) ? match_cost : change_cost
				+ MyMath.log2(alphaSize - 1));
	}

	public static int required_counts() {
		return 2;
	}

	public Params counts_to_params(Counts counts) {
		double sum = counts.counts[countIndex + matchIndex]
				+ counts.counts[countIndex + changeIndex];
		Params par = new Params();
		par.put("match_cost",
				-MyMath.log2(counts.counts[countIndex + matchIndex] / sum));
		par.put("change_cost",
				-MyMath.log2(counts.counts[countIndex + changeIndex] / sum));
		return par;
	}

	public void update_count_encA(Counts c, double w, char a, int i) {
	};

	public void update_count_encB(Counts c, double w, char a, int i) {
	};

	public void update_count_encBoth(Counts c, double w, char a, char b, int i,
			int j) {
		if (a == b) {
			c.inc(countIndex + matchIndex, w);
		} else {
			c.inc(countIndex + changeIndex, w);
		}
	}

	public double encode_params(double N) {
		return 0;
	}
}
