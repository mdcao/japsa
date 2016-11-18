/*  
 *  Copyright (c) David Powell <david@drp.id.au>
 *
 * 
 * This file is used by both FuzzyLZ and AlignCompress

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

package japsadev.misc.common;

/**
 * A model of 2 sequences (for alignments) with counts.
 * 
 * Characters from each sequence are encoded with a sequence specific model.
 * Matches/changes average the probabilites from these models.
 * 
 * Uses two parameters: match_cost, change_cost
 **/
public class Model_SeqAB implements Two_Seq_Model_Counts {
	double match_cost, change_cost;
	Seq_Model modelA, modelB;

	int countIndex;
	final private int matchIndex = 0, changeIndex = 1;

	public Model_SeqAB(Params p, Seq_Model modelA, Seq_Model modelB,
			int countIndex) {
		this.countIndex = countIndex;
		this.modelA = modelA;
		this.modelB = modelB;

		if (!p.exists("match_cost")) {
			set_default_costs();
		} else {
			match_cost = p.get("match_cost");
			change_cost = p.get("change_cost");
		}

		normalize_costs();
	}

	public String toString() {
		return this.getClass() + ": match_cost=" + match_cost + " change_cost="
				+ change_cost;
	}

	void set_default_costs() {
		match_cost = -MyMath.log2(9.0);
		change_cost = -MyMath.log2(1.0);
	}

	void normalize_costs() {
		double sum = MyMath.exp2(-match_cost) + MyMath.exp2(-change_cost);
		match_cost = match_cost + MyMath.log2(sum);
		change_cost = change_cost + MyMath.log2(sum);
	}

	public double encA(char a, int i) {
		return modelA.encodeLen(a, i);
	}

	public double encB(char a, int i) {
		return modelB.encodeLen(a, i);
	}

	public double encBoth(char a, char b, int i, int j) {
		double A_cost = encA(a, i);
		double B_cost = encB(b, j);
		if (a == b) {
			// Match
			// Do: P(match) * ( P(char a) + P(char b) ) / 2
			// System.err.println("enc match = " + ( match_cost +
			// MyMath.logplus(A_cost, B_cost) + 1));
			return match_cost + MyMath.logplus(A_cost, B_cost) + 1;
		} else {
			// Change
			// Do: P(change) * P(char a) * P(char b) * 0.5 * (1/(1-P(char b)) +
			// 1/(1-P(char a)))
			double aN = MyMath.exp2(-encA(b, i));
			double bN = MyMath.exp2(-encB(a, j));
			double norm = -MyMath.log2(1 / (1 - aN) + 1 / (1 - bN));
			// System.err.println("enc change = " + (change_cost + A_cost +
			// B_cost + 1 + norm));
			return change_cost + A_cost + B_cost + 1 + norm;
		}
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
		return Multinomial.MMLparameter_cost(
				new double[] { MyMath.exp2(-match_cost),
						MyMath.exp2(-change_cost) }, N);
	}
}
