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

package japsa.bio.misc.common;

import java.io.*;

/**
 * Despite its name, this class does _not_ implement a full 3 state FSM. It
 * implements a 3 state FSM for the purpose of Linear gapped costs for a DPA.
 * There are therefore 5 parameters: given the last operation was a diagonal
 * (match or change): diag, start gap given the last operation was a gap: diag,
 * start new gap, continue current gap
 * 
 * So, the start_fromD encodes the start of a gap _and_ the first character
 * 2**-diag_fromD + 2 * 2**-start_fromD = 1 (2 possible ways to start a gap: ins
 * or del) 2**-diag_fromI + 2**-start_fromI + 2**-cont_fromI = 1
 * 
 * The match/mismatch costs are part of the Two_Seq_Model
 */

public abstract class Mutation_3State extends Mutation_FSM {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private class FSM_Params implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		Two_Seq_Model_Counts s;

		double diag_fromD, start_fromD;

		double diag_fromI, start_fromI, cont_fromI;

		/*
		 * FSM_Params(Two_Seq_Model_Counts s, double diag_cost, double
		 * start_gap, double cont_gap) { this.s = s;
		 * 
		 * normalize_costs(diag_cost, start_gap, cont_gap); }
		 */

		FSM_Params(Two_Seq_Model_Counts s, Params p) {
			this.s = s;

			if (!p.exists("diag_fromD")) {
				set_default_costs();
			} else {
				diag_fromD = p.get("diag_fromD");
				start_fromD = p.get("start_fromD");
				diag_fromI = p.get("diag_fromI");
				start_fromI = p.get("start_fromI");
				cont_fromI = p.get("cont_fromI");

				normalize_costs();
			}
		}

		void set_default_costs() {
			/*
			 * diag_fromD = 0; start_fromD = 3;
			 * 
			 * diag_fromI = 0; start_fromI = 3; cont_fromI = 1;
			 */

			// These costs are equivalent to SW scores (provided local char cost
			// = 0.104808)
			// Provided we are using DNA at 2 bits per char
			diag_fromD = MyMath.logplus(0.709616, 1.824653);
			start_fromD = 4.235037;

			diag_fromI = MyMath.logplus(1.379387, 2.494424);
			start_fromI = 4.904808;
			cont_fromI = 1.304808;

			if (s instanceof Model_SeqAB) {
				((Model_SeqAB) s).match_cost = 0.709616;
				((Model_SeqAB) s).change_cost = 1.824653;
				((Model_SeqAB) s).normalize_costs();
			}

			normalize_costs();
		}

		void normalize_costs() {
			double sum;

			sum = MyMath.exp2(-diag_fromD) + 2 * MyMath.exp2(-start_fromD);
			diag_fromD = diag_fromD + MyMath.log2(sum);
			start_fromD = start_fromD + MyMath.log2(sum);

			sum = MyMath.exp2(-diag_fromI) + MyMath.exp2(-start_fromI)
					+ MyMath.exp2(-cont_fromI);
			diag_fromI = diag_fromI + MyMath.log2(sum);
			start_fromI = start_fromI + MyMath.log2(sum);
			cont_fromI = cont_fromI + MyMath.log2(sum);
		}

		public String toString() {
			return "diag_fromD=" + diag_fromD + " start_fromD=" + start_fromD
					+ "\n" + "diag_fromI=" + diag_fromI + " start_fromI="
					+ start_fromI + " cont_fromI=" + cont_fromI;
		}
	}

	FSM_Params p;

	double dval, hval, vval;
	Counts d_counts, h_counts, v_counts;

	final protected int diag_fromD = 0, start_fromD = 1, diag_fromI = 2,
			start_fromI = 3, cont_fromI = 4;

	public Mutation_3State(Two_Seq_Model_Counts s, Params par, int numCounts,
			int countIndex) {
		super(numCounts, countIndex);
		p = new FSM_Params(s, par);
		reset();
	}

	public Mutation_3State(FSM_Params p, int numCounts, int countIndex) {
		super(numCounts, countIndex);
		this.p = p;
		reset();
	}

	public void reset() {
		super.reset();
		dval = MyMath.Big_Double;
		hval = MyMath.Big_Double;
		vval = MyMath.Big_Double;

		if (d_counts == null)
			d_counts = new Counts(numCounts);
		if (h_counts == null)
			h_counts = new Counts(numCounts);
		if (v_counts == null)
			v_counts = new Counts(numCounts);
		d_counts.zero();
		h_counts.zero();
		v_counts.zero();

		// counts = null;
	}

	public static int required_counts() {
		return 5;
	};

	// counts_to_params - convert counts into parameters.
	// Call the Two_Seq_Model to convert its own.
	// Note the 0.5* in the start_fromD computation. This is cause there are 2
	// start_fromD
	// arcs out of each state, but we only have one count for them combined.
	public Params counts_to_params(Counts c) {
		Params par = new Params();
		double sum1 = c.counts[countIndex + diag_fromD]
				+ c.counts[countIndex + start_fromD];
		double sum2 = c.counts[countIndex + diag_fromI]
				+ c.counts[countIndex + start_fromI]
				+ c.counts[countIndex + cont_fromI];

		par.put("diag_fromD",
				-MyMath.log2(c.counts[countIndex + diag_fromD] / sum1));
		par.put("start_fromD",
				-MyMath.log2(0.5 * c.counts[countIndex + start_fromD] / sum1));
		par.put("diag_fromI",
				-MyMath.log2(c.counts[countIndex + diag_fromI] / sum2));
		par.put("start_fromI",
				-MyMath.log2(c.counts[countIndex + start_fromI] / sum2));
		par.put("cont_fromI",
				-MyMath.log2(c.counts[countIndex + cont_fromI] / sum2));
		par.join(p.s.counts_to_params(c));

		return par;
	};

	public double encode_params() {
		Counts c = get_counts();
		double[] probs;
		double dataLen;

		// First encode match/change paramters. Pass the number of these events
		// to encode_params
		double len = p.s.encode_params(c.counts[countIndex + diag_fromD]
				+ c.counts[countIndex + diag_fromI]);

		// 'probs' is the paramaters of the multinomial distribution.
		// 'dataLen' is the number of things in this multinomial.

		probs = new double[] { MyMath.exp2(-p.diag_fromD),
				2 * MyMath.exp2(-p.start_fromD) };
		dataLen = c.counts[countIndex + diag_fromD]
				+ c.counts[countIndex + start_fromD];

		len += Multinomial.MMLparameter_cost(probs, dataLen);

		probs = new double[] { MyMath.exp2(-p.diag_fromI),
				MyMath.exp2(-p.start_fromI), MyMath.exp2(-p.cont_fromI) };
		dataLen = c.counts[countIndex + diag_fromI]
				+ c.counts[countIndex + start_fromI]
				+ c.counts[countIndex + start_fromI];

		len += Multinomial.MMLparameter_cost(probs, dataLen);

		return len;
	}

	public double alignmentLength() {
		Counts c = get_counts();
		return c.counts[countIndex + diag_fromD]
				+ c.counts[countIndex + start_fromD]
				+ c.counts[countIndex + diag_fromI]
				+ c.counts[countIndex + start_fromI]
				+ c.counts[countIndex + cont_fromI];
	}

	public void init_val(double v) {
		dval = v;
		hval = vval = MyMath.Big_Double;
	};

	public void normalise(double v) {
		dval -= v;
		hval -= v;
		vval -= v;
	}

	/** Initialise diag counts */
	public void init_counts(Counts c) {
		d_counts.duplicate(c);
	}

	public String paramsToString() {
		return this.getClass() + ": " + p + "\n";
	}

	public String toString() {
		return "dval=" + dval + " hval=" + hval + " vval=" + vval + "\n"
				+ "d_counts=" + d_counts + "\n" + "h_counts=" + h_counts + "\n"
				+ "v_counts=" + v_counts + "\n";
	}

	/**
	 * add() - something extra must be encoded in this state. Update all three
	 * states
	 */
	public void add(double v, int cIndex) {
		dval += v;
		hval += v;
		vval += v;
		if (cIndex >= 0) {
			d_counts.inc(cIndex, 1);
			h_counts.inc(cIndex, 1);
			v_counts.inc(cIndex, 1);
		}
	}

	public static class One extends Mutation_3State implements TraceBack_Info {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		private static class TB_Data extends TraceBack_Data {
			int state;

			// TB_Data() {
			// super();
			// state = -1;
			// }

			TB_Data(TraceBack_Data td, int s) {
				super();
				i = td.i;
				j = td.j;
				state = s;
			}
		}

		TB_Data d_id, v_id, h_id;
		TB_Data d_from, v_from, h_from;

		public One(FSM_Params p, int numCounts, int countIndex) {
			super(p, numCounts, countIndex);
			d_id = v_id = h_id = null;
			d_from = v_from = h_from = null;
		}

		public One(Two_Seq_Model_Counts s, Params par, int numCounts,
				int countIndex) {
			super(s, par, numCounts, countIndex);
			d_id = v_id = h_id = null;
			d_from = v_from = h_from = null;
		}

		public void set_tbdata(TraceBack_Data id) {
			d_id = new TB_Data(id, 0);
			v_id = new TB_Data(id, 1);
			h_id = new TB_Data(id, 2);
		}

		public TraceBack_Data get_tbdata() {
			if (dval <= hval && dval <= vval) {
				return d_id;
			} else if (vval <= dval && vval <= hval) {
				return v_id;
			} else if (hval <= dval && hval <= vval) {
				return h_id;
			}
			Misc.my_assert(false, "Bad vals");
			return null;
		}

		public TraceBack_Data get_from(TraceBack_Data td) {
			TB_Data t = (TB_Data) td;
			Misc.my_assert(t.i == d_id.i && t.j == d_id.j,
					"Not my traceback data!");
			switch (t.state) {
			case 0:
				return d_from;
			case 1:
				return v_from;
			case 2:
				return h_from;
			default:
				Misc.my_assert(false, "Bad trackback data. t.state=" + t.state);
			}
			return null;
		}

		public Object clone() {
			return new Mutation_3State.One(p, numCounts, countIndex);
		}

		public double get_val() {
			return MyMath.min3(dval, vval, hval);
		};

		public Counts get_counts() {
			// Get counts from state with smallest val
			if (dval <= hval && dval <= vval)
				return d_counts;

			if (vval <= dval && vval <= hval)
				return v_counts;

			if (hval <= dval && hval <= vval)
				return h_counts;

			Misc.my_assert(false, "Bad vals!");
			return null;
		};

		public void calc(Mutation_FSM h, Mutation_FSM v, Mutation_FSM d,
				char a, char b, int i, int j) {
			Mutation_3State.One hcell = (Mutation_3State.One) h;
			Mutation_3State.One vcell = (Mutation_3State.One) v;
			Mutation_3State.One dcell = (Mutation_3State.One) d;

			Counts tcounts = (Counts) counts.clone();
			if (hcell != null) {
				double char_cost = p.s.encB(b, j);
				double w;
				// double count_inc;

				// From 'd' state
				w = dval + char_cost + p.start_fromD;
				tcounts.duplicate(d_counts);
				tcounts.inc(countIndex + start_fromD, 1);
				p.s.update_count_encB(tcounts, 1, b, j);
				hcell.or_h(w, tcounts, d_id);

				// From 'v' state
				w = vval + char_cost + p.start_fromI;
				tcounts.duplicate(v_counts);
				tcounts.inc(countIndex + start_fromI, 1);
				p.s.update_count_encB(tcounts, 1, b, j);
				hcell.or_h(w, tcounts, v_id);

				// From 'h' state
				w = hval + char_cost + p.cont_fromI;
				tcounts.duplicate(h_counts);
				tcounts.inc(countIndex + cont_fromI, 1);
				p.s.update_count_encB(tcounts, 1, b, j);
				hcell.or_h(w, tcounts, h_id);
			}

			if (vcell != null) {
				double char_cost = p.s.encA(a, i);
				double w;
				// double count_inc;

				// From 'd' state
				w = dval + char_cost + p.start_fromD;
				tcounts.duplicate(d_counts);
				tcounts.inc(countIndex + start_fromD, 1);
				p.s.update_count_encA(tcounts, 1, a, i);
				vcell.or_v(w, tcounts, d_id);

				// From 'v' state
				w = vval + char_cost + p.cont_fromI;
				tcounts.duplicate(v_counts);
				tcounts.inc(countIndex + cont_fromI, 1);
				p.s.update_count_encA(tcounts, 1, a, i);
				vcell.or_v(w, tcounts, v_id);

				// From 'h' state
				w = hval + char_cost + p.start_fromI;
				tcounts.duplicate(h_counts);
				tcounts.inc(countIndex + start_fromI, 1);
				p.s.update_count_encA(tcounts, 1, a, i);
				vcell.or_v(w, tcounts, h_id);
			}

			if (dcell != null) {
				double char_cost = p.s.encBoth(a, b, i, j);
				double w;
				// double count_inc;

				// From 'd' state
				w = dval + char_cost + p.diag_fromD;
				tcounts.duplicate(d_counts);
				tcounts.inc(countIndex + diag_fromD, 1);
				p.s.update_count_encBoth(tcounts, 1, a, b, i, j);
				dcell.or_d(w, tcounts, d_id);

				// From 'v' state
				w = vval + char_cost + p.diag_fromI;
				tcounts.duplicate(v_counts);
				tcounts.inc(countIndex + diag_fromI, 1);
				p.s.update_count_encBoth(tcounts, 1, a, b, i, j);
				dcell.or_d(w, tcounts, v_id);

				// From 'h' state
				w = hval + char_cost + p.diag_fromI;
				tcounts.duplicate(h_counts);
				tcounts.inc(countIndex + diag_fromI, 1);
				p.s.update_count_encBoth(tcounts, 1, a, b, i, j);
				dcell.or_d(w, tcounts, h_id);
			}
		}

		public void or(double d, Counts c) {
			or(d, c, null);
		}

		/**
		 * or() - a new transition into this cell. All new transitions start in
		 * the diag state, so just call or_d
		 */
		public void or(double d, Counts c, Mutation_FSM from) {
			TB_Data f = null;
			if (from != null)
				f = (TB_Data) (((One) from).get_tbdata());
			or_d(d, c, f);
		}

		public void or_d(double d, Counts c, TB_Data from_id) {
			if (d < dval) {
				d_counts.duplicate(c);
				dval = d;
				d_from = from_id;
			}
		}

		public void or_h(double d, Counts c, TB_Data from_id) {
			if (d < hval) {
				h_counts.duplicate(c);
				hval = d;
				h_from = from_id;
			}
		}

		public void or_v(double d, Counts c, TB_Data from_id) {
			if (d < vval) {
				v_counts.duplicate(c);
				vval = d;
				v_from = from_id;
			}
		}
	}

	public static class All extends Mutation_3State {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public All(FSM_Params p, int numCounts, int countIndex) {
			super(p, numCounts, countIndex);
		}

		public All(Two_Seq_Model_Counts s, Params par, int numCounts,
				int countIndex) {
			super(s, par, numCounts, countIndex);
		}

		public Object clone() {
			return new Mutation_3State.All(p, numCounts, countIndex);
		}

		public double get_val() {
			double v = MyMath.logplus(dval, hval);
			v = MyMath.logplus(v, vval);
			return v;
		};

		public Counts get_counts() {
			// Combine counts from each of the three states weighted by their
			// val
			Counts c = (Counts) d_counts.clone();
			double val = dval;

			c.combine_with_lens(val, h_counts, hval);
			val = MyMath.logplus(val, hval);

			c.combine_with_lens(val, v_counts, vval);

			return c;
		};

		public void calc(Mutation_FSM h, Mutation_FSM v, Mutation_FSM d,
				char a, char b, int i, int j) {
			Mutation_3State.All hcell = (Mutation_3State.All) h;
			Mutation_3State.All vcell = (Mutation_3State.All) v;
			Mutation_3State.All dcell = (Mutation_3State.All) d;

			if (hcell != null) {
				double char_cost = p.s.encB(b, j);
				double w;
				double count_inc;

				// From 'd' state
				w = dval + char_cost + p.start_fromD;
				count_inc = 1.0 / (1.0 + MyMath.exp2(w - hcell.hval));
				if (Double.isNaN(count_inc))
					count_inc = 0;
				hcell.or_h(w, d_counts);
				hcell.p.s.update_count_encB(hcell.h_counts, count_inc, b, j);
				hcell.h_counts.inc(countIndex + start_fromD, count_inc);

				// From 'v' state
				w = vval + char_cost + p.start_fromI;
				count_inc = 1.0 / (1.0 + MyMath.exp2(w - hcell.hval));
				if (Double.isNaN(count_inc))
					count_inc = 0;
				hcell.or_h(w, v_counts);
				hcell.p.s.update_count_encB(hcell.h_counts, count_inc, b, j);
				hcell.h_counts.inc(countIndex + start_fromI, count_inc);

				// From 'h' state
				w = hval + char_cost + p.cont_fromI;
				count_inc = 1.0 / (1.0 + MyMath.exp2(w - hcell.hval));
				if (Double.isNaN(count_inc))
					count_inc = 0;
				hcell.or_h(w, h_counts);
				hcell.p.s.update_count_encB(hcell.h_counts, count_inc, b, j);
				hcell.h_counts.inc(countIndex + cont_fromI, count_inc);
			}

			if (vcell != null) {
				double char_cost = p.s.encA(a, i);
				double w;
				double count_inc;

				// From 'd' state
				w = dval + char_cost + p.start_fromD;
				count_inc = 1.0 / (1.0 + MyMath.exp2(w - vcell.vval));
				if (Double.isNaN(count_inc))
					count_inc = 0;
				vcell.or_v(w, d_counts);
				vcell.p.s.update_count_encA(vcell.v_counts, count_inc, a, i);
				vcell.v_counts.inc(countIndex + start_fromD, count_inc);

				// From 'v' state
				w = vval + char_cost + p.cont_fromI;
				count_inc = 1.0 / (1.0 + MyMath.exp2(w - vcell.vval));
				if (Double.isNaN(count_inc))
					count_inc = 0;
				vcell.or_v(w, v_counts);
				vcell.p.s.update_count_encA(vcell.v_counts, count_inc, a, i);
				vcell.v_counts.inc(countIndex + cont_fromI, count_inc);

				// From 'h' state
				w = hval + char_cost + p.start_fromI;
				count_inc = 1.0 / (1.0 + MyMath.exp2(w - vcell.vval));
				if (Double.isNaN(count_inc))
					count_inc = 0;
				vcell.or_v(w, h_counts);
				vcell.p.s.update_count_encA(vcell.v_counts, count_inc, a, i);
				vcell.v_counts.inc(countIndex + start_fromI, count_inc);
			}

			if (dcell != null) {
				double char_cost = p.s.encBoth(a, b, i, j);
				double w;
				double count_inc;

				// From 'd' state
				w = dval + char_cost + p.diag_fromD;
				count_inc = 1.0 / (1.0 + MyMath.exp2(w - dcell.dval));
				if (Double.isNaN(count_inc))
					count_inc = 0;
				dcell.or_d(w, d_counts);
				dcell.p.s.update_count_encBoth(dcell.d_counts, count_inc, a, b,
						i, j);
				dcell.d_counts.inc(countIndex + diag_fromD, count_inc);

				// From 'v' state
				w = vval + char_cost + p.diag_fromI;
				count_inc = 1.0 / (1.0 + MyMath.exp2(w - dcell.dval));
				if (Double.isNaN(count_inc))
					count_inc = 0;
				dcell.or_d(w, v_counts);
				dcell.p.s.update_count_encBoth(dcell.d_counts, count_inc, a, b,
						i, j);
				dcell.d_counts.inc(countIndex + diag_fromI, count_inc);

				// From 'h' state
				w = hval + char_cost + p.diag_fromI;
				count_inc = 1.0 / (1.0 + MyMath.exp2(w - dcell.dval));
				if (Double.isNaN(count_inc))
					count_inc = 0;
				dcell.or_d(w, h_counts);
				dcell.p.s.update_count_encBoth(dcell.d_counts, count_inc, a, b,
						i, j);
				dcell.d_counts.inc(countIndex + diag_fromI, count_inc);
			}

		}

		/**
		 * or() - a new transition into this cell. All new transitions start in
		 * the diag state, so just call or_d
		 */
		public void or(double d, Counts c) {
			or_d(d, c);

			/*
			 * if (d<dval) { d_counts.duplicate(c); dval = d; hval = vval =
			 * Double.POSITIVE_INFINITY; }
			 */
		}

		public void or_d(double d, Counts c) {
			// Update the counts first
			d_counts.combine_with_lens(dval, c, d);
			dval = MyMath.logplus(dval, d);
		}

		public void or_h(double d, Counts c) {
			// Update the counts first
			h_counts.combine_with_lens(hval, c, d);
			hval = MyMath.logplus(hval, d);
		}

		public void or_v(double d, Counts c) {
			// Update the counts first
			v_counts.combine_with_lens(vval, c, d);
			vval = MyMath.logplus(vval, d);
		}
	}
}
