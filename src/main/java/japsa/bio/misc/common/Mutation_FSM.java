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

import java.io.Serializable;



public abstract class Mutation_FSM implements Has_Value, Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 4340122370640374306L;
	Counts counts;
	int numCounts;
	int countIndex;

	class FSM_Params {
		Two_Seq_Model s;

		FSM_Params(Two_Seq_Model s) {
			this.s = s;
		};
	}

	public Mutation_FSM() {
		this(0, 0);
	}

	public Mutation_FSM(int numCounts, int countIndex) {
		this.numCounts = numCounts;
		counts = new Counts(numCounts);
		this.countIndex = countIndex;
	}

	public abstract void init_val(double v);

	public abstract double get_val();

	public abstract void normalise(double v);

	public abstract void calc(Mutation_FSM h, Mutation_FSM v, Mutation_FSM d,
			char a, char b, int i, int j);

	public void reset() {
		counts.zero();
	}

	/** Initialise all counts */
	public void init_counts(Counts c) {
		counts.duplicate(c);
	}

	public static int required_counts() {
		return 0;
	};

	public abstract Params counts_to_params(Counts c);

	public Counts get_counts() {
		return counts;
	};

	public double alignmentLength() {
		System.err.println("WARNING: alignmentLength() not implemented in "
				+ this.getClass());
		return 0;
	};

	public abstract void add(double d, int cIndex);

	public abstract void or(double d, Counts c);

	public String paramsToString() {
		return this.getClass() + ": paramsToString() not defined" + "\n";
	}

	public double encode_params() {
		System.err.println("WARNING: encode_params() not implemented in "
				+ this.getClass());
		return 0;
	}

	public abstract Object clone();

	public static class TraceBack_Data {
		public int i, j;

		public TraceBack_Data() {
			i = -1;
			j = -1;
		}

		public TraceBack_Data(int i, int j) {
			this.i = i;
			this.j = j;
		}

		public String toString() {
			return "i=" + i + " j=" + j;
		};
	}

	public interface TraceBack_Info {
		public void set_tbdata(TraceBack_Data id);

		public TraceBack_Data get_tbdata();

		public TraceBack_Data get_from(TraceBack_Data td);

		public void or(double d, Counts c, Mutation_FSM from);
	}
}
