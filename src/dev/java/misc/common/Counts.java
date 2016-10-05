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

package misc.common;

public final class Counts {// implements Serializable {
	public double counts[];
	public int num;

	public Counts(int num) {
		counts = new double[num];
		this.num = num;
	}

	public double get(int i) {
		return counts[i];
	};

	public void zero() {
		for (int i = 0; i < num; i++) {
			counts[i] = 0;
		}
	}

	public void combine_with_lens(double myLen, Counts c, double otherLen) {
		double min = MyMath.min2(myLen, otherLen);
		double w1 = (min == myLen ? 1 : MyMath.exp2(min - myLen));
		double w2 = (min == otherLen ? 1 : MyMath.exp2(min - otherLen));
		scale(w1);
		linearWeight(c, w2);
		scale(1.0 / (w1 + w2));
	}

	public void linearWeight(Counts c, double w) {
		for (int i = 0; i < num; i++) {
			counts[i] += w * c.counts[i];
		}
	}

	public void inc(int index, double w) {
		counts[index] += w;
	}

	public void scale(double w) {
		for (int i = 0; i < num; i++)
			counts[i] *= w;
	}

	public void duplicate(Counts c) {
		for (int i = 0; i < num; i++)
			counts[i] = c.counts[i];
	}

	public Object clone() {
		Counts c = new Counts(num);
		c.duplicate(this);
		return c;
	}

	public String toString() {
		String r = new String("");
		for (int i = 0; i < num; i++) {
			r = Misc.sprintf("%s %d:%.3f", new Object[] { r, new Integer(i),
					new Double(counts[i]) });
		}
		return r;
	}
}
