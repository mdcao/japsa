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
 * An arbitrary order, adaptive Markov Model for an arbitrary dna of
 * characters
 * 
 * If order == -1 then we are a uniform model over the dna
 */
public class MarkovN implements Seq_Model {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	char[] chars;
	int[] charCounts;
	int[] countTotal;
	int order;

	StringBuffer past;

	public MarkovN(int order, char[] chars) {
		Misc.my_assert(order >= -1, "Bad order=" + order);

		this.chars = chars;
		this.order = order;

		if (order >= 0) {
			charCounts = new int[(int) Math.pow(chars.length, order + 1)];
			countTotal = new int[(int) Math.pow(chars.length, order)];
			for (int i = 0; i < charCounts.length; i++)
				charCounts[i] = 1;
			for (int i = 0; i < countTotal.length; i++)
				countTotal[i] = chars.length;
		}

		past = new StringBuffer();
	}

	int chars2Num(String c) {
		int res = 0, i;
		for (i = 0; i < c.length(); i++) {
			int j;
			for (j = 0; j < chars.length; j++)
				if (chars[j] == c.charAt(i))
					break;
			Misc.my_assert(j < chars.length, "Character '" + c.charAt(i)
					+ "' is unexpected");
			res = (res * chars.length) + j;
		}
		// System.out.println("char2Num("+c+")="+res);
		return res;
	}

	public double encodeLen(char a, int i) {
		if (past.length() < order || order < 0)
			return -MyMath.log2((double) 1.0 / chars.length);
		int n = charCounts[chars2Num(past.substring(i - order) + a)];
		int d = countTotal[chars2Num(past.substring(i - order))];
		return -MyMath.log2((double) n / d);
	}

	public double update(char a, int i) {
		if (past.length() < order || order < 0) {
			past.append(a);
			return -MyMath.log2((double) 1.0 / chars.length);
		}
		double res = encodeLen(a, i);
		charCounts[chars2Num(past.substring(i - order) + a)]++;
		countTotal[chars2Num(past.substring(i - order))]++;
		past.append(a);
		return res;
	}

	/**
	 * Minh Duc Added
	 */
	public double probability(char a, int i) {
		if (past.length() < order || order < 0)
			return 1.0 / chars.length;
		int n = charCounts[chars2Num(past.substring(i - order) + a)];
		int d = countTotal[chars2Num(past.substring(i - order))];
		return ((double) n) / d;
	}

	public static void main(String args[]) {
		String s = args[0];
		char[] a = { 'a', 't', 'g', 'c' };
		MarkovN m = new MarkovN(0, a);

		double tot = 0;
		for (int i = 0; i < s.length(); i++) {
			double r = m.update(s.charAt(i), i);
			tot += r;
			System.out.println(r);
		}
		System.out.println("Total entropy = " + tot + " bits/ch");
	}
}
