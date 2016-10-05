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

import java.io.Serializable;

public final class Params implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String names[];
	private double vals[];
	private int num, room;

	private void alloc(int r) {
		if (room >= r)
			return;

		String[] n = new String[r];
		double[] v = new double[r];

		if (room > 0) { // Copy over any data
			for (int i = 0; i < num; i++) {
				n[i] = names[i];
				v[i] = vals[i];
			}
		}

		room = r;
		names = n;
		vals = v;
	}

	public Params() {
		num = 0;
		room = 0;
		alloc(10);
	}

	public Params put(String s, double v) {
		int id = get_id(s);
		if (id < 0) {
			if (num == room)
				alloc(room * 2);
			names[num] = new String(s);
			id = num;
			num++;
		}
		vals[id] = v;

		return this;
	}

	private int get_id(String s) {
		for (int i = 0; i < num; i++) {
			if (s.equalsIgnoreCase(names[i]))
				return i;
		}
		return -1;
	}

	public boolean exists(String s) {
		return get_id(s) >= 0;
	}

	public double get(String s) {
		int id = get_id(s);
		Misc.my_assert(id >= 0, "Attempt Params.get() with non-existent key");
		return vals[id];
	}

	public int get_num() {
		return num;
	}

	public String get_name_by_id(int id) {
		return names[id];
	}

	public void join(Params p) {
		for (int i = 0; i < p.num; i++) {
			put(p.names[i], p.vals[i]);
		}
	}

	public String toString() {
		StringBuffer r = new StringBuffer();
		for (int i = 0; i < num; i++) {
			r.append(names[i] + "=" + vals[i] + "\n");
		}
		return r.toString();
	}

	public void fromString(String str) {
		int s = 0;
		int e;
		do {
			int p_comma, p_newline;
			p_comma = str.indexOf(',', s);
			p_newline = str.indexOf('\n', s);

			if (p_comma < 0)
				p_comma = str.length();
			if (p_newline < 0)
				p_newline = str.length();

			e = (p_comma < p_newline ? p_comma : p_newline);

			int m = str.indexOf('=', s);
			if (s < m && m < e) {
				String p = str.substring(s, m);
				String v = str.substring(m + 1, e);
				double d = Double.parseDouble(v);

				put(p, d);
			}
			s = e + 1;
		} while (s < str.length());
	}

	public static void main(String args[]) {
		Params p = new Params();
		p.fromString(args[0]);
		System.out.println("p=" + p);
	}

}
