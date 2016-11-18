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

import java.io.*;
import java.util.*;

public final class Misc {
	public static void my_assert(boolean e, String s) {
		if (!e)
			error(s);
	}

	public static void error(String s) {
		System.err.println("ERROR: " + s);

		// VNTRReadDepth o = (VNTRReadDepth) (new Object()); // Force a crash. want stack trace.
		throw new RuntimeException("Assertion wrong");
		// System.exit(1);
	}

	static String readString(InputStream in) {
		StringBuffer s = new StringBuffer();
		try {
			byte[] buf = new byte[1024];
			int n;
			while ((n = in.read(buf)) >= 0) {
				s.append(new String(buf, 0, n));
			}
		} catch (Exception e) {
			System.err.println("Unable to read stdin");
			System.exit(1);
		}
		return s.toString();
	}

	private static String pad(String s, int width, boolean pre0, boolean right) {
		if (s.length() >= width)
			return s;
		StringBuffer prefix = new StringBuffer();
		while (prefix.length() + s.length() < width)
			prefix.append((pre0 ? '0' : ' '));

		if (right || pre0)
			return prefix.toString() + s;
		else
			return s + prefix.toString();
	}

	static String numFmt(long num, int width, boolean pre0, boolean right) {
		return pad(Long.toString(num), width, pre0, right);
	}

	static String numFmt(float num, int prec, int width, boolean pre0,
			boolean right) {
		return numFmt((double) num, prec, width, pre0, right);
	}

	static String numFmt(double num, int prec, int width, boolean pre0,
			boolean right) {
		if (prec < 0)
			return pad(Double.toString(num), width, pre0, right);

		// Round num first.
		double mul = Math.pow(10, prec);
		num = Math.round(num * mul) / mul;

		long i = (long) num;
		StringBuffer res;
		if (i == 0 && num < 0)
			res = new StringBuffer("-0");
		else
			res = new StringBuffer(Long.toString(i));
		if (prec == 0)
			return pad(res.toString(), width, pre0, right);

		res.append(".");

		num = Math.abs(num - i);
		long dec = Math.round(num * Math.pow(10, prec));
		res.append(numFmt(dec, prec, true, false));
		return pad(res.toString(), width, pre0, right);
	}

	/**
	 * An extremely limited and poorly implemented sprintf function. Only
	 * handles %s %d %f with width and precision
	 */
	public static String sprintf(String fmt, Object[] objs) {
		StringBuffer res = new StringBuffer();
		int obj_i = 0;
		int i = 0;
		while (i < fmt.length()) {
			int p = fmt.indexOf('%', i);
			if (p < 0 || p + 1 == fmt.length()) {
				res.append(fmt.substring(i));
				break;
			}

			res.append(fmt.substring(i, p));

			if (fmt.charAt(p + 1) == '%') {
				res.append("%");
				i = p + 2;
				continue;
			}

			boolean rightAlign = true;
			boolean pre0 = false;
			int width = -1;
			int prec = -1;
			p++;
			while (true) {
				char c = fmt.charAt(p);
				if (c == 'd') {
					Object o = objs[obj_i++];
					if (o instanceof Integer) {
						res.append(numFmt(((Integer) o).intValue(), width,
								pre0, rightAlign));
					} else {
						System.err.println("Format %d not passed an integer");
						res.append(o);
					}
					p++;
					break;
				}
				if (c == 'f') {
					Object o = objs[obj_i++];
					if (o instanceof Double) {
						res.append(numFmt(((Double) o).doubleValue(), prec,
								width, pre0, rightAlign));
					} else if (o instanceof Float) {
						res.append(numFmt(((Float) o).floatValue(), prec,
								width, pre0, rightAlign));
					} else {
						System.err
								.println("Format %d not passed an double or float");
						res.append(o);
					}
					p++;
					break;
				}
				if (c == 's') {
					res.append(pad(objs[obj_i++].toString(), width, false,
							rightAlign));
					p++;
					break;
				}

				if (c == '.') {
					prec = -2;
					p++;
					continue;
				}

				if (c == '-') {
					rightAlign = false;
					p++;
					continue;
				}

				if (c >= '0' && c <= '9') {
					int j;
					for (j = p + 1; j < fmt.length() && c >= '0' && c <= '9'; j++) {
						c = fmt.charAt(j);
					}
					int v = Integer.parseInt(fmt.substring(p, j - 1));
					if (prec == -2)
						prec = v;
					else {
						width = v;
						pre0 = fmt.charAt(p) == '0';
					}
					p = j - 1;
					continue;
				}
				System.err.println("Unexpected character in % format '" + c
						+ "'");
				break;
			}
			i = p;
		}
		return res.toString();
	}

	/**
	 * Can use like this: sprintf( fmt, new
	 * VarArgs(i1).add(i2).add(i3).add(i4));
	 */
	public static String sprintf(String fmt, VarArgs v) {
		return sprintf(fmt, v.toArray());
	}

	// public static void printf(String fmt) {
	// System.out.print(fmt);
	// }

	public static void printf(String fmt, Object[] objs) {
		System.out.print(sprintf(fmt, objs));
	}

	// public static void printf(String fmt, VarArgs v) {
	// printf(fmt, v.toArray());
	// }

	// public static void printf(String fmt, Object v1) {
	// printf(fmt, new Object[] {v1});
	// }

	// public static void printf(String fmt, int v1) {
	// printf(fmt, new Object[] {new Integer(v1)});
	// }

	// public static void printf(String fmt, int v1, int v2) {
	// printf(fmt, new Object[] {new Integer(v1), new Integer(v2)});
	// }

	public static void printf(String fmt, int v1, int v2, int v3) {
		printf(fmt, new Object[] { new Integer(v1), new Integer(v2),
				new Integer(v3), });
	}

	public static void printf(String fmt, int v1, int v2, int v3, int v4) {
		printf(fmt, new Object[] { new Integer(v1), new Integer(v2),
				new Integer(v3), new Integer(v4) });
	}

	public static void printf(String fmt, double v1) {
		printf(fmt, new Object[] { new Double(v1) });
	}

	public static void printf(String fmt, double v1, double v2) {
		printf(fmt, new Object[] { new Double(v1), new Double(v2) });
	}

	public static void printf(String fmt, double v1, double v2, double v3) {
		printf(fmt, new Object[] { new Double(v1), new Double(v2),
				new Double(v3), });
	}

	public static void printf(String fmt, double v1, double v2, double v3,
			double v4) {
		printf(fmt, new Object[] { new Double(v1), new Double(v2),
				new Double(v3), new Double(v4) });
	}

	public static class VarArgs {
		Vector<Object> o = new Vector<Object>();

		public VarArgs(char v) {
			add(v);
		}

		public VarArgs(int v) {
			add(v);
		}

		public VarArgs(double v) {
			add(v);
		}

		public VarArgs(float v) {
			add(v);
		}

		public VarArgs(Object v) {
			add(v);
		}

		public VarArgs add(char v) {
			return add(new Character(v));
		}

		public VarArgs add(int v) {
			return add(new Integer(v));
		}

		public VarArgs add(double v) {
			return add(new Double(v));
		}

		public VarArgs add(float v) {
			return add(new Float(v));
		}

		public VarArgs add(Object v) {
			o.add(v);
			return this;
		}

		public Object[] toArray() {
			return o.toArray();
		}
	}

	public static void main(String args[]) {
		System.out.println(sprintf(args[0], new Object[] { new Double(5),
				new String("abc") }));
		// System.out.println( numFmt(10, 5, true) );
		// System.out.println( numFmt(10.1099, -2) );
		// System.out.println( numFmt(-10.9999, 3) );
	}

}
