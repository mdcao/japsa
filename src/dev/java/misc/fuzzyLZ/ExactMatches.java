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

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;

/**
 * This class implements a hashing scheme to find matches of a fixed size within
 * a String
 * 
 * @author David Powell
 * @version 24/5/2001
 */

public class ExactMatches implements Serializable {
	private static final long serialVersionUID = 1L;

	// public static class Convert implements Serializable {
	public static class Convert implements Serializable {
		private static final long serialVersionUID = 1L;

		public String conv(String s) {
			return s;
		}

		public char conv(char a) {
			return a;
		}
	}

	public static class Reverse_Complement_DNA extends Convert {
		private static final long serialVersionUID = 1L;

		public String conv(String s) {
			// return (new StringBuffer(s)).reverse().toString();
			char c[] = new char[s.length()];
			for (int i = 0; i < s.length(); i++) {
				char a = s.charAt(i);
				a = conv(a);
				c[s.length() - i - 1] = a;
			}
			return new String(c);
		};

		public char conv(char a) {
			switch (a) {
			case 'A':
				a = 'T';
				break;
			case 'T':
				a = 'A';
				break;
			case 'G':
				a = 'C';
				break;
			case 'C':
				a = 'G';
				break;
			case 'a':
				a = 't';
				break;
			case 't':
				a = 'a';
				break;
			case 'g':
				a = 'c';
				break;
			case 'c':
				a = 'g';
				break;
			default:
				System.err.println("WARNING: unknown DNA character " + a
						+ " in Reverse_Complement_DNA (char unchanged)");
				return a;
			}
			return a;
		}
	}

	// Internal representation of a linked-list of ints.
	public static class MyList implements Serializable {
		private static final long serialVersionUID = 1L;

		public static class L implements Serializable {
			private static final long serialVersionUID = 1L;
			L next;

			int val;

			L(int v) {
				next = null;
				val = v;
			}
		}

		L start, end;

		MyList(int v) {
			L n = new L(v);
			start = end = n;
		}

		void add(int v) {
			L n = new L(v);
			end.next = n;
			end = n;
		}
	}

	private static class MyHash implements Serializable {
		private static final long serialVersionUID = 1L;

		private static class HashChain implements Serializable {
			private static final long serialVersionUID = 1L;
			HashChain next;

			MyList intList;
		}

		int hSize = 1048573;

		HashChain hTable[];

		char[] str;

		int winSize;

		MyHash(char[] str, int winSize) {
			this.str = str;
			this.winSize = winSize;
			hTable = new HashChain[hSize];
		}

		private int hash(char[] s, int p) {
			int res = 0;
			for (int i = 0; i < winSize; i++) {
				res = (res << 1) + res + s[p + i];
			}
			return Math.abs(res) % hSize;
		}

		private boolean equal(char[] s1, int p1, char[] s2, int p2) {
			for (int i = 0; i < winSize; i++)
				if (s1[p1 + i] != s2[p2 + i])
					return false;
			return true;
		}

		MyList get(String s) {
			return get(s.toCharArray(), 0);
		}

		MyList get(char[] s, int ipos) {
			int i = hash(s, ipos);
			HashChain c = hTable[i];
			// Check if any of the strings at this hash bucket match the string
			// s
			while (c != null) {
				int pos = c.intList.start.val; // Get the string posSrc of the
				// first.
				if (equal(s, ipos, str, pos)) {
					// A match!
					return c.intList;
				}
				c = c.next;
			}
			return null;
		}

		void put(int pos) {
			MyList l = get(str, pos);
			if (l == null) {
				// First occurance of this string
				HashChain c = new HashChain();
				c.intList = new MyList(pos);

				int i = hash(str, pos);
				c.next = hTable[i];
				hTable[i] = c;
			} else {
				l.add(pos);
			}

		}
	}

	char[] str;

	int winSize;

	int strLen;

	MyHash h;

	ExactMatches(char[] str, int winSize) {
		this.str = str;
		this.winSize = winSize;

		strLen = str.length;

		if (FuzzyLZ.DEBUG >= 2)
			System.err.println("Constructing table of repeats...");

		h = new MyHash(str, winSize);
		for (int i = 0; i < strLen + 1 - winSize; i++) {
			h.put(i);
		}

		if (FuzzyLZ.DEBUG >= 2)
			System.err.println("Done constructing table of repeats.");
	}

	// Write our own serization handler. We will _not_ save the hash table.
	// Only the string, and recompute the hash table on reload
	private void writeObject(java.io.ObjectOutputStream out) throws IOException {
		out.writeObject(str);
		out.writeInt(winSize);
		out.writeInt(strLen);
	}

	private void readObject(java.io.ObjectInputStream in) throws IOException,
			ClassNotFoundException {
		str = (char[]) in.readObject();
		winSize = in.readInt();
		strLen = in.readInt();

		if (FuzzyLZ.DEBUG >= 2)
			System.err.println("Re-Constructing table of repeats...");

		h = new MyHash(str, winSize);
		for (int i = 0; i < strLen + 1 - winSize; i++) {
			h.put(i);
		}

		if (FuzzyLZ.DEBUG >= 2)
			System.err.println("Done re-constructing table of repeats.");
	}

	MyList get(char[] s, int pos) {
		return h.get(s, pos);
	}

	MyList get(String s) {
		return h.get(s);
	}

	public long count_hits(Convert convert) {
		long count = 0;
		for (int i = 0; i < strLen + 1 - winSize; i++) {
			String s = new String(str, i, winSize);
			MyList l = h.get(convert.conv(s));
			if (l != null) {
				for (MyList.L l2 = l.start; l2 != null && l2.val < i; l2 = l2.next) {
					count++;
				}
			}
		}
		return count;
	}

	public void dispAll(Convert convert) {
		for (int i = 0; i < strLen + 1 - winSize; i++) {
			String s = new String(str, i, winSize);
			MyList l = h.get(convert.conv(s));
			if (l != null) {
				for (MyList.L l2 = l.start; l2 != null && l2.val < i; l2 = l2.next) {
					System.out.println("String at " + i + " has a match at "
							+ l2.val);
					System.out.println("    :" + s + ":"
							+ new String(str, l2.val, winSize));
				}
			}
		}
	}

	@SuppressWarnings("unused")
	public static void main(String args[]) {
		int len = Integer.valueOf(args[0]).intValue();
		char[] input = null;

		try {
			StringBuffer s = new StringBuffer();
			byte[] buf = new byte[1024];
			int n;
			while ((n = System.in.read(buf)) >= 0) {
				s.append(new String(buf, 0, n));
			}
			input = s.toString().toCharArray();
		} catch (Exception e) {
			System.err.println("Unable to read stdin");
			System.exit(1);
		}

		System.out.println("String length = " + input.length);

		ExactMatches m = new ExactMatches(input, len);
		// System.out.println("Number of matches = "+m.count());

		if (true) {
			// m.dispAll(new Convert());
			m.dispAll(new Reverse_Complement_DNA());
		} else {
			try {
				File f = new File("ExactMatches.store");
				ObjectOutputStream oos = new ObjectOutputStream(
						new FileOutputStream(f));
				oos.writeObject(m);
				oos.close();
			} catch (Exception e) {
				System.err.println("Error writing file: " + e);
			}
		}
	}
}
