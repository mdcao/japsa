/*****************************************************************************
 * Copyright (c) 2010 Minh Duc Cao, Monash University.  All rights reserved. *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the name of Monash University nor the names of its contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

package japsadev.xm.hash;

import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

/**
 * @author Minh Duc Cao
 * 
 */
// A Collection of suffix array construction algorithm
public class SuffixArrayConstruction {

	static {
		// Load library
		System.loadLibrary("sarray");
	}

	/**
	 * s array M&M algorithm.
	 */
	protected static native void jsarray(String s, int a[], int n);

	/**
	 * Hybrid algorithm s: array of chars a:of size length(s) + 1
	 */
	protected static native void jbsarray(byte s[], int a[], int n);

	protected static native void jlcp(int a[], byte[] s, int b[], int n);

	/**
	 * Deep - Shallow
	 */
	protected static native void jds(byte s[], int a[], int n);

	/**
	 * Divsufsort s: array of chars a:of size length(s) + 1
	 */
	protected static native void jdivsufsort(byte s[], int a[], int n);

	public static int[] construct(byte[] seq) {
		int[] a = new int[seq.length + 1];
		jdivsufsort(seq, a, seq.length);
		return a;
	}

	public static int[] ssarray(byte s[]) {
		int[] a = new int[s.length];
		for (int i = 0; i < a.length; i++)
			a[i] = s[i];

		int ORIG = 1 << 31, // ~(~0>>1), /* sign bit */
		BUCK = 1 << 31;// ~(~0>>1);

		int h, i, j, l;
		int k = 0; /* initialized for lint */
		int[] p;

		int t;
		j = 4;// dna size

		p = new int[a.length];

		for (i = 0; i < a.length; i++)
			/* (0) initialize */
			p[i] = i | ORIG;// Turn the first bit to 1?? why

		for (h = 0;; h = h == 0 ? 1 : 2 * h) {

			for (i = 0; i < a.length; i++) { /* (1) link */
				for (j = p[i]; j >= 0; j = a[j])
					;
				t = (j & ~ORIG) - h;
				j = t < 0 ? t + a.length : t;
				l = a[j];
				a[j] = p[l];
				p[l] = j;
			}

			if (h == 0) { /* find k */
				for (k = 0; k < a.length; k++)
					if (p[k] < 0)
						break;
			}

			for (i = a.length; --k >= 0;) { /* (2) order */
				j = p[k];
				do
					p[--i] = j;
				// while(((j=al[j]) & ORIG) == 0);
				while ((j = a[j]) >= 0);
				p[i] |= BUCK;
			}

			for (i = 0; i < a.length; i++) { /* (3) reconstruct */
				// if((p[i] & BUCK) != 0)
				if (p[i] < 0)
					k++;
				a[p[i] & ~BUCK] = k;
			}

			for (i = 0, j = -1; i < a.length; i++, j = l) { /* (4) refine */
				t = (p[i] & ~BUCK) + h;

				l = a[t >= a.length ? t - a.length : t];

				if (l != j)
					p[i] |= BUCK;

			}

			for (i = 0, k = -1; i < a.length; i++) { /* (5) recode */
				if (p[i] < 0)
					k++;

				a[p[i] & ~BUCK] = k;
				p[i] |= ORIG; /* (0b) reinitialize */
			}
			if (++k >= a.length)
				break;
		}

		for (i = 0; i < a.length; i++) {
			a[i] = p[i] & ~ORIG;

		}
		return a;
	}

	/**
	 * 
	 */
	public SuffixArrayConstruction() {
		// TODO Auto-generated constructor stub
	}

	public static int[] suffixArray(byte[] st) {
		int[] a = new int[st.length + 1];
		jds(st, a, st.length);
		return a;
	}

	public static void print(byte[] s, int[] array, int[] lcp) {
		for (int i = 0; i < array.length; i++) {
			System.out.printf("%5d %5d ", lcp[i], array[i]);
			for (int j = array[i]; j < s.length; j++) {
				System.out.printf("%2d ", s[j]);
			}
			System.out.println();
		}
	}

	/**
	 * Compute the lcp of the suffix array a, on sequence s0
	 * 
	 * @param a
	 * @param s0
	 * @return
	 */

	static int[] lcp(int[] a, byte[] s0) {

		int[] lcp = new int[a.length];
		int i, h;
		int[] inv = new int[a.length];

		for (i = 0; i < a.length; i++)
			inv[a[i]] = i;

		h = 0; /* visit in string order */
		for (i = 0; i < a.length - 1; i++) { /* omit last, least suff */
			int x = inv[i]; /* i,j,x,h as in intro */
			int j = a[x - 1];
			int p1 = i + h;
			int p0 = j + h;

			while (p1 < s0.length && p0 < s0.length && s0[p1++] == s0[p0++])
				h++;
			lcp[x] = h;
			if (h > 0)
				h--;
		}
		lcp[0] = 0; /* least suffix has no predecessor */

		return lcp;
	}

	static int sufcheck(byte[] T, int[] SA) {
		int ALPHABET_SIZE = 4;
		int[] C = new int[ALPHABET_SIZE];
		int i = 0, p, t = 0;
		int c;
		int err = 0;
		int n = T.length;

		/* ranges. */
		if (err == 0) {
			for (i = 0; i <= n; ++i) {
				if ((SA[i] < 0) || (n < SA[i])) {
					err = -2;
					break;
				}
			}
		}

		/* first characters. */
		if (err == 0) {
			for (i = 1; i < n; ++i) {
				if (T[SA[i]] > T[SA[i + 1]]) {
					err = -3;
					break;
				}
			}
		}

		/* suffixes. */
		if (err == 0) {
			for (i = 0; i < ALPHABET_SIZE; ++i) {
				C[i] = 0;
			}
			for (i = 0; i < n; ++i) {
				++C[T[i]];
			}
			for (i = 0, p = 1; i < ALPHABET_SIZE; ++i) {
				t = C[i];
				C[i] = p;
				p += t;
			}

			for (i = 0; i <= n; ++i) {
				p = SA[i];
				if (0 < p) {
					c = T[--p];
					t = C[c];
				} else {
					p = n;
					c = -1;
					t = 0;
				}
				if (p != SA[t]) {
					err = -4;
					break;
				}
				if (0 <= c) {
					++C[c];
					if ((n < C[c]) || (T[SA[C[c]]] != c)) {
						C[c] = -1;
					}
				}
			}
		}

		return err;
	}

	/**
	 * @param args
	 */

	public static void main(String[] args) throws Exception {

		Sequence dna = SequenceReader.getReader(args[0]).nextSequence(null);// (filename)IOTools.read(args[0]);

		System.out.print("Initilise ... ");

		// byte[] s = (new BioCompDNA(s1,"aa")).toBytes();

		byte[] s = dna.toBytes();
		// for (int i = 0; i< s.length; i++)
		// System.out.print(s[i]);

		System.out.println();
		System.out.println("Done  ");
		long timeEnd, timeStart = System.currentTimeMillis();
		int[] a;// = new int[s.length + 1];

		/****************************************************************************
		 * { a = new int[s.length + 1]; System.out.println(
		 * "Testing for jbsarray ************************************** ");
		 * Runtime.getRuntime().runFinalization (); Runtime.getRuntime().gc ();
		 * Thread.currentThread ().yield ();
		 * 
		 * timeStart = System.currentTimeMillis(); jbsarray(s,a,s.length);
		 * timeEnd = System.currentTimeMillis();
		 * 
		 * System.out.println("Running time = " +(timeEnd - timeStart));
		 * System.out.println("Checking = " + (sufcheck(s,a) == 0));
		 * 
		 * } /
		 ****************************************************************************/
		{
			a = new int[s.length + 1];
			System.out
					.println("Testing for jdivsufsort ************************************** ");
			Runtime.getRuntime().runFinalization();
			Runtime.getRuntime().gc();

			// Thread.currentThread ().yield ();

			timeStart = System.currentTimeMillis();
			jdivsufsort(s, a, s.length);
			timeEnd = System.currentTimeMillis();

			System.out.println("Running time = " + (timeEnd - timeStart));
			System.out.println("Checking = " + (sufcheck(s, a) == 0));

		}

		/****************************************************************************/
		/*************************************************************
		 * char [] s1 = dna.getCharSequence(); byte[] s2 = new byte[s1.length];
		 * for (int i = 0; i < s1.length; i++) s2[i] = (byte)s1[i];
		 * /************************************************************* { a =
		 * new int[s.length + 1]; System.out.println(
		 * "Testing for jds ************************************** ");
		 * Runtime.getRuntime().runFinalization (); Runtime.getRuntime().gc ();
		 * Thread.currentThread ().yield ();
		 * 
		 * timeStart = System.currentTimeMillis(); jds(s,a,s.length); timeEnd =
		 * System.currentTimeMillis();
		 * 
		 * System.out.println("Running time = " +(timeEnd - timeStart));
		 * System.out.println("Checking = " + (sufcheck(s,a) == 0));
		 * //print(s,a,h);
		 * 
		 * } /
		 *************************************************************/

		int[] h = lcp(a, s);
		// jlcp(a,s,h,a.length);
		print(s, a, h);

	}

}
