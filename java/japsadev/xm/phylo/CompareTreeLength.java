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

package japsadev.xm.phylo;

import japsa.bio.phylo.PhylogenyTree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Iterator;
import java.util.Vector;

public class CompareTreeLength {

	public static double compare(PhylogenyTree treeA, PhylogenyTree treeB) {

		Iterator<PhylogenyTree> iterA = treeA.getLeafIterator();

		Vector<PhylogenyTree> a = new Vector<PhylogenyTree>();
		while (iterA.hasNext()) {
			a.add(iterA.next());
		}

		PhylogenyTree[] aA = new PhylogenyTree[a.size()];
		a.toArray(aA);

		// sort array A
		for (int i = 1; i < aA.length; i++) {
			for (int j = 0; j < i; j++) {
				if (aA[i].toString().compareTo(aA[j].toString()) < 0) {
					// swap
					PhylogenyTree tmp = aA[j];
					aA[j] = aA[i];
					aA[i] = tmp;
				}
			}
		}

		double[][] disA = new double[aA.length][aA.length];
		// double sumA = 0;
		for (int i = 1; i < aA.length; i++) {
			for (int j = 0; j < i; j++) {
				disA[i][j] = disA[j][i] = aA[i].distanceTo(aA[j]);
				// sumA += disA[i][j];
			}
		}

		Iterator<PhylogenyTree> iterB = treeB.getLeafIterator();

		Vector<PhylogenyTree> b = new Vector<PhylogenyTree>();
		while (iterB.hasNext()) {
			b.add(iterB.next());
		}

		PhylogenyTree[] aB = new PhylogenyTree[b.size()];
		b.toArray(aB);

		// sort array A
		for (int i = 1; i < aB.length; i++) {
			for (int j = 0; j < i; j++) {
				if (aB[i].toString().compareTo(aB[j].toString()) < 0) {
					// swap
					PhylogenyTree tmp = aB[j];
					aB[j] = aB[i];
					aB[i] = tmp;
				}
			}
		}

		double[][] disB = new double[aB.length][aB.length];

		// double sumB = 0;
		for (int i = 1; i < aB.length; i++) {
			for (int j = 0; j < i; j++) {
				disB[i][j] = disB[j][i] = aB[i].distanceTo(aB[j]);
			}
		}

		// for (int i = 0; i < aB.length; i++) {
		// System.out.println(aB[i]+ "   " + aA[i]);
		// }
		double sum = 0;

		for (int i = 1; i < aB.length; i++) {
			for (int j = 0; j < i; j++) {
				double v = (disB[i][j] - disA[i][j]);
				sum += v * v;
			}
		}

		return sum;
	}

	public static void main(String[] args) throws Exception {
		BufferedReader bf = new BufferedReader(new FileReader(args[0]));

		String line = "";
		String dndA = "";// ((n07:93.344,( n02:78.0547, (n04:71.9283, (
							// n08:55.1436,(n09:11.2462,
							// n01:11.2462):43.8974):16.7848):6.12634):15.2893
							// ):6.63599, ( n06:84.565, ( n03:73.6297, (
							// n05:22.5472, n00:22.5472 ):51.0826 ):10.9352
							// ):15.415);";
		String dndB = "";// ((((n00:0.18065,n05:0.17464):0.25103,n03:0.34336):0.16730,n06:0.31645):0.18090,(n04:0.28654,(((n01:0.09601,n09:0.10990):0.26280,n08:0.27821):0.17186,n02:0.31050):0.15744):0.16279,n07:0.33363)[0.5000];";

		while ((line = bf.readLine()) != null) {
			dndA = dndA + line.trim();
		}

		bf.close();

		bf = new BufferedReader(new FileReader(args[1]));
		while ((line = bf.readLine()) != null) {
			dndB = dndB + line.trim();
		}

		bf.close();

		PhylogenyTree treeA = PhylogenyTree.parseTree(dndA);
		PhylogenyTree treeB = PhylogenyTree.parseTree(dndB);

		String tmp = "";
		if (args.length > 2) {
			tmp = "  " + args[2];
		}

		System.out.println("## " + compare(treeA, treeB) + tmp);
	}

}
