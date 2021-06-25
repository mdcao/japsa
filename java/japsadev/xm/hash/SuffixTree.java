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
import japsa.util.IntIterator;
import japsa.xm.hash.PatternStore;

import java.util.Random;

/**
 * @author Minh Duc Cao
 * 
 */
public class SuffixTree implements PatternStore {
	public static int ALPHABET = 4;
	public int longestS = 0, longestNode = 0;
	int incrementLeaf = 1;

	int minHash = 11;

	byte[] seqs;
	SuffixNode root = new SuffixNode(-1, -1);

	byte[] refSeq;

	SuffixNode currentNode = null, palinNode = null, matchNode = null;

	int currentIndex = 0;

	public SuffixTree(byte[] aSeqs) {
		seqs = aSeqs;
		refSeq = seqs;// For building
	}

	public void setMinHash(int minH) {
		minHash = minH;
	}

	public void setRefSeq(byte[] ref) {
		refSeq = ref;
	}

	public void setIncrementLeaf(int incrementLeaf) {
		this.incrementLeaf = incrementLeaf;
		currentIndex = 0;
	}

	public void clear() {
		root = new SuffixNode(-1, -1);

	}

	public SuffixNode matchNode() {
		SuffixNode node = currentNode;
		if (longestNode != longestS) {
			node = currentNode.children[refSeq[currentIndex - longestNode]];
		}

		return node;
	}

	// Break this into 2: First one to find the Maximum match
	// Second one to put the new posSrc in

	public void next(int baseInd) {
		longestS = 0;
		currentNode = moveNext(root);
		addNode(currentIndex);
		currentIndex++;

		// return currentNode;

	}

	private SuffixNode moveNext(SuffixNode node) {
		if (currentIndex <= longestS)
			return node;

		node.leaves += incrementLeaf;
		// index == currentIndex - longestS;
		// leaves ++;

		// count ++;
		longestNode = longestS;

		if (node.children == null) {
			return node;// Case 1
		}

		if (node.children[refSeq[currentIndex - longestS]] == null) {// No child
																		// at
																		// this
																		// brand
			return node;// case 2

		} else {// the two prefixes are the same, travel down more
			SuffixNode child = node.children[refSeq[currentIndex - longestS]];
			// assert seqHash[child.start] == seqHash[index];
			int childIndex = child.start;

			int endIndex = (child.isLeaf()) ? 0 : child.end;

			while (childIndex >= endIndex
					&& refSeq[currentIndex - longestS] == seqs[childIndex]) {
				longestS++;
				childIndex--;
				if (currentIndex == longestS)
					return node;
			}

			if (childIndex >= endIndex) {// add a node here
				return node;// case 3

			} else {
				return moveNext(child);
			}
		}

	}

	/**
	 * Precond: moveNext has been call, longestS > 0, node = the node nearest to
	 * the longestS
	 * 
	 * @param node
	 * @param posSrc
	 */
	private void addNode(int pos) {
		if (longestNode == longestS) {
			if (currentNode.children == null)
				currentNode.children = new SuffixNode[ALPHABET];
			// assert node.children[seqHash[currentIndex - longestS]] == null
			currentNode.children[seqs[currentIndex - longestS]] = new SuffixNode(
					currentIndex - longestS, pos + 1);
		} else {
			SuffixNode node = new SuffixNode(currentIndex - longestNode,
					currentIndex - longestS + 1);

			node.children = new SuffixNode[ALPHABET];
			SuffixNode.internalCount++;

			SuffixNode child = currentNode.children[seqs[currentIndex
					- longestNode]];
			currentNode.children[seqs[currentIndex - longestNode]] = node;

			child.start -= (longestS - longestNode);
			// Set the previous child
			node.children[seqs[child.start]] = child;

			// Set the new branch=>leaf node
			node.children[seqs[currentIndex - longestS]] = new SuffixNode(
					currentIndex - longestS, pos + 1);
			node.leaves = 1 + child.leaves;

		}

	}

	public void printSummary() {

	}

	public void printTree() {
		root.print(1, true);
	}

	public void nextKey(int baseInd) {
		longestS = 0;
		currentNode = moveNext(root);
	}

	public void putCurrentValue(int val) {
		if (incrementLeaf == 1)
			addNode(val);
		currentIndex++;
	}

	/**
	 * Return the iterator of all previous candidate copy experts
	 * 
	 * @return
	 */
	public IntIterator copyIterator() {
		return new CopyIterator();
	}

	/**
	 * Return the iterator of all previous candidate palin experts
	 * 
	 * @return
	 */
	public IntIterator palinIterator() {
		return new CopyIterator();
	}

	/**
	 * Return the iterator of all previous candidate copy and palin experts They
	 * are in fact the combination of both copy and palin experts
	 * 
	 * @return
	 */
	public IntIterator iterator() {
		return new CopyIterator();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int i = 0;
		try {
			Sequence dna = SequenceReader.getReader(args[0]).nextSequence(null);// (filename)IOTools.read(args[0]);

			System.out.print("Initilise ... ");
			// int size = Integer.parseInt(args[0]);
			// byte[] bytes = new byte[size];
			// Random rnd = new Random();
			// for (i = 0; i < bytes.length; i++){
			// bytes[i] = (byte) rnd.nextInt(ALPHABET);
			// }

			byte[] bytes = dna.toBytes();
			System.out.println("done");

			// byte[] x = {0,1,1,1,0,0,1,1,0,1,0,1,1,0,1,0,1,0,1,0,1};

			// bytes = x;

			SuffixTree tree = new SuffixTree(bytes);
			long now = 0, start = System.currentTimeMillis();

			for (i = 0; i < bytes.length; i++) {
				// System.out.printf("%3d  ",i);
				// System.out.print(bytes[i]);
				// SuffixNode node =

				tree.nextKey(bytes[i]);

				tree.putCurrentValue(i);

				if (i % 1000000 == 0 && i > 1) {
					now = System.currentTimeMillis();
					System.out
							.printf("Milestone %8d , created %8d leaves and %8d internal in %8d ms %f\n",
									i, SuffixNode.count
											- SuffixNode.internalCount,
									SuffixNode.internalCount, (now - start),
									SuffixNode.internalCount * 1.0 / i);

					System.gc();
					Runtime.getRuntime().gc();

					System.gc();
					Runtime.getRuntime().gc();

				}
				// System.out.println(node + "  " + node.count);
				// System.out.println(tree.root.count);
				//
				// tree.printTree();
			}
			if (now == 0)
				now = System.currentTimeMillis();
			System.out.println("\nBuild tree done in " + (now - start) + " ms");

			// for (int i = 0; i < bytes.length; i++){
			// System.out.print(bytes[i]);
			// }

			// System.out.println();
			// tree.printTree();

			System.out
					.printf("Milestone %8d , created %8d leaves and %8d internal in %8d ms %f\n",
							bytes.length, SuffixNode.count
									- SuffixNode.internalCount,
							SuffixNode.internalCount, (now - start),
							SuffixNode.internalCount * 1.0 / bytes.length);

			// tree.printTree();

			/***********************************************
			 * tree.nextKey(); tree.printTree();
			 * 
			 * 
			 * tree.nextKey(); tree.printTree();
			 * 
			 * tree.nextKey(); tree.printTree();
			 * 
			 * tree.nextKey(); tree.printTree();
			 * 
			 * tree.nextKey(); tree.printTree();
			 * 
			 * tree.nextKey(); tree.printTree();
			 * 
			 * /
			 ***********************************************/

		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("End " + i);

		}
		System.out.println(i);
	}

	public Random rnd = new Random(1);

	class CopyIterator implements IntIterator {
		private int total;
		private int[] myArray = null;

		public CopyIterator() {
			if (longestS < minHash || currentNode == root) {
				total = 0;
				return;
			}

			SuffixNode node = currentNode;

			// System.out.println(currentNode + "  " + currentIndex + " " +
			// longestNode
			// + "      " + longestS);
			if (longestNode != longestS) {
				node = currentNode.children[refSeq[currentIndex - longestNode]];
			}

			// System.out.println(node );
			myArray = new int[node.leaves];
			total = 0;
			populate(node);
			// assert total = currentNode.leaves;
		}

		public boolean hasNext() {
			return total > 0;
		}

		public int sizeAvailable() {
			return total;
		}

		public int next() {
			int random_index = rnd.nextInt(total);
			// Shuffle the array

			if (random_index + 1 < total) {// Move the expert selected to end of
											// the list
				// so wont be selected again
				int tmp = myArray[random_index];
				myArray[random_index] = myArray[total - 1];
				myArray[total - 1] = tmp;
			}// else random_index +1 == total

			total--;

			return myArray[total];
		}

		private void populate(SuffixNode node) {
			if (node.isLeaf()) {
				myArray[total] = node.end - 1;

				total++;
			}

			if (node.children == null)
				return;
			for (int i = 0; i < ALPHABET; i++) {
				if (node.children[i] != null) {
					populate(node.children[i]);
				}

			}

		}
	}
}

// Leaf is the one that has end > start. End is actually the position

class SuffixNode {
	// Can compress start/end/leaves into 1
	// 30 bit ~1B
	// 100Tr ->27 bits
	// 3 numbers=> 54 bits=>two ints

	public static int count = 0;
	public static int internalCount = 0;

	int start, end;
	int leaves = 0;// Number of leaves can be reached from here

	SuffixNode[] children;// = new SuffixNode[4];

	public SuffixNode(int s, int e) {
		start = s;
		end = e;
		count++;
		leaves = (isLeaf()) ? 1 : 0;
	}

	public boolean isLeaf() {
		return (end > start);
	}

	public void print(int offSet, boolean newLine) {
		if (newLine)
			for (int i = 0; i < offSet; i++)
				System.out.print(' ');

		System.out.print("----" + this);
		boolean first = true;
		if (children != null) {
			for (int i = 0; i < children.length; i++) {

				if (children[i] != null) {
					children[i].print(offSet + 12, !first);
					first = false;
					// System.out.println();
				}
			}
		}

		System.out.println();
	}

	public String toString() {
		// return "(" + start + "," + count +"," + end +")";
		return "(" + start + "," + leaves + "," + end + ")";
	}

}
