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

package japsa.bio.phylo;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.Iterator;
import java.util.Stack;
import java.util.StringTokenizer;
import java.util.Vector;




/**
 * @author Minh Duc Cao www.caominhduc.org
 * 
 */

public class PhylogenyTree {
	public static double DEFAULT_LENGTH = .4999;

	private PhylogenyTree parent = null;
	private PhylogenyTree[] children = null;
	private double distance;// Distance to its parent
	private Sequence seq;

	// The index of this node in its parent list
	private int treeIndex;
	String name;

	private double height;// The height of the tree (the node)

	public PhylogenyTree(PhylogenyTree parent) {
		this.parent = parent;
		distance = DEFAULT_LENGTH;
		height = 0;
	}
	
	// Create an internal tree
	public PhylogenyTree(PhylogenyTree parent, PhylogenyTree left,
			PhylogenyTree right) {
		this(parent);

		children = new PhylogenyTree[2];
		setChild(left, 0);
		setChild(right, 1);
	}
	
	/**
	 * Create a new tree who where the the root is moved to the parent of the
	 * node A new node to be created (the root) and returned
	 * 
	 * @param node
	 *            : A leaf node whom the new root will be the parent of
	 * @return : The new root node
	 */

	public PhylogenyTree moveRootTo(PhylogenyTree node) {
		PhylogenyTree oldParen = node.parent;
		if (oldParen == null)
			return this;// Need not do anything

		PhylogenyTree newNode = oldParen.addChild(new PhylogenyTree(null), node
				.getIndex(), node.distance / 2);

		newNode.getChild(1).setName("Tmp");
		newNode.modifiedRoot(1);
		newNode.parent = null;

		return newNode;
	}
	// how far I am from the farther desentdent?
	public int distanceToBottom() {
		if (isLeaf())
			return 0;
		int ld = children[0].distanceToBottom();
		int rd = children[1].distanceToBottom();

		if (ld > rd)
			return ld + 1;
		else
			return rd + 1;
	}

	// How far I am from the most ancestor
	public int level() {
		if (parent == null)
			return 0;
		else
			return parent.level() + 1;
	}

	/**
	 * Find a leaf node with name
	 * 
	 * @param aName
	 * @return
	 */
	public PhylogenyTree findLeaf(String aName) {
		if (isLeaf())
			if (this.name.equals(aName))
				return this;
			else
				return null;

		PhylogenyTree target = children[0].findLeaf(aName);
		if (target != null)
			return target;
		else
			return children[1].findLeaf(aName);
	}

	/**
	 * Find the left most leaf of the tree
	 * 
	 * @return
	 */
	public PhylogenyTree findLeftMost() {
		if (isLeaf())
			return this;
		return children[0].findLeftMost();
	}

	// Change parent to the (index)th child
	private void modifiedRoot(int index) {

		if (parent == null)
			return;

		if (parent.parent == null) {
			PhylogenyTree sibling = parent.getChild(1 - this.treeIndex);
			sibling.distance += this.distance;
			setChild(sibling, index);

		} else {
			parent.modifiedRoot(this.treeIndex);
			setChild(parent, index);
		}
	}

	public PhylogenyTree clone() {
		if (isLeaf()) {
			PhylogenyTree aLeaf = new PhylogenyTree(null);
			aLeaf.height = this.height;
			aLeaf.distance = this.distance;
			aLeaf.name = this.name;
			aLeaf.seq = this.seq;
			return aLeaf;
		}
		PhylogenyTree aNode = new PhylogenyTree(null);
		aNode.children = new PhylogenyTree[2];
		aNode.setChild(children[0].clone(), 0);
		aNode.setChild(children[1].clone(), 1);
		aNode.distance = this.distance;
		aNode.height = this.height;
		// if (compareTree(aNode))
		return aNode;

		// System.err.println("Error " + this);
		// return null;
	}

	public PhylogenyTree(PhylogenyTree parent,Sequence aSeq, String aName) {
		this(parent);

		setSeq(aSeq);
		setName(aName);
	}

	public boolean compareUnRootedTree(PhylogenyTree aTree) {

		PhylogenyTree newTree = aTree.findLeftMost();
		PhylogenyTree myNewTree = this.findLeaf(newTree.getName());
		if (myNewTree == null)
			return false;

		newTree = aTree.moveRootTo(newTree);
		myNewTree = this.moveRootTo(myNewTree);

		System.out.println(newTree);
		System.out.println(" vs ");
		System.out.println(myNewTree);

		return myNewTree.compareRootedTree(newTree);
	}

	public boolean compareRootedTree(PhylogenyTree aTree) {
		if (isLeaf()) {
			if (aTree.isLeaf())
				return aTree.name.equals(name);
			else
				return false;
		}
		// Is not a leaf
		if (aTree.isLeaf())
			return false;

		// Both trees are internal
		if (getChild(0).compareRootedTree(aTree.getChild(0))) {
			if (getChild(1).compareRootedTree(aTree.getChild(1))) {
				return true;
			} else {
				System.out.println("Left == left but right != right");
				System.out.println(getChild(0) + "  ==  " + aTree.getChild(0));
				System.out.println(getChild(1) + "  !=  " + aTree.getChild(1));
				System.out
				.println("********************************************");
				return false;
			}
		}
		if (getChild(0).compareRootedTree(aTree.getChild(1))) {
			if (getChild(1).compareRootedTree(aTree.getChild(0))) {
				return true;
			} else {

				System.out.println("left == right but right != left");
				System.out.println(getChild(0) + "  ==  " + aTree.getChild(1));
				System.out.println(getChild(1) + "  !=  " + aTree.getChild(0));
				System.out.println("=========================================");

				return false;
			}

		}

		System.out.println(getChild(0) + "  !=  " + aTree.getChild(0));
		System.out.println(getChild(0) + "  !=  " + aTree.getChild(1));
		System.out
		.println("--------------------------------------------------");

		return false;
	}

	public static PhylogenyTree readFromFile(String fileName) throws Exception {
		BufferedReader bf = new BufferedReader(new FileReader(fileName));

		String line = null;
		String treeStr = "";
		while ((line = bf.readLine()) != null) {
			treeStr += line.trim();
		}
		bf.close();

		return parseTree(treeStr.trim());
	}

	static int LATEXDIS = 5;// Distance in latex between branches (very much the
	// height of the image)
	// the height = LATEXDIS * # of leaves
	static int CRT_LATEXY = 60;// starting point (highest y coordinator)
	static int EXTRA_WIDTH = 0; //

	static int INDEX = 0;// the first index of the node

	static String formatName(String in) {
		return in.replace('_', ' ');
	}

	int latexDrawTree(PrintStream out, boolean len) {
		NumberFormat nf = NumberFormat.getInstance();

		nf.setMaximumFractionDigits(2);

		if (isLeaf()) {
			// Draw the node
			int x = EXTRA_WIDTH;
			INDEX++;
			out.println("\\node[exnode,label=right:\\emph{"
					+ formatName(this.name) + "}]  at (" + x + "," + CRT_LATEXY
					+ ") (node" + INDEX + ") {};");
			// draw the name of the node
			int y = CRT_LATEXY;
			CRT_LATEXY -= LATEXDIS;
			return y;

		} else {
			int y0 = children[0].latexDrawTree(out, len);
			int ind0 = INDEX;
			int y1 = children[1].latexDrawTree(out, len);
			int ind1 = INDEX;

			INDEX++;
			int y = (y0 + y1) / 2;
			int x = level() * 5;
			out.println("\\node[innode]  at (" + x + "," + y + ") (node"
					+ INDEX + ") {};");
			// out.println("  ;");
			if (len) {
				out.println("\\node at (" + (x + 2.0) + "," + (y0 + 1.2)
						+ ") {" + nf.format(children[0].distance) + "};");
				out.println("\\node at (" + (x + 2.0) + "," + (y1 + 1.2)
						+ ") {" + nf.format(children[1].distance) + "};");
			}
			out.println("\\path[edge] (node" + INDEX + ")  -- (" + x + "," + y0
					+ ") -- (node" + ind0 + "); ");
			out.println("\\path[edge] (node" + INDEX + ")  -- (" + x + "," + y1
					+ ") -- (node" + ind1 + "); ");

			return y;
		}
	}

	public void drawTree(String fileName, boolean len) throws Exception {
		PrintStream out = new PrintStream(new FileOutputStream(fileName));
		out
		.println("\\tikzstyle{exnode}=[circle,fill=red!75,minimum size=4pt,inner sep=1pt];");
		out
		.println("\\tikzstyle{innode}=[circle,fill=green!75!blue,minimum size=3pt,inner sep=1pt];");
		out.println("\\tikzstyle{edge} = [draw,line width=1pt,-,blue!50];\n");

		EXTRA_WIDTH += distanceToBottom() * 5;

		latexDrawTree(out, len);

		out.close();
	}

	/**
	 * treeStr doesn't contain white spaces
	 * 
	 * @param treeStr
	 *            :
	 * @return
	 * @throws Exception
	 */
	public static PhylogenyTree parseTree(String treeStr){
		double defDistanceRead = 0.0;
		Stack<PhylogenyTree> t = new Stack<PhylogenyTree>();
		PhylogenyTree tree = null;

		StringTokenizer stk = new StringTokenizer(treeStr, "\n ,\t():;[]", true);
		String currentToken = "";

		while (stk.hasMoreElements()) {
			currentToken = stk.nextToken();

			if (";".equals(currentToken))
				break;
			if ("[".equals(currentToken))
				break;
			if ("]".equals(currentToken))
				break;

			if (" ".equals(currentToken) || "\t".equals(currentToken) || "\n".equals(currentToken)
					|| ",".equals(currentToken))
				continue;

			if ("(".equals(currentToken)) {// Start a new tree
				// PhylogenyContinuousTree tree = new
				// PhylogenyContinuousTree(null);
				// t.push(tree);

			} else if (")".equals(currentToken)) {// End a new tree
				PhylogenyTree right = t.pop();
				PhylogenyTree left = t.pop();
				tree = new PhylogenyTree(null, left, right);
				tree.distance = defDistanceRead;
				t.push(tree);				
			}else if (":".equals(currentToken)) {// End a new tree
				double dis = Double.parseDouble(stk.nextToken());
				t.peek().setDistance(dis);
			} else {// no bracket, supposed to be in format name
				tree = new PhylogenyTree(null);
				tree.distance = defDistanceRead;

				tree.name = currentToken;
				t.push(tree);// A leaf

			}
		}

		tree = t.pop();
		while (!t.isEmpty()) {
			tree = new PhylogenyTree(null, t.pop(), tree);
		}
		return tree;
	}

	/**
	 * 
	 *Add a newChild into the left (0) or right PreCond: this is an internal
	 * tree Return the new created child
	 */
	public PhylogenyTree addChild(PhylogenyTree aTree, int index, double dis) {

		children[index].distance -= dis;
		PhylogenyTree newChild = new PhylogenyTree(this, this.getChild(index),
				aTree);

		newChild.setDistance(dis);
		aTree.distance = newChild.height = this.height - dis;

		setChild(newChild, index);
		return newChild;
	}

	public void setHeight(double h) {
		height = h;
	}

	public PhylogenyTree removeGrandChild(int childIndex, int grandChildIndex) {

		PhylogenyTree child = children[childIndex];
		PhylogenyTree theGrandChild = child.children[grandChildIndex];

		this.setChild(child.children[1 - grandChildIndex], childIndex);

		// Restore the distance
		this.getChild(childIndex).distance = this.height
				- this.getChild(childIndex).height;

		return theGrandChild;
	}

	/**
	 * 
	 * @param subtree
	 * @param index
	 */
	public void setChild(PhylogenyTree subtree, int index) {
		// The subtree should have a height already
		children[index] = subtree;
		subtree.treeIndex = index;
		subtree.parent = this;

	}

	public PhylogenyTree getChild(int index) {
		if (children == null)
			return null;// Leaf
		else
			return children[index];
	}

	public boolean isLeaf() {
		return children == null;
	}

	public boolean isRoot() {
		return parent == null;
	}

	public static PhylogenyTree createLeaf(String aName, Sequence aSeq) {
		PhylogenyTree leaf = new PhylogenyTree(null);

		leaf.name = aName;
		leaf.seq = aSeq;

		return leaf;
	}

	/**
	 * Find the nearest common (conceptially) ancestor of the two trees (or two
	 * nodes) The real common ancestor is somewhere in the path Make it
	 * protected to avoid abusing it
	 * 
	 * @param aTree
	 * @return
	 */
	protected PhylogenyTree commonAncestor(PhylogenyTree aTree) {
		StringBuffer routA = new StringBuffer(), routB = new StringBuffer();
		PhylogenyTree treePtr = this;

		while (treePtr.parent != null) {// Not a root
			routA.append(treePtr.treeIndex);
		}

		treePtr = aTree;
		while (treePtr.parent != null) {// Not a root
			routB.append(treePtr.treeIndex);
		}

		// Assert aTree = root
		// Now go down the tree
		int indA = routA.length(), indB = routB.length();

		while (indA >= 0 && indB >= 0
				&& routA.charAt(indA) == routB.charAt(indB)) {

			int ind = Integer.parseInt(routA.charAt(indA) + "");
			treePtr = treePtr.getChild(ind);
			// dis -= iNode.distance * 2;

			indA--;
			indB--;
		}

		return treePtr;
	}

	public double distanceTo(PhylogenyTree aTree) {
		double dis = 0;

		StringBuffer routA = new StringBuffer(), routB = new StringBuffer();
		PhylogenyTree treePtr = this;

		while (treePtr.parent != null) {// Not a root
			routA.append(treePtr.treeIndex);
			dis += treePtr.distance;
			treePtr = treePtr.parent;
		}

		treePtr = aTree;
		while (treePtr.parent != null) {// Not a root
			routB.append(treePtr.treeIndex);

			dis += treePtr.distance;
			treePtr = treePtr.parent;
		}

		// Assert aTree = root
		// Now go down the tree
		int indA = routA.length() - 1, indB = routB.length() - 1;

		while (indA >= 0 && indB >= 0
				&& routA.charAt(indA) == routB.charAt(indB)) {

			int ind = Integer.parseInt(routA.charAt(indA) + "");
			treePtr = treePtr.getChild(ind);
			dis -= treePtr.distance * 2;

			indA--;
			indB--;
		}

		//
		// if (treePtr.parent == null) dis -= 1;

		return dis;
	}

	public Iterator<PhylogenyTree> getInternalIterator() {
		Vector<PhylogenyTree> v = new Vector<PhylogenyTree>();
		this.addInternalTrees(v);
		return v.iterator();
	}

	public Iterator<PhylogenyTree> getLeafIterator() {
		Vector<PhylogenyTree> v = new Vector<PhylogenyTree>();
		this.addLeafTrees(v);
		return v.iterator();
	}

	private void addInternalTrees(Vector<PhylogenyTree> v) {
		// Preorder
		if (!isLeaf()) {
			v.add(this);
			children[0].addInternalTrees(v);
			children[1].addInternalTrees(v);
		}
	}

	private void addLeafTrees(Vector<PhylogenyTree> v) {
		// Preorder
		if (!isLeaf()) {
			children[0].addLeafTrees(v);
			children[1].addLeafTrees(v);
		} else
			v.add(this);
	}

	public String toString() {
		if (isLeaf()) {
			return "" + name;
		}
		String tree = "(" + children[0] + ":" + children[0].getDistance() + ","
				+ children[1] + ":" + children[1].getDistance() + ")";
		// if (parent == null) return tree + ";";//If the root
		// else
		return tree;

	}

	/**
	 * Scale ALL distances by a factor f (this trees and all subtrees)
	 * 
	 * @param f
	 */
	public void scale(double f) {
		this.distance *= f;
		if (children != null) {
			if (children[0] != null)
				children[0].scale(f);
			if (children[1] != null)
				children[1].scale(f);
		}
	}

	public double sumHops() {
		if (this.isLeaf())
			return 0;
		else
			return this.getChild(0).sumHops() + this.getChild(0).distance
					+ this.getChild(1).sumHops() + this.getChild(1).distance;
	}

	public int numHops() {
		if (this.isLeaf())
			return 0;
		else
			return this.getChild(0).numHops() + 1
					+ this.getChild(1).numHops() + 1;
	}

	/************************* Get and Set ********************************/
	public PhylogenyTree getParent() {
		return parent;
	}

	public void setParent(PhylogenyTree parent) {
		this.parent = parent;
	}

	public PhylogenyTree[] getChildren() {
		return children;
	}

	public void setChildren(PhylogenyTree[] children) {
		this.children = children;
	}

	public double getDistance() {
		return distance;
	}

	public void setDistance(double distance) {
		this.distance = distance;
	}

	public Sequence getSeq() {
		return seq;
	}

	public void setSeq(Sequence seq) {
		this.seq = seq;
	}

	public int getIndex() {
		return treeIndex;
	}

	public void setIndex(int index) {
		this.treeIndex = index;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public static void main(String[] args) throws Exception {

	}

	public static void writeNexus(Sequence[] seqs, PrintStream ps)
			throws IOException {

		int length = seqs[0].length();
		int charPerLine = 60;//	

		//PrintStream ps = new PrintStream(new FileOutputStream(fileName));
		ps.println("#NEXUS\nBEGIN DATA;");
		ps.println("        Dimensions NTax=" + seqs.length + " NChar="
				+ length + ";");
		ps
		.println("        Format DataType=DNA Interleave=yes Gap=- Missing=?;");
		ps.println("        Matrix");
		int count = 0;

		while (true) {
			for (int i = 0; i < seqs.length; i++) {
				ps.printf("%s ", (seqs[i].getName() + "             ")
						.substring(0, 10));

				for (int x = count; x < count + charPerLine && x < length; x++) {
					if (x % 10 == 0 && x > count)
						ps.print(' ');
					ps.print(seqs[i].charAt(x));
				}
				ps.println();
			}

			ps.println();

			count += charPerLine;
			if (count >= length)
				break;
		}
		ps.println("         ;\nEND;");
		ps.flush();
		//ps.close();

	}

	public static void writePhylip(Sequence[] seqs, PrintStream ps)
			throws IOException {

		//PrintStream ps = new PrintStream(new FileOutputStream(fileName));
		int length = seqs[0].length();
		int charPerLine = length + 12;

		ps.println(seqs.length + "   " + length);
		int count = 0;

		while (true) {
			for (int i = 0; i < seqs.length; i++) {
				if (count == 0) {
					ps.printf("%s ", (seqs[i].getName() + "             ")
							.substring(0, 10));
				}

				for (int x = count; x < count + charPerLine && x < length; x++) {
					if (x % 10 == 0 && x > count)
						ps.print(' ');
					ps.print(seqs[i].charAt(x));
				}
				ps.println();
			}

			ps.println();

			count += charPerLine;
			if (count >= length)
				break;
		}
		ps.flush();
		//ps.close();
	}

	public static Sequence[] readPhylip(String fileName) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(fileName));	

		String line = br.readLine();
		StringTokenizer stoken = new StringTokenizer(line);
		int noSeq = Integer.parseInt(stoken.nextToken());
		int length = Integer.parseInt(stoken.nextToken());

		Sequence[] seqs = new Sequence[noSeq];
		int[] count = new int[noSeq];		

		int seqIndex = 0;		
		while ((line = br.readLine()) != null){
			line = line.trim();
			if (line.length() == 0)
				continue;

			int lineIndex = 0;
			if (seqs[seqIndex] == null){
				seqs[seqIndex] = new Sequence(Alphabet.DNA4(), length, line.substring(0,10).trim());
				lineIndex = 10;
			}	

			for (;lineIndex < line.length(); lineIndex++){
				char c = line.charAt(lineIndex);
				int b = seqs[seqIndex].alphabet().char2int(c);
				if (b >=0)
					seqs[seqIndex].setBase(count[seqIndex] ++ , (byte) b);				
			}
			//seqIndex ++;
			if (++seqIndex >= noSeq)
				seqIndex = 0;			
		}
		br.close();
		return seqs;
	}
}
