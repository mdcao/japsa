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

import japsa.xm.expert.Expert;

public class TrieNode {
	private static final int INITIAL_SIZE = 16;

	private TrieNode[] children;// Should have some children unless it is a leaf
	int count;
	private int[] values = null;
	int level;

	int maxLevel, minLevel;

	/**
	 * Create a normal node
	 * 
	 */
	public TrieNode() {
		children = new TrieNode[Expert.alphabet().size()];
		for (int i = 0; i < Expert.alphabet().size(); i++)
			children[i] = null;
		count = 0;
	}

	/**
	 * Create a node to store values, and depends of whether it is a leaf
	 * 
	 * @param isLeaf
	 */
	public TrieNode(boolean isLeaf) {
		if (!isLeaf) {
			children = new TrieNode[Expert.alphabet().size()];
			for (int i = 0; i < Expert.alphabet().size(); i++)
				children[i] = null;
		}
		count = 0;
		values = new int[INITIAL_SIZE];
	}

	public void addNewChild(int i) {
		children[i] = new TrieNode();
	}

	public void addNewChild(int i, TrieNode trieNode) {
		children[i] = trieNode;
	}

	public TrieNode getChild(int i) {
		return children[i];
	}

	public int getCount() {
		return count;
	}

	/**
	 * Only if the values array has been allocated
	 * 
	 * @param value
	 */
	public void addValue(int value) {
		if (count >= values.length) {// Full
			// reallocate
			int[] array = new int[values.length * 2];
			System.arraycopy(values, 0, array, 0, values.length);
			values = array;
		}
		// assert: enough room for adding
		values[count] = value;
		count++;
	}

	public int[] getValues() {
		return values;
	}

	public int getLevel() {
		return level;
	}

	public void setLevel(int level) {
		this.level = level;
	}

}
