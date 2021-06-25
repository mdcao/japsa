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

package japsadev.tools.work;

import japsa.bio.phylo.PhylogenyTree;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Random;



/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.work.treehet",
scriptDesc = "Make the branch of the tree with random model")
public class TreeHetLinageCmd  extends CommandLine{	
	public TreeHetLinageCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addStdInputFile();		
		addString("output", "-", "Name of the file for output, - for stdout");
		addString("model", "m", "Model prefix");

		addStdHelp();		
	} 
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {		
		CommandLine cmdLine = new TreeHetLinageCmd();
		args = cmdLine.stdParseLine(args);				
		String output = cmdLine.getStringVal("output");
		prefix = " #" + cmdLine.getStringVal("model");
		BufferedReader bf = SequenceReader.openFile(cmdLine.getStringVal("input"));	

		String line = null, str = "";
		while ((line = bf.readLine()) != null) {
			str = str + line.trim();
		}
		bf.close();
		PhylogenyTree tree = PhylogenyTree.parseTree(str);

		PrintStream ps = System.out;
		if(!"-".equals(output))
			ps = new PrintStream(new FileOutputStream(output));

		ps.println(ranTree(tree) + prefix + rand.nextInt(4)+";");
		ps.close();
	}
	
	static Random rand = new Random();
	static String prefix = "";

	public static String ranTree(PhylogenyTree tree){
		if (tree.isLeaf()) {
			return tree.getName();
		}
		int b1 = rand.nextInt(4);
		int b2 = rand.nextInt(4);
		
		String s = "(" + ranTree(tree.getChild(0)) + prefix + b1 + "," + ranTree(tree.getChild(1)) + prefix+b2 + ")";		
		return s;
	}
}
