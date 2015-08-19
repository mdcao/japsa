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

package japsa.tools.bio.phylo;

import japsa.bio.phylo.PhylogenyTree;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;



/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.phylo.normalise",
			scriptDesc = "Scale branches of a tree so that the sum of branch lengths is equal to a value")
public class NormaliseTreeCmd  extends CommandLine{	
	public NormaliseTreeCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addStdInputFile();
		addDouble("sum",1.0,"Sum of branches after normalising");
		addString("output", "-", "Name of the file for output, - for stdout");
				
		addStdHelp();		
	} 
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {		
		CommandLine cmdLine = new NormaliseTreeCmd();
		args = cmdLine.stdParseLine(args);
				
		String output = cmdLine.getStringVal("output");		
		double nSum   = cmdLine.getDoubleVal("sum");
		
		BufferedReader bf = SequenceReader.openFile(cmdLine.getStringVal("input"));	
		
		String line = null, str = "";
		while ((line = bf.readLine()) != null) {
			str = str + line.trim();
		}

		bf.close();

		PhylogenyTree tree = PhylogenyTree.parseTree(str);
		double sum = tree.sumHops();
		double num = tree.numHops();

		tree.scale(nSum * num / sum);

		PrintStream ps = System.out;
		if(!"-".equals(output))
			ps = new PrintStream(new FileOutputStream(output));
		
		ps.println(tree.toString() + ";");

		ps.close();
	}
}


/*RST*
-------------------------------------------------------------
*jsa.phylo.normalise*: Normalise branch length of a phylogeny 
-------------------------------------------------------------

*jsa.phylo.normalise* scales the branches of a phylogeny so that their sum
equates to a value.
 
*jsa.phylo.normalise* is included in the `Japsa package <http://mdcao.github.io/japsa/>`_. 
Please see check the installation_ page for instructions.  

.. _installation: ../install.html

<usage>


*RST*/
