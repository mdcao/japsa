/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
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
 * 3. Neither the names of the institutions nor the names of the contributors*
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

/*                           Revision History                                
 * 30/11/2015 - Minh Duc Cao: Revised
 ****************************************************************************/

package japsadev.tools;


import java.io.BufferedReader;
import java.util.HashSet;
import java.util.Iterator;

import japsa.bio.phylo.PhylogenyTree;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.pcgene",
	scriptDesc = "Pan- and core- genes analyses of bacterial population genomics"
	)
public class PanCoreGeneCmd extends CommandLine{
	public PanCoreGeneCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		this.addString("config", "sample.cfg", "Configure file");
		this.addString("stage", null,  "Stage, one of: status, assemble, annotate, clean, analyseVar",true);
		this.addBoolean("force", false,  "Force to overwrite existing results");
		this.addInt("thread", PanCoreGenes.threads , "Number of threads");		

		//addBoolean("reverse",false,"Reverse sort order");	

		addStdHelp();		
	} 

	public static void main(String[] args) throws Exception {
		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new PanCoreGeneCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String config = cmdLine.getStringVal("config");
		String stage  = cmdLine.getStringVal("stage");
		PanCoreGenes.threads = cmdLine.getIntVal("thread");

		//0. Read meta-data and configurations from the config file
		//PanCoreGenes.readSampleList(config);
		PanCoreGenes.readConfig(config);		

		if ("status".equals(stage)){			
			PanCoreGenes.checkStatus();
			System.out.println("Not done yet");			
		}else if ("assemble".equals(stage)){
			PanCoreGenes.assemble();
		}else if ("annotate".equals(stage)){
			PanCoreGenes.annotate();
		}else if ("group".equals(stage)){
			PanCoreGenes.groupGenes();
		}else if ("alignment".equals(stage)){
			PanCoreGenes.alignment();
		}else if ("var".equals(stage)){
			PanCoreGenes.runFreeBayes();	
		}else if ("analyseVar".equals(stage)){
			PanCoreGenes.analyseVar();	
		}else if ("analyseRefVar".equals(stage)){
			for (String arg:args)
				PanCoreGenes.analyseRefVar(arg);	
		}else if ("associate".equals(stage)){
			PanCoreGenes.associate(args[0], args[1]);	
		}
		else if ("genePresence".equals(stage)){
			PanCoreGenes.genePresence();
		}else if ("coreSNPs".equals(stage)){
			PanCoreGenes.getSNPSites();
		}else if ("clean".equals(stage)){
			//run alignment to pan-genes
			//TODO
			System.out.println("Not done yet");
		}else if ("breseq".equals(stage)){
			PanCoreGenes.runBredSeq();
		}else if ("test".equals(stage)){
			System.out.println("Running test");
			PanCoreGenes.runBlast();
		}else{
			Logging.error("Unknown commands "+stage);
		}
	}




	private static void findLGT(String treeFile) throws Exception{				
		PhylogenyTree phylo = PhylogenyTree.readFromFile(treeFile);

		String cdHitFile = "GeneSequences/" + "output.clstr";

		BufferedReader bf = SequenceReader.openFile(cdHitFile);
		String line;

		String geneName = "";
		HashSet<String> geneSampleSet = new HashSet<String>();

		while ((line = bf.readLine())!=null){			
			if (line.startsWith(">")){
				if (geneSampleSet.isEmpty())
					continue;

				//processing
				PhylogenyTree subTree= LGT(phylo,geneSampleSet);
				if (subTree == null){
					System.out.println(geneSampleSet);
				}
				int count = subTree.getLeaves().size();
				if (count == geneSampleSet.size()){
					System.out.println(geneName + "\tgood\t" + count);
				}else if (count < geneSampleSet.size()){
					System.out.println(geneName + "\tbad1\t" + count + " vs " + geneSampleSet.size() + "\t" + subTree + "\t" + geneSampleSet);
				}else{//>
					System.out.println(geneName + "\tbad2\t" + count + " vs " + geneSampleSet.size() + "\t" + subTree + "\t" + geneSampleSet);					
				}

				//clear for the next gene
				geneName = "";
				geneSampleSet.clear();
				continue;
			}			
			String [] toks = line.trim().split(" ");


			if (toks.length < 3){
				bf.close();
				Logging.exit("Not expected 1", -1);
			}


			toks[1] = toks[1].substring(1, toks[1].length()-3);			
			if (toks[2].equals("*")){
				geneName = toks[1];				
			}			

			toks = toks[1].split("_NODE_");			
			geneSampleSet.add(toks[0]);

			//System.out.println(geneName + " " + toks[0] + " " + line);
		}//while
		//processing again

		PhylogenyTree subTree= LGT(phylo,geneSampleSet);
		int count = subTree.getLeaves().size();
		if (count == geneSampleSet.size()){
			System.out.println(geneName + "\tgood\t" + count);
		}else if (count < geneSampleSet.size()){
			System.out.println(geneName + "\tbad1\t" + count + " vs " + geneSampleSet.size());
		}else{//>
			System.out.println(geneName + "\tbad2\t" + count + " vs " + geneSampleSet.size());					
		}
		//phylo.toString()

	}

	public static PhylogenyTree LGT(PhylogenyTree phylo, HashSet<String> geneSampleSet){
		Iterator<PhylogenyTree> iter = phylo.getLeafIterator();
		PhylogenyTree myTree = null;

		//get the common ancestor of all of these taxa in the set
		while (iter.hasNext()){
			PhylogenyTree aTree = iter.next();
			if (geneSampleSet.contains(aTree.getName())){
				if (myTree == null){
					myTree = aTree;
				}else{
					while(!treeContain(myTree,aTree.getName())){
						if (myTree.isRoot())
							return myTree;
						else 
							myTree = myTree.getParent();
					}					
				}
			}
		}
		return myTree;
	}

	public static boolean treeContain(PhylogenyTree tree, String sample){
		Iterator<PhylogenyTree> iter = tree.getLeafIterator();
		while (iter.hasNext()){
			if (sample.equals(iter.next().getName()))
				return true;
		}
		return false;		
	}
}
/*RST*

Design:
input = a config file, defining: samples, folder
commands = 
  status: check status of data, what have been done, what commands available
  assemble: assemle all the genome
  annotate: annotate the assembled genomes
  extract:
  group:
  alignment:
  coreSNPs:



 *RST*/

