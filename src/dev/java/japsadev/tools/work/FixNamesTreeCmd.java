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
 * 28/05/2014 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsadev.tools.work;


import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import japsa.bio.phylo.PhylogenyTree;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * Filter a bam filem based on some criteria. Input file in bam format assumed
 * to be sorted and indexed
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.dev.nametaxa", 
	scriptDesc = "Fix names"
	)
public class FixNamesTreeCmd extends CommandLine{
	//CommandLine cmdLine;
	public FixNamesTreeCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("mode", null, "Mode: name, compare, trip",true);		

		addStdHelp();
	}
	public static void main(String [] args) throws Exception{
		FixNamesTreeCmd cmdLine = new FixNamesTreeCmd ();
		args = cmdLine.stdParseLine(args);

		String mode = cmdLine.getStringVal("mode");
		/**********************************************************************/
		if (mode.equals("name")){
			//String base  = cmdLine.getStringVal("base");		
			BufferedReader br = SequenceReader.openFile("index2Names");

			//HashMap<String, String> key12Species = new HashMap<String, String>();
			HashMap<String, String> key2Name = new HashMap<String, String>();
			HashMap<String, String> name2Name = new HashMap<String, String>();
			//HashMap<String, String> key2NewName = new HashMap<String, String>();

			String line;

			while ((line = br.readLine())!=null){
				String [] toks = line.trim().split(" ");			
				String name = toks[3];
				String [] xxxx = name.split("_");
				String key  = (xxxx[0] + "_" + xxxx[2]);			
				key2Name.put(key,toks[1]+ "_" + xxxx[1] + "_" + xxxx[2]);
				name2Name.put(name,toks[1]+ "_" + xxxx[1] + "_" + xxxx[2]);			
			}
			br.close();
			PhylogenyTree tree = PhylogenyTree.readFromFile(args[0]);

			ArrayList<PhylogenyTree> leaves = tree.getLeaves();
			for (PhylogenyTree leaf:leaves){
				String lName = leaf.getName();
				String name = key2Name.get(lName);
				if (name != null){
					Logging.info(lName + " found in key ");
					leaf.setName(name);
				}else{
					name = name2Name.get(lName);
					if (name != null){
						Logging.info(lName + " found in name ");
						leaf.setName(name);
					}else{
						Logging.info(lName + " not found");
					}
				}
			}
			System.out.println(tree.toString() + ";");
		}

		if (mode.equals("trip")){
			PhylogenyTree tree0 = PhylogenyTree.readFromFile(args[0]);
			HashSet<String> leafSet0 = new HashSet<String>();
			ArrayList<PhylogenyTree> leaves0 = tree0.getLeaves();
			for (PhylogenyTree leaf:leaves0){
				String name = leaf.getName().split("_")[0];				
				if (!leafSet0.add(name)){
					PhylogenyTree parent = leaf.getParent();
					if (parent != null){
						int gIndex = leaf.getIndex();

						PhylogenyTree grantParent = parent.getParent();
						if (grantParent != null){
							int cIndex = parent.getIndex();
							PhylogenyTree removed = grantParent.removeGrandChild(cIndex, gIndex);
							if (removed == leaf){
								Logging.info("Yay, removed " + name);
								continue;
							}else{
								Logging.info("Some thing not right " + name);
							}
						}else{
							Logging.warn("GRANT not found of " + name);
						}
					}else{
						Logging.warn("PARENT not found " + name);						
					}					

					//Logging.warn("Once more  " + name);

					if (tree0.getChild(0) == leaf){
						tree0 = tree0.getChild(1);
						tree0.setParent(null);
						Logging.info("Yay, removed with round" + name);
					}else if (tree0.getChild(1) == leaf){
						tree0 = tree0.getChild(0);
						tree0.setParent(null);
						Logging.info("Yay, removed with round" + name);
					}else
						Logging.warn("CAN not remove " + name);
				}	



			}			

			System.out.println(tree0 + ";");
		}		

		if (mode.equals("compare")){
			PhylogenyTree tree0 = PhylogenyTree.readFromFile(args[0]);
			HashSet<String> leafSet0 = new HashSet<String>();

			ArrayList<PhylogenyTree> leaves0 = tree0.getLeaves();			
			for (PhylogenyTree leaf:leaves0){
				leafSet0.add(leaf.getName());
			}

			PhylogenyTree tree1 = PhylogenyTree.readFromFile(args[1]);
			HashSet<String> leafSet1 = new HashSet<String>();

			ArrayList<PhylogenyTree> leaves1 = tree1.getLeaves();

			for (PhylogenyTree leaf:leaves1){				
				leafSet1.add(leaf.getName());

				String name = leaf.getName();
				if (!leafSet0.contains(name)){
					//remove
					PhylogenyTree parent = leaf.getParent();
					if (parent != null){
						int gIndex = leaf.getIndex();

						PhylogenyTree grantParent = parent.getParent();
						if (grantParent != null){
							int cIndex = parent.getIndex();
							PhylogenyTree removed = grantParent.removeGrandChild(cIndex, gIndex);
							if (removed == leaf){
								Logging.info("Yay, removed " + name);
								continue;
							}else{
								Logging.info("Some thing not right " + name);
							}
						}else{
							Logging.warn("GRANT not found of " + name);
						}
					}else{
						Logging.warn("PARENT not found " + name);						
					}					

					//Logging.warn("Once more  " + name);

					if (tree1.getChild(0) == leaf){
						tree1 = tree1.getChild(1);
						tree1.setParent(null);
						Logging.info("Yay, removed with round" + name);
					}else if (tree1.getChild(1) == leaf){
						tree1 = tree1.getChild(0);
						tree1.setParent(null);
						Logging.info("Yay, removed with round" + name);
					}else
						Logging.warn("CAN not remove " + name);
				}					

			}


			for (PhylogenyTree leaf:leaves0){				
				String name = leaf.getName();
				if (!leafSet1.contains(name)){
					//remove
					PhylogenyTree parent = leaf.getParent();
					if (parent != null){
						int gIndex = leaf.getIndex();

						PhylogenyTree grantParent = parent.getParent();
						if (grantParent != null){
							int cIndex = parent.getIndex();
							PhylogenyTree removed = grantParent.removeGrandChild(cIndex, gIndex);
							if (removed == leaf){
								Logging.info("Yay, removed " + name);
								continue;
							}else{
								Logging.info("Some thing not right " + name);
							}
						}else{
							Logging.warn("GRANT not found of " + name);
						}
					}else{
						Logging.warn("PARENT not found " + name);						
					}					

					//Logging.warn("Once more  " + name);
					if (tree0.getChild(0) == leaf){
						tree0 = tree0.getChild(1);
						tree0.setParent(null);
						Logging.info("Yay, removed with round" + name);
					}else if (tree0.getChild(1) == leaf){
						tree0 = tree0.getChild(0);
						tree0.setParent(null);
						Logging.info("Yay, removed with round" + name);
					}else
						Logging.warn("CAN not remove " + name);
				}						

			}//
			System.out.println(tree0 + ";");
			System.out.println(tree1 + ";");
			//System.out.println(tree1);
		}
	}
}

