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

/**************************     REVISION HISTORY    **************************
 * File: Deploy.java
 * 15/11/2013 - Minh Duc Cao: Created
 *
 ****************************************************************************/

package japsa.util.deploy;


import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * This class is used to deploy tools: create a makefile to generate scripts
 * 
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 */
public class GenDocs {
	/**
	 * 
	 * @param masterScript
	 * @throws IOException
	 */
	public static void genDocs(ArrayList<Object> toolList) 
		throws IOException{		

		//First check what tools has documentation
		ArrayList<CommandLine> cmdTools = new ArrayList<CommandLine>();
		ArrayList<ArrayList<String>> cmdDocs = new ArrayList<ArrayList<String>>();
		HashSet<String> cmdSet = new HashSet<String>();

		for (Object obj : toolList) {
			if (obj instanceof CommandLine){
				CommandLine cmdObj = (CommandLine) obj;
				Class<?> tool =  obj.getClass();				

				Deployable annotation = (Deployable) tool.getAnnotation(Deployable.class);
				if (annotation == null)
					continue;//for

				String scriptName = annotation.scriptName();
				ArrayList<String> docs = new ArrayList<String>();

				File sourceFile = new File("src/main/java/" + tool.getCanonicalName().replace('.', '/') + ".java");
				if (!sourceFile.isFile())
					continue;		

				BufferedReader bf = new BufferedReader (new FileReader(sourceFile));
				String line = "";
				boolean addDocs = false;
				while ((line = bf.readLine())!= null){
					if (line.trim().startsWith("/*RST*")){
						//Start reading
						addDocs = true;
						continue;
					}					
					if (line.trim().startsWith("*RST*/")){
						break;
					}

					if (addDocs)
						docs.add(line.replaceAll("\\*#/", "\\*/"));					
				}

				bf.close();
				if (docs.size() == 0){
					System.out.println("No document found for " + scriptName);
					continue;
				}

				cmdDocs.add(docs);
				cmdTools.add(cmdObj);				
				cmdSet.add(scriptName);
			}
		}

		//Generating index file
		System.out.println("Generating documentation :");
		String docDir = "docs/source/tools";

		PrintStream toolIndex = new PrintStream(new FileOutputStream(docDir + "/index.rst"));

		toolIndex.println("=============\nList of tools\n=============\n");
		toolIndex.println("This chapter presents the list of tools provided by Japsa.\n"
			+ " We are in the process of documenting 40+ tools, so stay tuned.\n");

		toolIndex.println(".. toctree::\n   :maxdepth: 1\n");		


		for (int i = 0; i < cmdTools.size();i++){
			CommandLine cmd = cmdTools.get(i);
			ArrayList<String> docs = cmdDocs.get(i);			
			Deployable annotation = (Deployable) cmd.getClass().getAnnotation(Deployable.class);			
			String scriptName = annotation.scriptName();

			SequenceOutputStream outOS = SequenceOutputStream.makeOutputStream(docDir + "/" + scriptName + ".rst");

			for (String doc:docs){
				if (doc.trim().equals("<usage>")){
					outOS.print("~~~~~~~~\nSynopsis\n~~~~~~~~\n\n");
					outOS.print("*" + scriptName + "*: " + annotation.scriptDesc()+"\n\n");
					outOS.print("~~~~~\nUsage\n~~~~~\n::\n\n   "+cmd.usage() + "\n\n");
					outOS.print("~~~~~~~\nOptions\n~~~~~~~\n"+cmd.options().replace("_","\\_")+"\n\n");						

					if (annotation.seeAlso().length() > 0){
						boolean seeAlso = false;

						StringBuilder sb = new StringBuilder();

						String [] toks = annotation.seeAlso().split(",");
						for (String tok:toks){
							tok = tok.trim();
							if (tok.length()<=0) 
								continue;
							
							
							if (cmdSet.contains(tok)){
								sb.append(".. _" + tok+": " + tok+".html\n");
								tok = tok + "_";
							}

							if (!seeAlso){
								outOS.print("~~~~~~~~\nSee also\n~~~~~~~~\n\n" + tok);
								seeAlso = true;
							}else{
								outOS.print(", " + tok);
							}								

															
						}
						if (sb.length() > 0)
							outOS.print("\n\n" + sb.toString() + "\n");							
					}
					if (annotation.citation().length() > 0){
						outOS.print("~~~~~~~~\nCitation\n~~~~~~~~\n\n" + annotation.citation() + "\n\n");							
					}						
				}else
					outOS.print(doc);

				outOS.println();
			}

			outOS.close();
			toolIndex.println("   " + scriptName + ".rst");
			System.out.println(" " + docDir + "/" + scriptName + ".rst created");


		}
		toolIndex.close();	
	}

	public static void main(String[] args) throws IOException {
		genDocs(Deploy.tools);
	}
}
