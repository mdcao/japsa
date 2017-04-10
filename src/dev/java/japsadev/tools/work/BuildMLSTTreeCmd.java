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
import java.io.IOException;
import java.util.HashMap;

import japsa.bio.amra.MLSTyping;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;
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
	scriptName = "jsa.dev.treest", 
	scriptDesc = "Build tree of MLST profiles"
	)
public class BuildMLSTTreeCmd extends CommandLine{
	//CommandLine cmdLine;
	public BuildMLSTTreeCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("base", null, "base");		
		
		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		BuildMLSTTreeCmd cmdLine = new BuildMLSTTreeCmd ();
		args = cmdLine.stdParseLine(args);
		
		/**********************************************************************/
		String base  = cmdLine.getStringVal("base");		
		BufferedReader br = SequenceReader.openFile("index2Names");
		
		HashMap<String, MLSTyping> mlstMap = new HashMap<String, MLSTyping>();
		HashMap<String, SequenceBuilder> profileMap = new HashMap<String, SequenceBuilder>();
		
		String line;		
		while ( (line = br.readLine())!=null){
			String [] toks = line.trim().split(" ");
			
			MLSTyping mlst = mlstMap.get(toks[1]);
			if (mlst == null){
				mlst = new MLSTyping(base + "/" + toks[1]);
				mlstMap.put(toks[1], mlst);
			}
			//mlst			
			
			String profile = toks[1] + "#" + toks[4];
			if (!profileMap.containsKey(profile)){
				Logging.info("Profile " + profile);
				SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA4(), 1000, toks[3]);
				String [] genesNames = toks[4].split("\\|");
				for (int x =0; x < 7;x++){
					int sepPos = genesNames[x].lastIndexOf('_');
					int alleleNo = Integer.parseInt(genesNames[x].substring(sepPos+1));					
					int alleleIndex = mlst.alleleNo2AlleleIndex(x, alleleNo);
					Sequence seq = mlst.alleles(x).get(alleleIndex);
					if (!genesNames[x].equals(seq.getName())){
						Logging.exit("Error at " + genesNames[x] + " vs " + seq.getName(), 1);
					}
					Logging.info("Found " + seq.getName());
					sb.append(seq);
				}
				profileMap.put(profile, sb);
			}			
		}
		
		br.close();
		
		HashMap<String, SequenceOutputStream> outMap = new HashMap<String, SequenceOutputStream>();
		for (String key:profileMap.keySet()){
			String [] toks = key.split("#");
			String species = toks[0];
			SequenceOutputStream myOut = outMap.get(species);
			if (myOut == null){
				myOut = SequenceOutputStream.makeOutputStream("ST_"+species + ".fasta");
				outMap.put(species, myOut);
			}
			
			profileMap.get(key).writeFasta(myOut);
		}//for
		
		for (SequenceOutputStream myOut:outMap.values()){
			myOut.close();
		}		
	}	
}
