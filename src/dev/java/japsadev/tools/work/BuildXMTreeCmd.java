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


import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Filter a bam filem based on some criteria. Input file in bam format assumed
 * to be sorted and indexed
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.dev.treexm", 
	scriptDesc = "Build tree of XM"
	)
public class BuildXMTreeCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(BuildXMTreeCmd.class);

	public BuildXMTreeCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("output", null, "output",true);		

		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		BuildXMTreeCmd cmdLine = new BuildXMTreeCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		//String base  = cmdLine.getStringVal("base");		
		BufferedReader br = SequenceReader.openFile("index2Names");

		HashMap<String, Integer> id2Index = new HashMap<String, Integer>();
		ArrayList<String> idList = new ArrayList<String>();
		ArrayList<String> nameList = new ArrayList<String>();

		LOG.info("point 1");
		String line;
		int noTaxa = 0;
		while ((line = br.readLine())!=null){
			String [] toks = line.trim().split(" ");
			String ID  = toks[0];
			
			String name = toks[3];
			String [] xxxx = name.split("_");

			idList.add(ID);
			nameList.add(xxxx[0] + "_" + xxxx[2]);

			id2Index.put(ID,noTaxa);
			noTaxa ++;			
		}		
		br.close();

		int [] lengths = new int[noTaxa];
		double [] myInfo = new double[noTaxa];

		double [][] disMeasure1 = new double[noTaxa][];
		double [][] disMeasure2 = new double[noTaxa][];

		LOG.info("point 2");
		int index = 0;		
		while (index < noTaxa){
			File file = new File("rdistance."+index+".out");
			if (!file.exists()){
				break;//no more
			}

			LOG.info("Read " + index);

			BufferedReader reader = new BufferedReader(new FileReader(file));
			line = reader.readLine();
			String [] toks = line.trim().split("\t");
			if (toks.length != 3){
				break;
			}
			int myLength = Integer.parseInt(toks[1]);
			lengths[index] = myLength;
			myInfo [index] = Double.parseDouble(toks[2]);			

			boolean good =  true;


			disMeasure1[index] = new double[index];
			disMeasure2[index] = new double[index];

			int count = 0;
			while ((line = reader.readLine())!=null){
				toks = line.trim().split("\t");
				if (toks.length != 4){
					good = false;
					break;
				}
				int mateIndex = id2Index.get(toks[0]);
				int mateLength = Integer.parseInt(toks[1]);

				double  ij = Double.parseDouble(toks[2]);
				double  ji = Double.parseDouble(toks[3]);

				disMeasure1[index][mateIndex] = (ji + ij) /(myInfo[index] + myInfo[mateIndex]);
				disMeasure2[index][mateIndex] = (ji * myLength + ij * mateLength) /(myInfo[index] * myLength + myInfo[mateIndex] * mateLength);								

				count ++;
			}		
			reader.close();

			if (count < index)
				good = false;

			if (!good)
				break;

			index++;
		}

		System.out.println(index);

		PrintStream out = 
			new PrintStream(new BufferedOutputStream(
				new FileOutputStream(
					"dis1"+cmdLine.getStringVal("output"))));


		out.println(" " + index);
		for (int s = 0; s < index; s++) {
			out.printf("%-12s ", nameList.get(s));
			for (int x = 0; x < index; x++) {
				if (x < s)
					out.printf(" %10f ", disMeasure1[s][x]);
				else if (x==s)
					out.printf(" %10f ", 0.0);
				else
					out.printf(" %10f ", disMeasure1[x][s]);
			}
			out.println();			
		}	
		out.close();

		out = new PrintStream(new BufferedOutputStream(
			new FileOutputStream(
				"dis2"+cmdLine.getStringVal("output"))));


		out.println(" " + index);
		for (int s = 0; s < index; s++) {
			out.printf("%-12s ", nameList.get(s));
			for (int x = 0; x < index; x++) {
				if (x < s)
					out.printf(" %10f ", disMeasure2[s][x]);
				else if (x==s)
					out.printf(" %10f ", 0.0);
				else
					out.printf(" %10f ", disMeasure2[x][s]);
			}
			out.println();			
		}	
		out.close();
	}	
}
