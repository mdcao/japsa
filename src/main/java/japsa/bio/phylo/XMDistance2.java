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
 * 08/10/2012 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsa.bio.phylo;

import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import japsa.xm.ExpertModel;
import japsa.xm.ExpertModelDriver;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;


@Deployable(scriptName = "jsa.phylo.distance2",
            scriptDesc = "Generate a distance matrix from genomes (potentially not alignable")
public class XMDistance2 {
	//public static boolean adapt = false;
	public static void main(String[] args) throws Exception {
		/*********************** Setting up script ****************************/
		Deployable annotation = XMDistance2.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options]");
		/**********************************************************************/		
		cmdLine.addStdInputFile();		
		cmdLine.addString("output", "output", "Name of the file for output (distances in phylip format)");
		
						
		cmdLine.addInt("hashSize", 11, "Hash size");
		cmdLine.addInt("context", 15, "Length of the context");
		cmdLine.addInt("limit", 200, "Expert Limit");
		cmdLine.addDouble("threshold", 0.15, "Listen threshold");
		cmdLine.addInt("chance", 20, "Chances");
		cmdLine.addBoolean("binaryHash", false, "Use binary hash or not");
		cmdLine.addString("offsetType", "counts",
						"Way of update offset/palindrome expert: possible value count, subs");
		cmdLine.addBoolean("optimise", false,
						"Running in optimise mode, just report the entropy,recommended for long sequence");
		cmdLine.addInt("checkPoint", 1000000, "Frequency of check point");
		cmdLine.addString("hashType", "hash",
						"Type of Hash table: hash=hashtable, sft=SuffixTree,sfa = SuffixArray");
		cmdLine.addBoolean("selfRep", true,
					"Propose experts from the sequence to compressed?");
	
		args = cmdLine.stdParseLine(args);		
		/**********************************************************************/
		
		String input = cmdLine.getStringVal("input");
		
		System.out.println(input);
		FastaReader sin = new FastaReader(input);
		
		ArrayList<Sequence> seqs = new ArrayList<Sequence>(100);
		Sequence seq;
		while ((seq = sin.nextSequence(Alphabet.DNA4())) != null){
			seqs.add(seq);
		}
		/**************************************************/
		sin.close();
		

		double[][] mtx = new double[seqs.size()][seqs.size()];		
		
		Sequence [] mS = new Sequence[2];
		Sequence [] sS = new Sequence[1];
		ExpertModel eModel = ExpertModelDriver.getExpertModel(cmdLine);
		
		mtx[0][0] = 0.0;
		double [] e_i =new double[seqs.size()];
		
		for (int i = 0; i < seqs.size();i++){
			sS[0] = seqs.get(i);
			e_i[i] = eModel.encode_optimise(sS);
			
			System.out.println("Single " +  i + " : " + e_i[i]);
		}
		
		for (int i = 1; i < seqs.size();i++){
			mtx[i][i] = 0.0;
			for (int j = 0; j < i; j++){				
				mS[0] = seqs.get(i);
				mS[1] = seqs.get(j);
				double e_ij = eModel.encode_optimise(mS);
				
				System.out.println("BG " +  i +"," + j + " : " + e_ij);
				
				mS[1] = seqs.get(i);
				mS[0] = seqs.get(j);
				double e_ji = eModel.encode_optimise(mS);
				System.out.println("BG " +  j +"," + i + " : " + e_ji);
				
				
				mtx[i][j] = mtx[j][i] = (e_ij + e_ji)/(e_i[i] + e_i[j]);				
				System.out.println(i + "  -  " + j + " : " + mtx[i][j]);
			}
		}
		
		
		
		PrintStream out = 
			new PrintStream(new BufferedOutputStream(new FileOutputStream(
					cmdLine.getStringVal("output"))));		
		printMtx(seqs, mtx, out);		
		/***********************************************
		BufferedOutputStream out = 
			new BufferedOutputStream(new FileOutputStream(
					cmdLine.getStringVal("output")));		
		StringBuilder sb = new StringBuilder();
		sb.append(' ');
		sb.append(seqHash.size());
		sb.append('\n');
		Formatter formatter = new Formatter(sb);

		for (int s = 0; s < seqHash.size(); s++) {
			formatter.format("%-12s ", seqHash.get(s));
			for (int x = 0; x < seqHash.size(); x++) {
				formatter.format("%8.4f", mtx[s][x]);
			}
			sb.append("\n");							
		}
		out.write(sb.toString().getBytes());
		/***********************************************/
		out.close();
	}

	
	public static void printMtx(ArrayList<Sequence> dnaSeqs, double[][] mtx,
			PrintStream out) {		
		out.println(" " + dnaSeqs.size());
		for (int s = 0; s < dnaSeqs.size(); s++) {
			out.printf("%-12s ", dnaSeqs.get(s).getName());
			for (int x = 0; x < dnaSeqs.size(); x++) {
				out.printf(" %10f ", mtx[s][x]);
			}
			out.println();			
		}
	}
}
