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

package japsa.tools.bio.phylo;

import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.util.CommandLine;
import japsa.util.JapsaMath;
import japsa.util.deploy.Deployable;

import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;


@Deployable(scriptName = "jsa.phylo.xmas",
            scriptDesc = "Generate a distance matrix from aligned sequences")
public class XMasCmd  extends CommandLine{	
	public XMasCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addStdInputFile();		
		addString("output", "output", "Name of the file for output (distances in phylip format)");
		addBoolean("adapt", false, "Use adaptive");
		
		
		addStdHelp();		
	} 
	public static boolean adapt = false;
	public static void main(String[] args) throws Exception {		 		
		CommandLine cmdLine = new XMasCmd();
		args = cmdLine.stdParseLine(args);
		
		/**********************************************************************/
		FastaReader sin = new FastaReader(new FileInputStream(cmdLine.getStringVal("input")));
		ArrayList<Sequence> seqs = new ArrayList<Sequence>(100);
		Sequence seq;
		while ((seq = sin.nextSequence(Alphabet.DNA4())) != null){
			seqs.add(seq);
		}
		/**************************************************/
		sin.close();
		adapt = cmdLine.getBooleanVal("adapt");

		double[][] mtx = buildMtx(seqs);
		
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

	public static double[][] buildMtx(ArrayList<Sequence> seqs) {		

		double[][] mtx = new double[seqs.size()][seqs.size()];
		double myCosts[] = new double[seqs.size()];

		int [] counts = new int[4];
		double [] f = new double[4];
		for (int i = seqs.size() - 1; i >= 0; i--) {
			//1. first compute the markov 
			/*****************************************************************/
			counts[3]  = counts[2] = counts[1] = counts[0] = 1;
			Sequence seq = seqs.get(i); 
			for (int j = seq.length() - 1; j >=0; j--){
				counts[seq.getBase(j)] ++;				
			}
			double sum = 4.0 + seq.length();
			f[0] = -JapsaMath.log2(counts[0]/sum);
			f[1] = -JapsaMath.log2(counts[1]/sum);
			f[2] = -JapsaMath.log2(counts[2]/sum);
			f[3] = -JapsaMath.log2(counts[3]/sum);

			for (int j = seq.length() - 1; j >=0; j--){
				myCosts[i] += f[seq.getBase(j)];
			}

			myCosts[i] /= seq.length();
			/*****************************************************************/			
		}

		for (int i = seqs.size() - 1; i > 0; i--) {		
			mtx[i][i] = 0.0;
			for (int j = 0; j < i; j++) {

				double r;
				if (adapt)
					r = condEntropyAdaptive(seqs.get(j).toBytes(), seqs.get(i).toBytes());
				else
					r = condEntropy(seqs.get(j).toBytes(), seqs.get(i).toBytes());

				mtx[i][j] = mtx[j][i] = -10.0 
				* JapsaMath.log2((myCosts[i] + myCosts[j] - r)
						/ ((myCosts[i] + myCosts[j])));
			}
		}
		return mtx;
	}

	static double condEntropyAdaptive(byte[] s1, byte[] s2) {
		// return I(x|y) + I(y|x)
		double r = 0;

		//q2[a,b] = Pr(x1=b|x2=a); sum2[a] = Pr(x2=a) -> predict x1
		//q1[a,b] = Pr(x2=b|x1=a); sum1[a] = Pr(x1=a) -> predict x2
		
		double[][] q1 = 
			  { { 1.0, 1.0, 1.0, 1.0 },
				{ 1.0, 1.0, 1.0, 1.0 },
				{ 1.0, 1.0, 1.0, 1.0 }, 
				{ 1.0, 1.0, 1.0, 1.0 } };

		double[][] q2 = 
			  { { 1.0, 1.0, 1.0, 1.0 }, 
				{ 1.0, 1.0, 1.0, 1.0 },
				{ 1.0, 1.0, 1.0, 1.0 }, 
				{ 1.0, 1.0, 1.0, 1.0 } };

		double[] sum1 =   {4, 4, 4, 4 }, sum2 = { 4, 4, 4, 4 };
				
		for (int i = 0; i < s1.length; i++) {			
			r -=  JapsaMath.log2(q2[s2[i]][s1[i]] / sum2[s2[i]])
				+ JapsaMath.log2(q1[s1[i]][s2[i]] / sum1[s1[i]]);
			
			q1[s1[i]][s2[i]] += 1.0;
			sum1[s1[i]] += 1.0;				

			q2[s2[i]][s1[i]] += 1.0;
			sum2[s2[i]] += 1.0;			
		}
		return r / s1.length;
	}
	
	
	static double condEntropy(byte[] s1, byte[] s2) {
		// return I(x|y) + I(y|x)
		//double r = 0;
		// First pass
		double[][] q1 = { { 1.0, 1.0, 1.0, 1.0 },
				{ 1.0, 1.0, 1.0, 1.0 },
				{ 1.0, 1.0, 1.0, 1.0 }, 
				{ 1.0, 1.0, 1.0, 1.0 } };

		double[][] q2 = { { 1.0, 1.0, 1.0, 1.0 }, 
				{ 1.0, 1.0, 1.0, 1.0 },
				{ 1.0, 1.0, 1.0, 1.0 }, 
				{ 1.0, 1.0, 1.0, 1.0 } };

		double[] sum1 =   {4, 4, 4, 4 }, sum2 = { 4, 4, 4, 4 };


		// first pass		
		for (int i = 0; i < s1.length; i++) {			
			q1[s1[i]][s2[i]] += 1.0;
			sum1[s1[i]] += 1.0;				

			q2[s2[i]][s1[i]] += 1.0;
			sum2[s2[i]] += 1.0;			
		}
		/***********************************************************/
		for (int i = 0; i< q1.length; i++){			
			for (int j = 0; j< q1.length; j++){
				q1[i][j] /= sum1[i];
				q2[i][j] /= sum2[i];

				q1[i][j] = - JapsaMath.log2(q1[i][j]);
				q2[i][j] = - JapsaMath.log2(q2[i][j]);				
			}			
		}
		/***********************************************************/
		// second pass
		double r1 = 0, r2 = 0;

		for (int i = 0; i < s1.length; i++) {
			r2 += q1[s1[i]][s2[i]];
			r1 += q2[s2[i]][s1[i]];
			//r2 -= JapsaMath.log2(q1[s1[i]][s2[i]]);
			//r1 -= JapsaMath.log2(q2[s2[i]][s1[i]]);			
		}
		/***********************************************************/

		
		return (r1 + r2) / s1.length;
	}
	
	static double condEntropyF(byte[] s1, byte[] s2) {
		// return I(x|y) + I(y|x)
		//double r = 0;
		// First pass
		double[][] q1 = { { 1.0, 1.0, 1.0, 1.0 },
				{ 1.0, 1.0, 1.0, 1.0 },
				{ 1.0, 1.0, 1.0, 1.0 }, 
				{ 1.0, 1.0, 1.0, 1.0 } };

		double[][] q2 = { { 1.0, 1.0, 1.0, 1.0 }, 
				{ 1.0, 1.0, 1.0, 1.0 },
				{ 1.0, 1.0, 1.0, 1.0 }, 
				{ 1.0, 1.0, 1.0, 1.0 } };

		double[] sum1 =   {4, 4, 4, 4 }, sum2 = { 4, 4, 4, 4 };


		// first pass		
		for (int i = 0; i < s1.length; i++) {			
			q1[s1[i]][s2[i]] += 1.0;
			sum1[s1[i]] += 1.0;				

			q2[s2[i]][s1[i]] += 1.0;
			sum2[s2[i]] += 1.0;			
		}
		/***********************************************************/
		double r = 0;
		//double [][] l = new double[q1.length][q1.length];
		for (int i = 0; i< q1.length; i++){			
			for (int j = 0; j< q1.length; j++){
				//q1[i][j] /= sum1[i];
				//q2[i][j] /= sum2[i];
				r -= JapsaMath.log2(q1[i][j]/sum1[i]) * (q1[i][j] );
				r -= JapsaMath.log2(q2[j][i]/sum2[j]) * (q2[j][i] );				
			}			
		}
		
		return r / s1.length;
	}

	/********************************************************************/

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
