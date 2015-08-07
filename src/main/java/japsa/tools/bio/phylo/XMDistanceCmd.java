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
import japsa.tools.xm.ExpertModelCmd;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;
import japsa.xm.ExpertModel;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;


@Deployable(scriptName = "jsa.phylo.xmdist",
scriptDesc = "Generate a distance matrix from genomes (potentially not alignable")
public class XMDistanceCmd  extends CommandLine{	
	public XMDistanceCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addStdInputFile();		
		addString("output", "output", "Name of the file for output (distances in phylip format)");


		addInt("hashSize", 11, "Hash size");
		addInt("context", 15, "Length of the context");
		addInt("limit", 200, "Expert Limit");
		addInt("thread", 1, "Number of threads");
		addDouble("threshold", 0.15, "Listen threshold");
		addInt("chance", 20, "Chances");
		addBoolean("binaryHash", false, "Use binary hash or not");
		addString("offsetType", "counts",
			"Way of update offset/palindrome expert: possible value count, subs");
		addBoolean("optimise", false,
			"Running in optimise mode, just report the entropy,recommended for long sequence");
		addInt("checkPoint", 1000000, "Frequency of check point");
		addString("hashType", "hash",
			"Type of Hash table: hash=hashtable, sft=SuffixTree,sfa = SuffixArray");
		addBoolean("selfRep", true,
			"Propose experts from the sequence to compressed?");	

		addStdHelp();		
	} 
	//public static boolean adapt = false;

	static double [] resultSingle;
	static double [][] resultBG;
	static ArrayList<Sequence> seqs;

	public static void main(String[] args) throws Exception {		 		
		CommandLine cmdLine = new XMDistanceCmd();
		args = cmdLine.stdParseLine(args);
		

		String input = cmdLine.getStringVal("input");
		int thread = cmdLine.getIntVal("thread");

		System.out.println(input);
		FastaReader sin = new FastaReader(input);

		seqs = new ArrayList<Sequence>(100);
		Sequence seq;
		while ((seq = sin.nextSequence(Alphabet.DNA4())) != null){
			seqs.add(seq);
		}
		/**************************************************/
		sin.close();		
		resultBG = new double[seqs.size()][seqs.size()];		
		resultSingle =new double[seqs.size()];


		ExecutorService executor = Executors.newFixedThreadPool(thread);
		//Sequence [] mS = new Sequence[2];
		//Sequence [] sS = new Sequence[1];
		//ExpertModel eModel = ExpertModelDriver.getExpertModel(cmdLine);

		resultBG[0][0] = 0.0;
		for (int i = 0; i < seqs.size();i++){
			executor.execute(new CompressSingle(cmdLine,i));			
			//sS[0] = seqs.get(i);
			//resultSingle[i] = eModel.encode_optimise(sS);

			//System.out.println("Single " +  i + " : " + resultSingle[i]);
		}

		for (int i = 1; i < seqs.size();i++){
			synchronized (resultBG){
				resultBG[i][i] = 0.0;
			}
			for (int j = 0; j < i; j++){				
				executor.execute(new CompressSingle(cmdLine,i,j));				
				//mS[0] = seqs.get(i);
				//mS[1] = seqs.get(j);
				//double e_ij = eModel.encode_optimise(mS);

				//System.out.println("BG " +  i +"," + j + " : " + e_ij);

				//mS[1] = seqs.get(i);
				//mS[0] = seqs.get(j);
				//double e_ji = eModel.encode_optimise(mS);
				//System.out.println("BG " +  j +"," + i + " : " + e_ji);


				//resultBG[i][j] = resultBG[j][i] = (e_ij + e_ji)/(resultSingle[i] + resultSingle[j]);				
				//System.out.println(i + "  -  " + j + " : " + resultBG[i][j]);
			}
		}
		executor.shutdown();

		boolean finished = executor.awaitTermination(3, TimeUnit.DAYS);


		double [][] mtx = new double[seqs.size()] [seqs.size()];

		Logging.info("ALL DONE " + finished);
		for (int i = 0; i < seqs.size();i++){
			mtx[i][i] = 0;
			for (int j = i+1; j < seqs.size(); j++){
				mtx[i][j] = mtx[j][i] = (resultBG[i][j] + resultBG[j][i])/(resultSingle[i] + resultSingle[j]);	
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

	static class CompressSingle implements Runnable {
		ExpertModel eModel;		

		int index1;
		int index2;		

		private CompressSingle(CommandLine cmdLine) throws Exception{
			eModel = ExpertModelCmd.getExpertModel(cmdLine);
		}

		public CompressSingle(CommandLine cmdLine, int i1) throws Exception{
			this(cmdLine);
			this.index1 = i1;
			this.index2 = -1;
		}

		public CompressSingle(CommandLine cmdLine, int i1, int i2) throws Exception{
			this(cmdLine);
			this.index1 = i1;
			this.index2 = i2;
		}


		/* (non-Javadoc)
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run() {
			try {
				if (index2 < 0){
					Logging.info("Thread Single " + index1 + " started");
					Sequence [] mS = new Sequence[1];
					mS[0] = seqs.get(index1);
					double score = eModel.encode_optimise(mS);
					synchronized(resultSingle){
						resultSingle[index1] = score;
					}
					Logging.info("Thread Single " + index1 + " done!");
				}else{
					Logging.info("Thread GB " + index1 + " - " + index2 + " started");
					Sequence [] mS = new Sequence[2];
					mS[0] = seqs.get(index1);
					mS[1] = seqs.get(index2);
					double e_ij = eModel.encode_optimise(mS);
					synchronized(resultBG){
						resultBG[index1][index2] = e_ij; 
					}
					Logging.info("Thread GB " + index1 + " - " + index2 + " done");


					Logging.info("Thread GB2 " + index2 + " - " + index1 + " started");
					mS[0] = seqs.get(index2);
					mS[1] = seqs.get(index1);					
					double e_ji = eModel.encode_optimise(mS);

					synchronized(resultBG){
						resultBG[index2][index1] = e_ji; 
					}
					Logging.info("Thread GB2 " + index2 + " - " + index1 + " done");
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
}
