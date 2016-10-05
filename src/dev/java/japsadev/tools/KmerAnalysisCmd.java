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
 * 07/09/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.tools;

import java.io.IOException;
import java.util.ArrayList;

import japsa.bio.alignment.ProbFSM.Emission;
import japsa.bio.alignment.ProbFSM.ProbOneSM;
import japsa.bio.alignment.ProfileDP;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.seq.nanopore.Fast5DetailReader;
import japsa.seq.nanopore.Fast5DetailReader.DetectionEvents;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.dev.np", scriptDesc = "Kmer Analysis of nanopore")
public class KmerAnalysisCmd extends CommandLine{	
	public KmerAnalysisCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("type", "kmer",  "kmer, event");
		addString("input",null,"Input file");
		addInt("kmer",10,"K");
		addString("fastq",null,"File containing reads");		
		
		
		addStdHelp();		
	} 

	/**
	 * @param args
	 * @throws Exception 
	 * @throws OutOfMemoryError 
	 */
	public static void main(String[] args) throws OutOfMemoryError, Exception {

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new KmerAnalysisCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
				
		String input = cmdLine.getStringVal("input");					
		int k = cmdLine.getIntVal("kmer");
		String type = cmdLine.getStringVal("type");	
		String fq   = cmdLine.getStringVal("fastq");
		
		
		long time = System.currentTimeMillis();
		if (type.equals("kmer"))
			kmerAnalysis(input,k);
		else if (type.equals("events"))
			eventAlignment(input,k);
		else if (type.equals("params")){
			System.out.println("Param with DP");
			ArrayList<Sequence> seqs = SequenceReader.readAll(input, Alphabet.DNA4());
			Sequence seq = seqs.get(0);//take only the first
			seqs = SequenceReader.readAll(fq, Alphabet.DNA4());
			
			
			for (int i = 0; i < seqs.size();i++){				
				Sequence rseq = seqs.get(i);
				ProfileDP dp = new ProfileDP(seq, -10 ,-10);
				
				dp.setTransitionProbability(0.88, 0.05 , 0.07);
				dp.setMatchProbability(0.86);
				
				//ProbThreeSM s3SM = new ProbThreeSM(seq);								
				ProbOneSM s1SM = new ProbOneSM(seq);
				
				System.out.println(rseq.length() * seq.length());
				
				for (int x = 0; x < 10; x ++){				
					
					double cost;
					ProfileDP.EmissionState state = dp.align(rseq);
					cost = state.getScore();
					System.out.println(rseq.getName() + " " + x + ":" + cost + " " + cost*1.0/seq.length() + " " + state.getCountMB() + " mismatches " + state.getCountMG() + " matches " + state.getCountDel() + " del " + state.getCountIns() + " ins");
					System.out.println("DPProfile cost " + state.getScore());
					dp.setTransitionProbability(state.getCountMB() + state.getCountMG(), state.getCountIns(), state.getCountDel());
					dp.setMatchProbability(state.getCountMG() * 1.0 / (state.getCountMB() + state.getCountMG()));
					
					long timeNow = System.currentTimeMillis();
					System.out.println("Time " + (timeNow - time));
					time = timeNow;
					System.out.println("---------------------------------------------------------");
					
					
					System.gc();					
					Emission retState;
					
					//s3SM.resetCount();					
					//retState = s3SM.align(rseq);
					//cost = retState.myCost;
					//System.out.println("Three-State cost " + cost);
					//s3SM.updateCount(retState);
					//s3SM.reEstimate();
								
					//s3SM.showProb();
					//timeNow = System.currentTimeMillis();
					//System.out.println("Time " + (timeNow - time));
					//time = timeNow;
					//System.out.println("---------------------------------------------------------");
					
					
					s1SM.resetCount();					
					retState = s1SM.align(rseq);
					cost = retState.myCost;
					System.out.println("One-State cost " + cost);
					s1SM.updateCount(retState);
					s1SM.reEstimate();
								
					s1SM.showProb();
					timeNow = System.currentTimeMillis();
					System.out.println("Time " + (timeNow - time));
					time = timeNow;
					System.out.println("---------------------------------------------------------");
					
				}//for x				
				System.out.println("=====================================================================");
			}//for i
		}		
	}
	
	public static void eventAlignment(String input, int k) throws OutOfMemoryError, Exception {
		Fast5DetailReader reader = new Fast5DetailReader(input);
		reader.readData();
		
		DetectionEvents events = reader.getEvents();
		double [] mean = events.getMean();
		int min = 1, max = 50;
		double eps = 0.7;
		
		for (int x = min;x <= max; x++){
			int lastBad = max;
			//for (int i = max; i < mean.length;i++ ){
			for (int i = 1000; i < 5000;i++ ){
				if (Math.abs(mean[i] - mean[i-x]) > eps)
					lastBad = i;
				
				if (i - lastBad > 3){
					System.out.println(i+"  " + x);
				}				
			}		
		}				
		
		reader.close();
	}
	
	public static void kmerAnalysis(String input, int k) throws IOException {
		ArrayList<Sequence> seqs = SequenceReader.readAll(input, Alphabet.DNA4());
		StringBuilder sb = new StringBuilder(k);
		sb.setLength(k);
		
		for (Sequence seq:seqs){			 
			for (int i = 0; i< seq.length() - k; i ++){
				//get the kmer
				for (int j = 0; j < k; j++){
					sb.setCharAt(j, seq.charAt(i+j));					
				}				
				System.out.println(sb.toString()+"\t"+i);				
			}
		}
	}
}
