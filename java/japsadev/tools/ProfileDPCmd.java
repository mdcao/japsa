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

/****************************************************************************
 *                           Revision History                                
 * 5 Sep 2016 - Minh Duc Cao: Started                                 
 *  
 ****************************************************************************/
package japsadev.tools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import japsa.bio.alignment.ProfileDP;
import japsa.bio.alignment.ProfileDP.EmissionState;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */

/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.profileDP",
	scriptDesc = "Using a 1-state machine for alignment"
	)

public class ProfileDPCmd extends CommandLine{	
	public ProfileDPCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addInt("length", 20, "Length");
		addDouble("iProb", 0.10, "Probability of insertion");
		addDouble("dProb", 0.10, "Probability of deletion");
		addDouble("mProb", 0.10, "Probability of mutation");
		
		
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new ProfileDPCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		int length   = cmdLine.getIntVal("length");

		double iProb = cmdLine.getDoubleVal("iProb");
		double dProb = cmdLine.getDoubleVal("dProb");
		double mProb = cmdLine.getDoubleVal("mProb");

		Alphabet dna = Alphabet.DNA4();
		Random rnd = new Random(1);
		Sequence seq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rnd);

		//SequenceOutputStream out = SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("out"));
		//seq.print(out);
		/*****************************************************/
		//Sequence seq = SequenceReader.getReader(args[0]).nextSequence(dna);
		//ProfileDP dp = new ProfileDP(seq, 20, seq.length() - 20);
		ProfileDP genDp = new ProfileDP(seq, -1, seq.length() *2);

		//datOutGen = SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("prefix") + "gen.dat");
		//datOutEst = SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("prefix") + "est.dat");

		genDp.setTransitionProbability(1 - iProb - dProb, iProb,  dProb);
		genDp.setMatchProbability(1 - mProb);		

		ProfileDP dp = new ProfileDP(seq, -1, seq.length() *2);

		System.out.println("Length = " + length + " Ins = " + iProb + " Del = " + dProb + " Mis = " + mProb);
		System.out.printf("%8.4f %8.4f %8.4f\n",dp.getMatCost(),dp.getInsCost(),dp.getDelCost());
		System.out.printf("%8.4f %8.4f\n",dp.getMatchCost(),dp.getMisMatchCost());

		int numSeq = 20;
		ArrayList<Sequence> seqs = new ArrayList<Sequence>(numSeq); 
		for (int i = 0; i < numSeq; i++ ){	
			System.out.printf("%3d  ",i);
			Sequence genSeq = genDp.generate(1, rnd);
			seqs.add(genSeq);
			//genSeq.print(out);			
			/***************************************************************
			//Emission alignScore = dp.align(genSeq);			
			Emission alignScore = genDp.align(genSeq);
			System.out.println(alignScore.score + "  " + genSeq.length() * 2);	

			IntArray iP = new IntArray();
			IntArray iS = new IntArray();

			Emission tmp = alignScore;
			do{
				iP.add(tmp.profilePos);
				iS.add(tmp.seqPos);
				tmp = tmp.bwdState;
			}while (tmp != null);

			for (int x = iP.size()-1; x> 0; x--){
				int p = iP.get(x) + 1;
				int s = iS.get(x) + 1;

				if (iP.get(x) == iP.get(x-1)){
					datOutEst.print("I " + p + " " + s + " " + genSeq.charAt(s) + "\n");
				}else if (iS.get(x) == iS.get(x-1)){
					datOutEst.print("D " + p + " " + s + " " + seq.charAt(p) + "\n");
				}else if (seq.getBase(p) == genSeq.getBase(s)){
					datOutEst.print("= " + p + " " + s + " " + seq.charAt(p) + "\n");
				}else
					datOutEst.print("X " + p + " " + s + " " + seq.charAt(p) + " " + genSeq.charAt(s) + "\n");
			}
			datOutEst.print("EST: " + alignScore.countMG + " " +alignScore.countMB + " " + alignScore.countIns + " " + alignScore.countDel + " " + alignScore.score + "\n");
			/***************************************************************/
		}		
		//out.close();
		/***************************************************************
		datOutEst.close();
		datOutGen.close();
		/***************************************************************/
		for (int x = 0; x < 5;x++){
			int countIns = 0, countDel = 0, countMG = 0, countMB = 0;
			for (int i = 0; i < numSeq; i++ ){
				System.out.printf("%3d  ",i);
				EmissionState retState = dp.align(seqs.get(i));
				countIns += retState.getCountIns();
				countDel += retState.getCountDel();
				countMG += retState.getCountMG();
				countMB += retState.getCountMB();
			}

			double sum = 3.0 + countMG + countMB + countIns + countDel;
			double insP = (countIns + 1.0) /sum;
			double delP = (countDel + 1.0) /sum;
			double matP = (countMG + countMB + 1.0) /sum;
			double matchP = (countMG + 1.0) / (countMG + countMB + 2.0);
			double misMatchP = 1 - matchP;
			System.out.printf("Total: %3d %3d %3d %3d %8.4f %8.4f %8.4f\n", countMG, countMB, countIns, countDel,  insP, delP,  misMatchP);
			dp.setTransitionProbability(matP, insP, delP);
			dp.setMatchProbability(matchP);
		}
		/***************************************************************/
		
		
	}	
}
/*RST*



 
  
  
  
  
  
  
  
  
*RST*/