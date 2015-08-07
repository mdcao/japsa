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
 * 10/12/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.sim;

import japsa.bio.alignment.ProbFSM.Emission;
import japsa.bio.alignment.ProbFSM.ProbThreeSM;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.util.ArrayList;
import java.util.Random;

/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.sim.testpfsm", 
scriptDesc = "Testing estimation of parameters using a three state finite machine")
public class SimProbFSM {
	//public static SequenceOutputStream datOutGen , datOutEst; 
	public static void main(String[] args) throws Exception{
		/*********************** Setting up script ****************************/
		Deployable annotation = SimProbFSM.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options]", annotation.scriptDesc());		
		/**********************************************************************/		

		cmdLine.addInt("length", 1000, "Length");
		cmdLine.addInt("num", 20, "Number of sequences");
		
		cmdLine.addDouble("miProb", 0.1, "Model probability of insertion");
		cmdLine.addDouble("mdProb", 0.1, "Model probability of deletion");
		cmdLine.addDouble("mmProb", 0.1, "Model probability of mutation");		
		cmdLine.addDouble("meiProb", 0.2, "Model probability of extending insertion");
		cmdLine.addDouble("medProb", 0.2, "Model probability of extending deletion");
		
		
		cmdLine.addDouble("eiProb", 0.025, "Estimate probability of insertion");
		cmdLine.addDouble("edProb", 0.025, "Estimate probability of deletion");
		cmdLine.addDouble("emProb", 0.05, "Estimate probability of mutation");		
		cmdLine.addDouble("eeiProb", 0.1, "Estimate probability of extending insertion");
		cmdLine.addDouble("eedProb", 0.1, "Estimate probability of extending deletion");		
		
		args = cmdLine.stdParseLine_old(args);	
		/**********************************************************************/
		int length   = cmdLine.getIntVal("length");
		int numSeq = cmdLine.getIntVal("num");

		double miProb = cmdLine.getDoubleVal("miProb");
		double mdProb = cmdLine.getDoubleVal("mdProb");
		double mmProb = cmdLine.getDoubleVal("mmProb");
		double meiProb = cmdLine.getDoubleVal("meiProb");
		double medProb = cmdLine.getDoubleVal("medProb");
		
		double eiProb = cmdLine.getDoubleVal("eiProb");
		double edProb = cmdLine.getDoubleVal("edProb");
		double emProb = cmdLine.getDoubleVal("emProb");
		double eeiProb = cmdLine.getDoubleVal("eeiProb");
		double eedProb = cmdLine.getDoubleVal("eedProb");		
		

		Alphabet dna = Alphabet.DNA4();
		Random rnd = new Random(1);
		Sequence mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rnd);			

		ProbThreeSM genDp = new ProbThreeSM(mSeq);

		genDp.getMatState().setTransitionProb(1 - miProb - mdProb, miProb,  mdProb);
		genDp.getMatState().setCopyProb(1 - mmProb);		
		
		genDp.getInsState().setTransitionProb(1 - meiProb, meiProb,  0);
		genDp.getInsState().setCopyProb(1 - mmProb);		
		
		genDp.getDelState().setTransitionProb(1 - medProb, 0, medProb);
		genDp.getDelState().setCopyProb(1 - mmProb);

		
		System.out.println("Length = " + length + " Ins = " + miProb + " Del = " + mdProb + " Mis = " + mmProb);
		
		ArrayList<Sequence> seqs = new ArrayList<Sequence>(numSeq); 
		for (int i = 0; i < numSeq; i++ ){	
			System.out.printf("%3d  ",i);
			Sequence genSeq = genDp.generate(rnd);
			seqs.add(genSeq);
		}
		
		ProbThreeSM eDp = new ProbThreeSM(mSeq);

		eDp.getMatState().setTransitionProb(1 - eiProb - edProb, eiProb,  edProb);
		eDp.getMatState().setCopyProb(1 - emProb);		
		
		eDp.getInsState().setTransitionProb(1 - eeiProb, eeiProb,  0);
		eDp.getInsState().setCopyProb(1 - emProb);		
		
		eDp.getDelState().setTransitionProb(1 - eedProb, 0, eedProb);
		eDp.getDelState().setCopyProb(1 - emProb);
		
		int itNum = 10;//number of iteration
		
		
		for (int x = 0; x < itNum;x++){				
			eDp.resetCount();
			double totCost = 0;
			for (int i = 0; i < numSeq; i++ ){
				System.out.printf("%3d  \n",i);
				Emission retState = eDp.align(seqs.get(i));
				double cost = retState.myCost;
				System.out.println(eDp.updateCount(retState) + " states and " + cost + " bits");
				//predDp.	backward(retState);
				totCost += cost;
			}
			eDp.reEstimate();
			System.out.println("----------------------------------------------\n Total cost = " + totCost);			
			eDp.showProb();
			System.out.println("=============================================");
		}
		genDp.showProb();
		/***************************************************************/
		
	}
}