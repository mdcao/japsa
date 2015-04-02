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

package japsa.bio.alignment;

import japsa.bio.alignment.ProbFSM.Emission;
import japsa.bio.alignment.ProbFSM.ProbThreeSM;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.seq.emalign", 
scriptDesc = "Get the best alignment of 2 sequences using Expectation-Maximisation on Finite State Machine")
public class AlignmentEM {
	//public static SequenceOutputStream datOutGen , datOutEst; 
	public static void main(String[] args) throws Exception{
		/*********************** Setting up script ****************************/
		Deployable annotation = AlignmentEM.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options] <seq1> <seq2>", annotation.scriptDesc());		
		/**********************************************************************/		
		
		args = cmdLine.stdParseLine(args);	
		/**********************************************************************/

		Alphabet dna = Alphabet.DNA();
		if (args.length <2){
			System.err.println("Two sequence files are required" +
			cmdLine.usage());
			System.exit(-1);
		}
		
		Sequence mSeq = SequenceReader.getReader(args[0]).nextSequence(dna);
		Sequence sSeq = SequenceReader.getReader(args[1]).nextSequence(dna);

		ProbThreeSM eDp = new ProbThreeSM(mSeq);
		
		int itNum = 10;//number of iteration
		
		
		for (int x = 0; x < itNum;x++){				
			eDp.resetCount();
			
			Emission retState = eDp.align(sSeq);
			double cost = retState.myCost;
			System.out.println(eDp.updateCount(retState) + " states and " + cost + " bits " + sSeq.length() + "bp"  );			
			eDp.reEstimate();
			System.out.println("----------------------------------------------\n Total cost = " + cost);			
			eDp.showProb();
			System.out.println("=============================================");
		}
		
		/***************************************************************/
		
	}
}