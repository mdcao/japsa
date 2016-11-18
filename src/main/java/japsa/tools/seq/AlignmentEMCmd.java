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

package japsa.tools.seq;

import japsa.bio.alignment.ProbFSM;
import japsa.bio.alignment.ProbFSM.Emission;
import japsa.bio.alignment.ProbFSM.ProbThreeSM;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.seq.emalign", 
	scriptDesc = "Get the best alignment of 2 sequences using Expectation-Maximisation on Finite State Machine"
)

public class AlignmentEMCmd extends CommandLine{	
	public AlignmentEMCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options] seq1 seq2");
		setDesc(annotation.scriptDesc());
		
		addInt("iteration",5,"Number of iteration");			
		
		addStdHelp();		
	} 
	
	//public static SequenceOutputStream datOutGen , datOutEst; 
	public static void main(String[] args) throws Exception{				 		
		CommandLine cmdLine = new AlignmentEMCmd();		
		args = cmdLine.stdParseLine(args);	
		
		int itNum = cmdLine.getIntVal("iteration");
		
		
		Alphabet dna = Alphabet.DNA();
		if (args.length <2){
			System.err.println("Two sequence files are required\n" +
			cmdLine.errorString());
			System.exit(-1);
		}
		
		SequenceReader readFile = SequenceReader.getReader(args[0]);
		SequenceReader barcodeFile = SequenceReader.getReader(args[1]);
		
		//Sequence mSeq = SequenceReader.getReader(args[0]).nextSequence(dna);
		//Sequence sSeq = SequenceReader.getReader(args[1]).nextSequence(dna);

		SequenceOutputStream out = SequenceOutputStream.makeOutputStream("-");
		Sequence readSeq = null, barcodeSeq = null;
		while ((readSeq = readFile.nextSequence(dna)) != null){		
			while ((barcodeSeq = barcodeFile.nextSequence(dna)) != null){
				//forward run
				ProbFSM eDp = new ProbThreeSM(readSeq);
				Emission retState = null;				
				
				for (int x = 0; x < itNum;x++){				
					eDp.resetCount();					
					retState = eDp.alignGenerative(barcodeSeq);
					double cost = retState.myCost;
					System.out.println(eDp.updateCount(retState) + " states and " + cost + " bits " + barcodeSeq.length() + "bp"  );			
					eDp.reEstimate();
					System.out.println("----------------------------------------------\n Total cost = " + cost);			
					eDp.showProb();
					System.out.println("=============================================");
				}				
				
				eDp.printAlignment(retState, barcodeSeq, out);
				
				barcodeSeq = Alphabet.DNA.complement(barcodeSeq);
				eDp = new ProbThreeSM(readSeq);
				retState = null;				
				
				for (int x = 0; x < itNum;x++){				
					eDp.resetCount();					
					retState = eDp.alignGenerative(barcodeSeq);
					double cost = retState.myCost;
					System.out.println(eDp.updateCount(retState) + " states and " + cost + " bits " + barcodeSeq.length() + "bp"  );			
					eDp.reEstimate();
					System.out.println("----------------------------------------------\n Total cost = " + cost);			
					eDp.showProb();
					System.out.println("=============================================");
				}
				eDp.printAlignment(retState, barcodeSeq, out);
			}	
		}
		
		out.close();	
		
		/***************************************************************/	
	}
}

/*RST*
----------------------------------------------
*jsa.seq.emalign* Align two sequences using EM
----------------------------------------------


<usage> 

*RST*/