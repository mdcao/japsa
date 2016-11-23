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

package japsadev.tools;

import java.io.IOException;
import java.util.ArrayList;

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
	scriptName = "jsa.dev.barcodeMDC", 
	scriptDesc = "Barcode"
)

public class NanoporeBarcodeCmd extends CommandLine{	
	public NanoporeBarcodeCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options] seq1 seq2");
		setDesc(annotation.scriptDesc());
		
		addInt("iteration",5,"Number of iteration");			
		
		addStdHelp();		
	} 
	
	//public static SequenceOutputStream datOutGen , datOutEst; 
	public static void main(String[] args) throws Exception{				 		
		CommandLine cmdLine = new NanoporeBarcodeCmd();		
		args = cmdLine.stdParseLine(args);	
		
		int itNum = cmdLine.getIntVal("iteration");
		
		
		Alphabet dna = Alphabet.DNA();
		if (args.length <2){
			System.err.println("Two sequence files are required\n" +
			cmdLine.errorString());
			System.exit(-1);
		}
		
		SequenceReader readFile = SequenceReader.getReader(args[0]);		
		ArrayList<Sequence>  barcodeSeqs = SequenceReader.readAll(args[1], dna);
		
		
		//Sequence mSeq = SequenceReader.getReader(args[0]).nextSequence(dna);
		//Sequence sSeq = SequenceReader.getReader(args[1]).nextSequence(dna);

		SequenceOutputStream out = SequenceOutputStream.makeOutputStream("-");
		Sequence readSeq = null;
		while ((readSeq = readFile.nextSequence(dna)) != null){
			if (readSeq.length() < 240)
				continue;
			double score = Double.MAX_VALUE;
			
			int length = readSeq.length();
			
			Sequence head = readSeq.subSequence(0, 120);
			head.setName("head");
			Sequence tail = readSeq.subSequence(length - 120, length);
			tail.setName("tail");
			
			for (int i = 0; i < barcodeSeqs.size()/2;i++){
				Sequence fcode = barcodeSeqs.get(i*2);
				Sequence rcode = barcodeSeqs.get(i*2 + 1);
				
				Sequence fcodeRev = Alphabet.DNA.complement(fcode);
				Sequence rcodeRev = Alphabet.DNA.complement(rcode);
				
				fcodeRev.setName("fRev");
				rcodeRev.setName("rRev");
				
				out.print(head.getName() + " vs " + fcode.getName());
				out.println();
				double score1 = align(head, fcode, itNum, out);
				
				out.print(head.getName() + " vs " + rcode.getName());
				out.println();
				double score2 = align(head, rcode, itNum, out);
				
				out.print(tail.getName() + " vs " + fcode.getName());
				out.println();
				double score3 = align(tail, fcode, itNum, out);
				
				out.print(tail.getName() + " vs " + rcode.getName());
				out.println();
				double score4 = align(tail, rcode, itNum, out);				
				
				out.print(head.getName() + " vs " + fcode.getName());
				out.println();
				double score5 = align(head, fcodeRev, itNum, out);
				
				out.print(head.getName() + " vs " + rcode.getName());
				out.println();
				double score6 = align(head, rcodeRev, itNum, out);
				
				out.print(tail.getName() + " vs " + fcode.getName());
				out.println();
				double score7 = align(tail, fcodeRev, itNum, out);
				
				out.print(tail.getName() + " vs " + rcode.getName());
				out.println();
				
				double score8 = align(tail, rcodeRev, itNum, out);
				out.print("####################################################################");
				out.println();
				
				out.print("XXX " + readSeq.getName() + " " + fcode.getName() + 
						"\t" + score1 +
						"\t" + score2 + 
						"\t" + score3 +
						"\t" + score4 +
						"\t" + score5 +
						"\t" + score6 +
						"\t" + score7 +
						"\t" + score8 +
						"\n");
				
				if (score1 < score)
					score = score1;
				if (score2 < score)
					score = score2;
				if (score3 < score)
					score = score3;
				if (score4 < score)
					score = score4;
				if (score5 < score)
					score = score5;
				if (score6 < score)
					score = score6;
				if (score7 < score)
					score = score7;
				if (score8 < score)
					score = score8;
			}	
			out.print("YYY " + readSeq.getName() + " " + score + "\n");
		}		
		out.close();
		
		/***************************************************************/	
	}
	public  static double align(Sequence r, Sequence b, int itNum, SequenceOutputStream out) throws IOException{
		ProbFSM eDp = new ProbThreeSM(r);
		
		Emission retState = null;
		double score = Double.MAX_VALUE;
		double diff = 0.1;
		
		for (int x = 0; x < itNum;x++){				
			eDp.resetCount();					
			retState = eDp.alignGenerative(b);
			double cost = retState.myCost;
			if (cost < score)
				score = cost;
			out.print(eDp.updateCount(retState) + " states and " + cost + " bits " + b.length() + "bp"  );
			out.println();
			eDp.reEstimate();
			out.print("----------------------------------------------\n Total cost = " + cost);
			out.println();
			eDp.showProb();
			out.print("=============================================");
			out.println();
		}				
		
		eDp.printAlignment(retState, b, out);
		return score;
	}
}

/*RST*
----------------------------------------------
*jsa.seq.emalign* Align two sequences using EM
----------------------------------------------


<usage> 

*RST*/