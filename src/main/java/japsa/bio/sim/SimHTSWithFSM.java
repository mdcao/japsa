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


import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import japsa.bio.alignment.ProbFSM.ProbThreeSM;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.DoubleArray;
import japsa.util.JapsaMath;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;


/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.sim.htsSim", 
scriptDesc = "Simulation of HTS sequencing with a probabistic finite (3) state machine")
public class SimHTSWithFSM {
	//public static SequenceOutputStream datOutGen , datOutEst; 
	public static void main(String[] args) throws Exception{
		/*********************** Setting up script ****************************/
		Deployable annotation = SimHTSWithFSM.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options]", annotation.scriptDesc());		
		/**********************************************************************/		
		cmdLine.addStdInputFile();		
		cmdLine.addString("output", null, "Name of output fastq file", true);
		
		cmdLine.addInt("lm", 4000, "Read length means");
		cmdLine.addInt("ls", 1500, "Read length standard deviation");		
		
		cmdLine.addDouble("mm", 0.15, "Probability of mismatches");		
		
		cmdLine.addDouble("io", 0.1, "Probability of insertion opening");
		cmdLine.addDouble("do", 0.1, "Probability of deletion opening");
		
		cmdLine.addDouble("ie", 0.2, "Probability of insertion extension");
		cmdLine.addDouble("de", 0.2, "Probability of deletion extension");
		//Note that this model forbids going directly from an insertion to deletion
		//match = 1 - mismatch - del - ins
		
		cmdLine.addInt("seed", 0, "Random seed, <=0 for current time");		
		cmdLine.addDouble("cov", 50, "Coverage");		
		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/		
		
		String output = cmdLine.getStringVal("output");
		String inFile = cmdLine.getStringVal("input");		
		
		double mmProb = cmdLine.getDoubleVal("mm");//
		
		double mdProb = cmdLine.getDoubleVal("do");//m->d		
		double miProb = cmdLine.getDoubleVal("io");//m->i
		
				
		double iiProb = cmdLine.getDoubleVal("ie");
		double ddProb = cmdLine.getDoubleVal("de");
		
		
		int lm   = cmdLine.getIntVal("lm");
		int ls   = cmdLine.getIntVal("ls");
		int seed = cmdLine.getIntVal("seed");
		double cov = cmdLine.getDoubleVal("cov");		
		
		Random rand = seed> 0?( new Random(seed)): (new Random());
		
		Alphabet alphabet = Alphabet.DNA();
				
		
		DoubleArray dArray = new DoubleArray();
		double sumLength = 0;		

		//Reading all chromosomes in
		SequenceReader reader = SequenceReader.getReader(inFile);
		Sequence seq;
		
		ArrayList<Sequence> genome = new ArrayList<Sequence>(); 
		while ((seq = reader.nextSequence(alphabet))!=null){
			genome.add(seq);
			dArray.add(seq.length());
			sumLength += seq.length();			
		}
		reader.close();
		
		//get the accumulative (to select chr)
		double [] frac = dArray.toArray();
		
		frac[0] /= sumLength;		
		for (int i = 1; i < frac.length;i++){
			frac[i] /= sumLength;
			frac[i] += frac[i-1];
		}
		frac[frac.length - 1] = 1.0;
		
		
		ProbThreeSM genDp = new ProbThreeSM(null);

		genDp.getMatState().setTransitionProb(1 - miProb - mdProb, miProb,  mdProb);
		genDp.getMatState().setCopyProb(1 - mmProb);		
		
		genDp.getInsState().setTransitionProb(1 - iiProb, iiProb,  0);
		genDp.getInsState().setCopyProb(1 - mmProb);		
		
		genDp.getDelState().setTransitionProb(1 - ddProb, 0, ddProb);
		genDp.getDelState().setCopyProb(1 - mmProb);
		
		
		//A rough estimation of quality
		double qual = (miProb > mdProb)?1-miProb:1-mdProb;
		
		qual *= (1-mmProb);
		System.out.println(qual);
		//convert to fred
		int qualPhred = (int) Math.floor(JapsaMath.prob2phred(1.0 - qual) + 0.5);
		//essentially, this is round
		
		//wont get this, but just incase
		if (qualPhred > 40)
			qualPhred = 40;
		
		if (qualPhred < 0)
			qualPhred = 0;
		
		char qualPhredChar = (char) (33 + qualPhred); 
				
		System.out.println(qualPhredChar);
		
		SequenceOutputStream outStream = 
				 SequenceOutputStream.makeOutputStream(output);
		
		double numBase = 0;
		cov *= sumLength;
		
		long readID = 0;
		while (numBase < cov){			
			
			int length = (int) (rand.nextGaussian() * ls + lm);
			if (length < 50)
				continue;
			
			double toss = rand.nextDouble();
			int index = 0;
			
			while (frac[index] < toss )
				index ++;
			//assert frac[index] >= toss
			if (genome.get(index).length() < length)
				continue;
			
			int start = rand.nextInt(genome.get(index).length() - length);
			genDp.setModelSequence(genome.get(index).subSequence(start, start + length));
			seq = genDp.generate(rand);
			
			readID ++;
			outStream.print("@R");
			outStream.print(readID);
			outStream.print(sep);
			outStream.print(genome.get(index).getName().replace(sep, '_'));
			outStream.print(sep);
			outStream.print(start+1);
			outStream.print(sep);
			if (rand.nextDouble() < .5){
				outStream.print('0');
				seq = Alphabet.DNA.complement(seq);
			}else{
				outStream.print('1');
			}
			outStream.print(sep);
			outStream.print('0');
			outStream.print(sep);
			outStream.print('0');
			outStream.print(sep);
			outStream.print(seq.length());
			outStream.print('\n');
				
			for (int i = 0; i < seq.length();i++){
				outStream.print(seq.charAt(i));
			}
			outStream.print("\n+\n");
			for (int i = 0; i < seq.length();i++){
				outStream.print(qualPhredChar);
			}
			outStream.print('\n');			
			
			numBase += seq.length();
		}
		outStream.close();						
	}
	
	static final char sep = '#';
	static final String sepSTR = "#";
	
	static int relax = 20;
	static void eval(String bamFile){
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFile));	
						
		SAMRecordIterator samIter = samReader.iterator();
		
		int countRead = 0;
		int TP = 0, FP = 0, FN =0 , dup = 0;
		while (samIter.hasNext()){
			countRead ++;
			SAMRecord samRecord = samIter.next();
			
			String readName = samRecord.getReadName();
			String[] toks = readName.split(sepSTR);
			
			
			if (samRecord.getReadUnmappedFlag()){
				FP ++;
				continue;
			}
			
			if (!samRecord.getReferenceName().equals(toks[1])){
				FP ++;
				continue;
			}

			int rightPos = (samRecord.getReadPairedFlag() && (!samRecord.getFirstOfPairFlag()))?
					Integer.parseInt(toks[4]):Integer.parseInt(toks[2]);
					
			int alignPos = samRecord.getReferencePositionAtReadPosition(1);
			
			if (Math.abs(alignPos - rightPos) < relax)
				TP ++;
			else 
				FP ++;					
		}
		System.out.println("Accuracy " + (TP * 100.0 / (TP+FP)) + "%");
	}
}