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
 * 10/05/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsadev.tools.misc;

import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.JapsaFileFormat;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.util.Random;



/**
 * Simulate genomes with short tandem repeat variations. Used in STRViper
 * @author minhduc
 * 
 */
public class STRVariationSimulation{
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/		 
		String scriptName = "str.simATGenome";
		String desc = "Simulate genomes with short tandem variations.\n";// + Sequence.AUTHOR;		
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [params]");
		/**********************************************************************/
		cmdLine.addString("input",  null, "Name of input file",true);		
		cmdLine.addString("output", null, "Prefix output file",true);		
		cmdLine.addString("myCost",  null, "Name of  file containing myCost values");		
		
		cmdLine.addBoolean("help", false, "Display usage message and exit");

		/**********************************************************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")){
			System.out.println(desc + cmdLine.usageMessage());			
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usageMessage());
			System.exit(-1);
		}	
		/**********************************************************************/
		
		
		String inFile = cmdLine.getStringVal("input");
		String scoreFile = cmdLine.getStringVal("myCost");
		String outFile = cmdLine.getStringVal("output");		
		
		//First read in the reference sequences		
		BufferedReader sIn = new BufferedReader(new FileReader(scoreFile));
		
		//BioCompFileFormat fileFormat =  
		//	new BioCompFileFormat(SequenceReader.openFile(inFile));
		
		JapsaFileFormat reader = new JapsaFileFormat(inFile);
		
		
		JapsaAnnotation anno = reader.readAnnotation();
		//fileFormat.getAnnotationIterator().next();
		Sequence seq = anno.getSequence();
		//fileFormat.getSequenceIterator().next();
		
		
		
		double SNPs = 1.5 * 4.9*1000000/119146348, indels = 1.5 * 810467.0/119146348;
		int numSeqs = 3;//maximum distance						
		int length = seq.length();		
		
		Random rnd = new Random(1);
		
		byte[][] seqByte = new byte[numSeqs][length + length / 5];		
		byte[][] seqByteNV = new byte[numSeqs][length + length / 5];		
		
		// make plenty of space
		System.out.println(length);

		JapsaAnnotation[] annos = new JapsaAnnotation[numSeqs];		
		JapsaAnnotation[] annosNV = new JapsaAnnotation[numSeqs];
		
		for (int i = 0; i < numSeqs; i++) {
			annos[i] = new JapsaAnnotation();
			annos[i].addDescription("Simulated genome " + (i+1) + " SNP = " + ((i + 1) * SNPs / numSeqs) + " Indels = " + ((i + 1) * indels / numSeqs));
			
			annosNV[i] = new JapsaAnnotation();
			annosNV[i].addDescription("Simulated genome " + (i+1) + " SNP = " + ((i + 1) * SNPs / numSeqs) + " Indels = " + ((i + 1) * indels / numSeqs));
			
		}
		
		int [] currentInx = new int[numSeqs];
		int [] currentInxNV = new int[numSeqs];
		
		int [] numSNPs = new int[numSeqs];
		int [] numIndels = new int[numSeqs];
		//currentIdx[i] = 0;
		
		int featureIdx = 0;
		JapsaFeature currentFeature = anno.getFeature(featureIdx);
		int index = 0;
		for (; index < seq.length();){
			if (currentFeature == null || index < currentFeature.getStart()){
				for (int seqIdx = 0; seqIdx < numSeqs; seqIdx++){
					char nucleotide = Character.toUpperCase(seq.charAt(index));
					if (nucleotide != 'A' && nucleotide != 'C' && nucleotide != 'G' && nucleotide != 'T'){
						//Generate a random
						seqByte[seqIdx][currentInx[seqIdx]] = (byte) (rnd.nextInt(4));						
						seqByteNV[seqIdx][currentInxNV[seqIdx]] = seqByte[seqIdx][currentInx[seqIdx]];//copy the previous symbol
						
						currentInx[seqIdx] ++;
						currentInxNV[seqIdx] ++;
						
					}else{//
						double val = rnd.nextDouble();				
						if (val <= SNPs * (seqIdx+1) / numSeqs){//SNPs
							//An SNP
							//Generate a random number between 0-2, then plus 1 and plus the index of the previous char
							//to avoid generating the same nucleotide
							seqByte[seqIdx][currentInx[seqIdx]] = (byte) ( (1 + rnd.nextInt(3) + Alphabet.DNA4().char2int(nucleotide)) % 4);
							seqByteNV[seqIdx][currentInxNV[seqIdx]] = seqByte[seqIdx][currentInx[seqIdx]];//copy the previous symbol
							
							currentInx[seqIdx] ++;
							currentInxNV[seqIdx] ++;
							
							numSNPs[seqIdx] ++;
							System.out.println(seqIdx+": substitution " + (index +1) + " at " + currentInx[seqIdx]);
						}else if (val <= SNPs * (seqIdx+1) / numSeqs + indels * (seqIdx + 1) / numSeqs){//indel
							numIndels[seqIdx] ++;
							val = rnd.nextDouble();
							if (val > 0.5 && currentInx[seqIdx] > 1){
								System.out.println(seqIdx+": deletion " + (index +1) + " at " + currentInx[seqIdx]);
								//A deletion
								//currentInx[seqIdx] --;								
							}else{//an insertion
								seqByte[seqIdx][currentInx[seqIdx]] = (byte)  (rnd.nextInt(4));
								seqByteNV[seqIdx][currentInxNV[seqIdx]] = seqByte[seqIdx][currentInx[seqIdx]];//copy the previous symbol
								
								currentInx[seqIdx] ++;
								currentInxNV[seqIdx] ++;							
								
								
								System.out.println(seqIdx+": insertion " + (index +1) + " at " + currentInx[seqIdx]);
								seqByte[seqIdx][currentInx[seqIdx]] = (byte) (Alphabet.DNA4().char2int(nucleotide));
								seqByteNV[seqIdx][currentInxNV[seqIdx]] = seqByte[seqIdx][currentInx[seqIdx]];//copy the previous symbol
								currentInx[seqIdx] ++;
								currentInxNV[seqIdx] ++;
							}//insert vs delete								
						}else{//Direct copy
							seqByte[seqIdx][currentInx[seqIdx]] = (byte) (Alphabet.DNA4().char2int(nucleotide));
							seqByteNV[seqIdx][currentInxNV[seqIdx]] = seqByte[seqIdx][currentInx[seqIdx]];//copy the previous symbol
							
							currentInx[seqIdx] ++;
							currentInxNV[seqIdx] ++;
							//System.out.println(seqIdx+": copy " + index + " at " + currentInx[seqIdx]);
						}
					}//if
				}//for				
				index ++;	
			}else{
				//start of a STR
				//Read the myCost
				String scoreLine = sIn.readLine();
				String [] toks = scoreLine.split("\\t");
				int period = Integer.parseInt(toks[2]);
				double varScore = Double.parseDouble(toks[4]);				
				
				for (int seqIdx = 0; seqIdx < numSeqs; seqIdx++){
					//scale such as 
					int diff = (int) ((seqIdx +1)* varScore * (rnd.nextDouble() + .5) * 1.5);
				
					if (diff <= 0){//include varScore <=0
						//no change
						diff = 0;						
					}else if (diff < currentFeature.getLength() / period && rnd.nextDouble() <= 0.5){
						//contraction
						diff = -diff;					
					}//else expantion
					
					System.out.println(seqIdx+":D = " +diff + " P = " + period + " L = " + currentFeature.getLength() + " " + varScore);
					//assert: diff = 0: no variation
					//diff < 0 : contraction
					//diff > 0 : expantion
					JapsaFeature aFeature =  currentFeature.cloneFeature();
					aFeature.setStart(currentInx[seqIdx]);
					aFeature.setEnd(currentFeature.getStart() + currentFeature.getLength() -1 + diff * period);
					aFeature.addDesc("@DIF:"+diff);
					aFeature.addDesc("@VR:"+period);
					
					JapsaFeature aFeatureNV =  currentFeature.cloneFeature();
					aFeatureNV.setStart(currentInxNV[seqIdx]);
					aFeatureNV.setEnd(currentFeature.getEnd());
					aFeatureNV.addDesc("@DIF:0");
					aFeatureNV.addDesc("@VR:"+period);
					
					annos[seqIdx].add(aFeature);
					annosNV[seqIdx].add(aFeatureNV);
					
					int i = currentFeature.getLength() - 1;
					int j = aFeature.getLength() - 1;
					
					for (;i >=0 && j >=0;){
						char nucleotide = Character.toUpperCase(seq.charAt(index + i));				
						seqByte[seqIdx][currentInx[seqIdx] + j] = (byte) (Alphabet.DNA4().char2int(nucleotide));
						if (seqByte[seqIdx][currentInx[seqIdx] + j] < 0) 
							seqByte[seqIdx][currentInx[seqIdx] + j] = (byte) rnd.nextInt(4);
						i--;j--;
					}					
					while(j >= 0){
						//expansion
						if (i < 0)
							i = period - 1;
						char nucleotide = Character.toUpperCase(seq.charAt(index + i));						
						seqByte[seqIdx][currentInx[seqIdx] + j] = (byte) (Alphabet.DNA4().char2int(nucleotide));
						if (seqByte[seqIdx][currentInx[seqIdx] + j] < 0) 
							seqByte[seqIdx][currentInx[seqIdx] + j] = (byte) rnd.nextInt(4);
						
						i--;j--;
					}
					currentInx[seqIdx] += aFeature.getLength();
					
					for (i = currentFeature.getLength() - 1;i >=0 ;){
						char nucleotide = Character.toUpperCase(seq.charAt(index + i));				
						seqByteNV[seqIdx][currentInxNV[seqIdx] + i] = (byte) (Alphabet.DNA4().char2int(nucleotide));
						if (seqByteNV[seqIdx][currentInxNV[seqIdx] + i] < 0) 
							seqByteNV[seqIdx][currentInxNV[seqIdx] + i] = (byte) rnd.nextInt(4);
						i--;
					}
					currentInxNV[seqIdx] += aFeatureNV.getLength();					
					
				}//for
				index += currentFeature.getLength();
				
				featureIdx ++;
				if (featureIdx < anno.numFeatures())
					currentFeature = anno.getFeature(featureIdx);
				else
					currentFeature = null;				
			}
		}	
		sIn.close();
		reader.close();
		
		//System.out.println("Insert reps done");
		for (int ni = 0; ni < numSeqs; ni++) {			
			Sequence nSeq = new Sequence(Alphabet.DNA4(), seqByte[ni], currentInx[ni], seq.getName());			
			
			SequenceOutputStream out = SequenceOutputStream.makeOutputStream(outFile + (ni+1) 	+ ".bio");
			
			annos[ni].addDescription(numSNPs[ni] + " SNPs " + numIndels[ni] + " indels ");
			JapsaAnnotation.write(nSeq, annos[ni], out);
			out.close();

			out = SequenceOutputStream.makeOutputStream(outFile + (ni+1) + ".fas");
			nSeq.writeFasta(out);
			out.close();			
			System.out.println("Sequence " + (ni+1) + " : " + numSNPs[ni] + " SNPs " + numIndels[ni] + " indels ");
			////////////////////////////////////////////////////////////////////////////////////////////////////////
			nSeq = new Sequence(Alphabet.DNA4(), seqByteNV[ni], currentInxNV[ni], seq.getName());			
			
			out = SequenceOutputStream.makeOutputStream(outFile + (ni+1) 	+ ".NV.bio");
			annosNV[ni].addDescription(numSNPs[ni] + " SNPs " + numIndels[ni] + " indels ");
			JapsaAnnotation.write(nSeq, annosNV[ni], out);
			out.close();

			out = SequenceOutputStream.makeOutputStream(outFile + (ni+1) + ".NV.fas");
			nSeq.writeFasta(out);
			out.close();			
			System.out.println("Sequence " + (ni+1) + " : " + numSNPs[ni] + " SNPs " + numIndels[ni] + " indels ");
		}
		System.out.println("Write  done");
		

	}

}
