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
 * 10/11/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsadev.tools.misc;



import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.JapsaFileFormat;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.util.Arrays;
import java.util.Random;




/**
 * Simulate the genome of Ler-0 with STR variations
 * @author minhduc
 * 
 */
public class SimulateLer0Genome{
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/		 
		String scriptName = "japsa.seq.simLer0";
		String desc = "Simulate a Ler-0 genome from an existing genome.\n";// + Sequence.AUTHOR;		
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [params]");
		/**********************************************************************/
		cmdLine.addString("input",  null, "Name of input file (genome of Col-0)",true);
		cmdLine.addString("var",    null, "Name of variations input file (genome of Col-0)",true);		
		cmdLine.addString("output", null, "Prefix output file",true);
		
		cmdLine.addInt("seed", 0, "Random seed, 0 for a random seed");
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

		String inFile = cmdLine.getStringVal("input");	//input in str.combio format	
		String outFile = cmdLine.getStringVal("output");
		String scoreFile = cmdLine.getStringVal("var");
		
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(outFile);
		
		//First read in the reference sequences		
		
		int seed = cmdLine.getIntVal("seed");
		Random rnd;
		if (seed <= 0)
			rnd = new Random();
		else
			rnd = new Random(seed);	
		
		
		BufferedReader sIn = new BufferedReader(new FileReader(scoreFile));		
		//BioCompFileFormat fileFormat =  
		//	new BioCompFileFormat(SequenceReader.openFile(inFile));
	//	
	//	Iterator<JapsaAnnotation> annoItr = fileFormat.getAnnotationIterator();
	//	Iterator<Sequence> seqItr = fileFormat.getSequenceIterator();
		
		JapsaFileFormat reader = new  JapsaFileFormat(inFile);
		
		
		for (int xx = 1; xx <= 5; xx++){
			JapsaAnnotation anno = reader.readAnnotation();			
			Sequence seq = anno.getSequence();
			
			//		seqItr.next();
			//JapsaAnnotation anno = annoItr.next();
			
			JapsaAnnotation newAnno = new JapsaAnnotation();			
			int fIndex = 0;
			JapsaFeature currentFeature = anno.getFeature(fIndex);
			
			//byte [] mSeq = japsa.seq.getSequence();
			byte [] newSeq = new byte[seq.length() + seq.length() / 2];
			int  [] indexes = new int[newSeq.length];//mapping back the new japsa.seq -> old
			
			Arrays.fill(indexes, -1);
			
			//note posSrc'es are 0-based index			
			BufferedReader in = SequenceReader.openFile(xx+".sdi");
			String line;
			int mPos = 0, newPos = 0;
			
			while (mPos < seq.length()){//
				//assert nextSTR == null || nextSTR.getStart >= mPos
				line = in.readLine();
				if (line == null){
					if (currentFeature != null){
						System.err.println("I am not done at " + mPos + " at " + currentFeature.getStart());
					}
					while (mPos < seq.length()){
						newSeq[newPos] = (byte) Alphabet.DNA4().char2int(seq.charAt(mPos));
						indexes[newPos] = mPos;
						
						if (newSeq[newPos] < 3){
							newSeq[newPos] = (byte) rnd.nextInt(4);
						}
						newPos ++;
						mPos ++;
					}					
					break;//while (mPos < mSeq.length)
				}
				
				String [] toks = line.trim().split("\\t");
				//do a check
				if (!("Chr"+xx).equals(toks[0])){
					System.err.println("Wrong info 1 : " + line);
					System.exit(1);					
				}//do a check				
								
				//note STR ois 1-based index
				int pos = Integer.parseInt(toks[1]) - 1;//convert to 0-based index
				int length = Integer.parseInt(toks[2]);				
								
				//copy up to this event
				//assert nextSTR == null || nextSTR.getStart >= mPos
				while (mPos < pos){
					if (currentFeature != null && mPos >= currentFeature.getStart() - 1){
						//infact, mPos == nextFeature.getStart()			
						
						String scoreLine = sIn.readLine();
						String [] ts = scoreLine.split("\\t");
						int period = Integer.parseInt(ts[2]);
						
						double varScore = Double.parseDouble(ts[4]);
						//scale such as 
						int diff = (int) (3 * varScore * (rnd.nextDouble() + .5) * 1.5);//one way to scale
					
						if (diff <= 0){//include varScore <=0
							//no change
							diff = 0;						
						}else if (diff < currentFeature.getLength() / period && rnd.nextDouble() <= 0.5){
							//contraction
							diff = -diff;					
						}//else expantion						
						//System.out.println(seqIdx+":D = " +diff + " P = " + period + " L = " + nextFeature.getLength() + " " + varScore);
						
						
						JapsaFeature aFeature =  currentFeature.cloneFeature();
						aFeature.setStart(newPos + 1);
						aFeature.setEnd(newPos + currentFeature.getLength() + diff * period);
						aFeature.addDesc("@DIF:"+diff);
						aFeature.addDesc("@VR:"+period);						
						newAnno.add(aFeature);
						
						
						int i = currentFeature.getLength() - 1;
						int j = aFeature.getLength() - 1;
						
						for (;i >=0 && j >=0;){
							char nucleotide = Character.toUpperCase(seq.charAt(mPos + i));				
							newSeq[newPos + j] = (byte) (Alphabet.DNA4().char2int(nucleotide));
							if (newSeq[newPos + j] < 0) 
								newSeq[newPos + j] = (byte) rnd.nextInt(4);
							
							indexes[newPos + j] = mPos + i;
							
							i--;j--;
							
						}
						while(j >= 0){
							//expansion
							if (i < 0)
								i = period - 1;
							char nucleotide = Character.toUpperCase(seq.charAt(mPos + i));						
							newSeq[newPos + j] = (byte) (Alphabet.DNA4().char2int(nucleotide));
							
							indexes[newPos + j] = mPos;
							
							if (newSeq[newPos + j] < 0) 
								newSeq[newPos + j] = (byte) rnd.nextInt(4);
							
							i--;j--;
						}//while j						
						mPos += currentFeature.getLength();
						newPos += aFeature.getLength();
						
						fIndex ++;
						if (fIndex < anno.numFeatures())
							currentFeature = anno.getFeature(fIndex);
						else
							currentFeature = null;
						
						
						continue;	//while mPos < posSrc					
					}//if
					
					newSeq[newPos] = (byte)Alphabet.DNA4().char2int(seq.charAt(mPos));
					indexes[newPos] = mPos;
					if (newSeq[newPos] < 0){
						newSeq[newPos] = (byte) rnd.nextInt(4);
					}
					newPos ++;
					mPos ++;
					
				}//while mPos < posSrc
				
				if (mPos > pos) {
					System.out.println("Ignore1 " + (pos + 1) + "  " + line);
					continue;
				}
				//assert mPos == posSrc && posSrc < currentFeature.start
				
				if (currentFeature != null && mPos + (toks[3].equals("-")?0:toks[3].length()) >= currentFeature.getStart()-1){
					System.out.println("Ignore2 " + (pos + 1) + " " + line);
					continue;//ignore this because it overlaps with STR
				}
								
				if (length == 0){
					//point mutation
					char nucleotide = toks[3].charAt(0);
					if (nucleotide != seq.charAt(mPos)){
						System.err.println("Wrong info 2: " + line + " : Expect " + seq.charAt(mPos) + " see " + nucleotide);
						System.err.println(currentFeature.getStart());
						System.err.println("mPos is " + mPos);
						
						System.exit(1);
					}
					nucleotide = toks[4].charAt(0);										
					newSeq[newPos] = (byte)Alphabet.DNA4().char2int(nucleotide);
					if (newSeq[newPos] < 0)
						newSeq[newPos] = (byte) rnd.nextInt(4);
					indexes[newPos] = mPos;
					
					JapsaFeature aFeature = new JapsaFeature(newPos + 1,1);
					aFeature.setType("SNP");
					aFeature.setID("S"+(mPos+1));					
					aFeature.addDesc(line);
					newAnno.add(aFeature);
					
					mPos ++;
					newPos ++;					
				}else if (length > 0){
					//insertion
					
					JapsaFeature aFeature = new JapsaFeature(newPos + 1,toks[4].length());
					aFeature.setType("INS");
					aFeature.setID("I"+(mPos+1));
					aFeature.addDesc(line);
					newAnno.add(aFeature);
					
					if ("-".equals(toks[3])){
						//straight insert
						if (length != toks[4].length()){
							System.err.println("Wrong info 3: " + line + " : Expect length " + length + " see " + toks[4].length());
							System.exit(1);
						}
						if (length != toks[4].length()){
							System.err.println("Wrong info 4: " + line + " : Expect length " + length + " see " + toks[4].length());
							System.exit(1);
						}
						for (int t = 0; t < toks[4].length(); t++){
							char c = toks[4].charAt(t);
							newSeq[newPos+t] = (byte)Alphabet.DNA4().char2int(c);
							if (newSeq[newPos+t] < 0) 
								newSeq[newPos+t] = (byte) rnd.nextInt(4);
							
							indexes[newPos+t] = mPos;
						}
						newPos += length;
						
					}else{
						if (length + toks[3].length() != toks[4].length()){
							System.err.println("Wrong info 5: " + line + " : Expect length " + (length + toks[3].length()) + " see " + toks[4].length());
							System.exit(1);
						}
						for (int t = 0; t < toks[4].length(); t++){
							indexes[newPos+t] = mPos;
							char c = toks[4].charAt(t);
							newSeq[newPos+t] = (byte) Alphabet.DNA4().char2int(c);
							if (newSeq[newPos+t] < 0) 
								newSeq[newPos+t] = (byte) rnd.nextInt(4);														
						}
						newPos += toks[4].length();
						mPos += toks[3].length();
					}				
				}else{					
					//do deletion
					JapsaFeature aFeature = new JapsaFeature(newPos + 1,0);
					aFeature.setType("DEL");
					aFeature.setID("D"+(mPos+1));
					aFeature.addDesc(line);
					newAnno.add(aFeature);
					
					mPos+= toks[3].length();
					if (!"-".equals(toks[4])){
						//TODO: critical check below
						//was aFeature.setLength(toks[4].length());						
						aFeature.setEnd(newPos + toks[4].length());
						for (int t = 0; t < toks[4].length(); t++){
							char c = toks[4].charAt(t);
							newSeq[newPos+t] =(byte) Alphabet.DNA4().char2int(c);
							indexes[newPos+t] = mPos;
							if (newSeq[newPos+t] < 0) 
								newSeq[newPos+t] = (byte) rnd.nextInt(4);							
						}
						newPos += toks[4].length(); 
					}		
				}				
			}//while
			Sequence newSequence = new Sequence(Alphabet.DNA4(), newSeq, newPos, "S"+xx);
			newAnno.sortFeatures();
			JapsaAnnotation.write(newSequence, newAnno, out);
			System.out.println("#" + newAnno.numFeatures());
			
			for (int k = 0; k < newPos; k++){
				System.out.println((k +1) + "    " + (indexes[k] +1));			
			}
		}//for
		
		out.close();	
		sIn.close();
		reader.close();
		
		
	}
}


