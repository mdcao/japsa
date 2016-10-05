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
 * 18/10/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.tools.misc;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;

import java.io.IOException;
import java.util.Random;



/**
 * Generate a sequence with short tandem repeats
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
public class TandemRepeatSim {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		/*********************** Setting up script ****************************/		 
		String scriptName = "japsa.sim trep";
		String desc = "A program to simulate tandem repeats\n";		
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [options]");
		/**********************************************************************/		
		
		cmdLine.addDouble ("at", 0.4, "AT content, assume %A ~ %T, %C ~ %G ");
		cmdLine.addInt    ("length", 1000000, "Length of the sequence");
		cmdLine.addInt    ("unit", 3, "Repeat unit length");
		cmdLine.addDouble ("meanLen",20, "Mean number of repeat units");
		cmdLine.addDouble ("stdLen",2, "Standard deviation of repeat units");
		cmdLine.addDouble ("repeatRatio", 0.05, "Propotion of repeat DNA");
		cmdLine.addString ("output", "output", "Name of output file");
		cmdLine.addDouble ("indel", 0.02, "rate of indels");
		cmdLine.addDouble ("subs", 0.05, "rate of substitutions");		
				
		/**********************************************************************/
		cmdLine.addStdHelp();
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
		
		double at = cmdLine.getDoubleVal("at");
		double ratio = cmdLine.getDoubleVal("repeatRatio");
		
		double rIndel = cmdLine.getDoubleVal("indel");
		double rSubs = cmdLine.getDoubleVal("subs");
		if (at < 0 || at > 1){
			System.err.println("AT content has to be between 0 and 1");
			System.exit(-1);
		}
		if (ratio < 0 || ratio > 1){
			System.err.println("Repeat ration has to be between 0 and 1");
			System.exit(-1);
		}
		
		//length of sequence
		int length = cmdLine.getIntVal("length");
		
		//unit length
		int unit = cmdLine.getIntVal ("unit");
		
		//mean and std of repeats in units 
		double meanLen = cmdLine.getDoubleVal("meanLen");
		double stdLen = cmdLine.getDoubleVal("stdLen");
		
		//number of repeats
		int numReps = (int) (length * ratio / unit / meanLen);		
		
		//Generate the sequence
		Alphabet alphabet = Alphabet.DNA4();
		double [] freqs = new double[4];
		freqs[0] = freqs[3] = at / 2.0;
		freqs[1] = freqs[2] = (1.0 - at) / 2.0;
		
		Sequence seq = Sequence.random(alphabet, length, freqs);
		
		//System.out.println(numReps);	
		
		//Generate repeats		
		Random rand = new Random();
		byte [] rUnits = new byte[unit];
		
		JapsaAnnotation anno = new JapsaAnnotation(seq);
		
		//the generation of repeat start here
		for (int n = 0; n < numReps; n++){
			
			//generate the length of repeats
			int rLen = (int) (rand.nextGaussian() * stdLen + meanLen);
			//make sure rLen >= 2
			if (rLen < 2)
				rLen = 2;
			
			//convert to nucleotide
			rLen *= unit;				
			
			//generate the position of repeats, making sure no repeats overlap
			int pos = rand.nextInt(length / numReps - rLen);
			pos += (length / numReps) * n;
			
			TandemRepeat rep = new TandemRepeat (seq.getName(), pos + 1, pos + rLen);//
			String repStr = "";
			//get repeat units
			int i;
			for (i = 0; i < unit;i++){
				rUnits[i] = seq.getBase(pos + i);
				repStr = repStr + alphabet.int2char(rUnits[i]);
			}			
			
			rep.addDesc("@U:"+repStr);
			
			repStr = repStr + "\t" + repStr+' ';
			
			
			rep.setID("P"+(pos+1));
			rep.setPeriod(unit);
			rep.setUnitNo(rLen / unit);			
			
			anno.add(rep);
			String change = "C:";			
			
			i = 0;//index to unit
			
			for (int p = unit; p < rLen; p++ ){
				//toss the coin
				double chance = rand.nextDouble();
				
				if (chance < rIndel / 2){
					//insertion
					change = change + "I"+(p+1)+":";
					repStr += seq.charAt(pos + p);
										
				}else if (chance < rIndel){
					//deletion
					change = change + "D"+p+":";
					p --;
					i ++;
					if (i >= unit){
						i = 0;
						repStr += ' ';
					}
					
				}else if(chance < rIndel + rSubs){
					//substitution
					byte b = (byte) ((rUnits[i] + 1 + rand.nextInt(alphabet.size() - 1)) % alphabet.size());
					seq.setBase(pos + p, b);
					change = change + "S"+(p+1)+":";
					repStr += seq.charAt(pos + p);
					i ++;
					if (i >= unit){
						i = 0;
						repStr += ' ';
					}							
				}else{
					//normal
					seq.setBase(pos + p, rUnits[i]);
					repStr += seq.charAt(pos + p);
					
					//advance i
					i ++;
					if (i >= unit){
						i = 0;
						repStr += ' ';
					}
				}				
			}//for		

			rep.addDesc(repStr);
			rep.addDesc(change+"L"+rLen);
		}		
		
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("output")+".jsa"); 
		JapsaAnnotation.write(seq, anno, out) ;
		out.close();
		
		out = SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("output")+".fas"); 
		seq.writeFasta(out);
		out.close();
	}

}
