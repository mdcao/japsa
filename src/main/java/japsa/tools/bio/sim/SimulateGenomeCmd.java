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
 * 10/06/2016: Revisit
 ****************************************************************************/
package japsa.tools.bio.sim;

import japsa.seq.AbstractSequence;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

/**
 * Simulate a haploid genome from a reference
 * 
 * @author minhduc
 * 
 */
@Deployable(
		scriptName = "jsa.sim.genome", 
		scriptDesc = "Simulate genomes with variation from an existing genome"
		)
public class SimulateGenomeCmd extends CommandLine{

	public SimulateGenomeCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("input", null, "Name of input file", true);
		addString("output", null, "Prefix output file", true);
		addString("logFile", "-", "Log file, - for stadard error");
		addString("sv", "", "List of deletions/insertions eg chr1:1000-1020:-,chr2:1300-1340:+");		
		addDouble("snp", 0.0025, "SNP rate");
		addDouble("indel", 0.00028, "indel rate");
		//Defaults parameters from doi:10.1038/nature09534:
		//5.9M SNPs, 650K indels (1-50) on 2.3Gb  
		addDouble("ext", 0.5, "indel extension rate");
		addInt("seed", 0, "Random seed, 0 for a random seed");

		addStdHelp();		

	}
	/**
	 * @param args
	 */

	public static void main(String[] args) throws IOException {		 		
		CommandLine cmdLine = new SimulateGenomeCmd();		
		args = cmdLine.stdParseLine(args);

		String inFile = cmdLine.getStringVal("input");
		String outFile = cmdLine.getStringVal("output");
		String logFile = cmdLine.getStringVal("logFile");

		String svOption = cmdLine.getStringVal("sv");

		double snp = cmdLine.getDoubleVal("snp");
		double indel = cmdLine.getDoubleVal("indel");
		double ext = cmdLine.getDoubleVal("ext");
		int seed = cmdLine.getIntVal("seed");
		


		ArrayList<StructualVarition> svs = new ArrayList<StructualVarition>();
		String [] toks = svOption.trim().split(",");

		for (int i =0; i < toks.length;i++){
			svs.add(StructualVarition.parseSV(toks[i]));
		}

		int svsIndex = 0;

		
		//generate a random seed if need to
		seed = seed(seed);
		Random rnd = new Random(seed);
		
		SequenceReader reader = SequenceReader.getReader(inFile);

		SequenceOutputStream logOS = 	logFile.equals("-")? 
				(new SequenceOutputStream(System.err)) 
				: 
				(SequenceOutputStream.makeOutputStream(logFile));
				

		SequenceOutputStream outFasta = SequenceOutputStream
				.makeOutputStream(outFile);

		//SequenceOutputStream outJsa = SequenceOutputStream
		//		.makeOutputStream(outFile + ".jsa");

		AbstractSequence seq = null;

		logOS.print("#Seed " + seed + "\n");

		StructualVarition sv = null;
		if (svsIndex < svs.size()){
			sv = svs.get(svsIndex);
		}

		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			int length = seq.length();
			//Step 0: Introduce any structural variation first			
			int currentIndex = 0;			
			SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA(), length  + length/5, seq.getName());
			
			logOS.print("#Start with " + seq.getName() + "  " + seq.length() +"\n");
			
			while (sv !=null && sv.chr.equals(seq.getName())){
				int start = sv.start;
				int end = sv.end;
				for (;currentIndex < start - 1 && currentIndex < seq.length();currentIndex++){
					sb.append(seq.getBase(currentIndex));
				}
				
				logOS.print("#  " + sv.chr + ":" + sv.start + "-" + sv.end + ":" + ((sv.svType == StructualVarition.DELETION)?"-":"+") + "\n");
				
				if (sv.svType == StructualVarition.DELETION){
					for (int x = currentIndex;x < end-1;x++){
						sb.append(seq.getBase(x));
					}
					//one more time
					for (int x = currentIndex;x < end-1;x++){
						sb.append(seq.getBase(x));
					}
					
				}else if (sv.svType == StructualVarition.DUPLICATION){
					//doing nothing
				}else{
					Logging.error("Dont know what to do");
				}				
				currentIndex = end - 1;
				
				svsIndex ++;
				if (svsIndex < svs.size()){
					sv = svs.get(svsIndex);
				}else
					sv = null;
				
			}
			for (;currentIndex < seq.length();currentIndex++){
				sb.append(seq.getBase(currentIndex));
			}			
			
			seq = sb;
			//step 1: introduce white noise						
			logOS.print("#Restart with " + seq.getName() + "  " + seq.length() +"\n");
			byte[] seqByte = new byte[length + length / 5];			
			int currentInx = 0;
			int numSNPs = 0;
			int numIndels = 0;

			int index = 0;
			for (; index < seq.length();) {
				byte base = seq.getBase(index);
				if (base >=4){
					seqByte[currentInx] = (byte) (rnd.nextInt(4));
					currentInx++;
				}else{
					double val = rnd.nextDouble();
					if (val <= snp) {// SNPs
						// An SNP
						// Generate a random number between 0-2, then plus 1 and
						// plus the index of the previous char
						// to avoid generating the same nucleotide						
						seqByte[currentInx] = (byte) ((1 + rnd.nextInt(3) + base) % 4);
						currentInx++;
						numSNPs++;
						logOS.print("SNP " + (index + 1) + " at " + currentInx + "\n");
					} else if (val <= snp + indel) {// indel
						numIndels++;
						val = rnd.nextDouble();
						if (val >= 0.5) {
							int size = 1;
							while (rnd.nextDouble() < ext) {
								size++;
								index++;
							}
							logOS.print("DEL " + (index + 1 - size) + " at " + currentInx + " of " + size + "\n");
							numIndels++;
							// A deletion
							// currentInx[seqIdx] --;
						} else {// an insertion
							seqByte[currentInx] = base;
							currentInx++;

							int size = 0;
							do {
								seqByte[currentInx] = (byte) (rnd.nextInt(4));
								currentInx++;
								size++;
							} while (rnd.nextDouble() < ext);
							logOS.print("INS " + (index + 1) + " at " + currentInx + " of " + size + "\n");
							numIndels++;
						}// insert vs delete
					} else {// Direct copy
						seqByte[currentInx] = base;
						currentInx++;
					}
				}// if
				index++;
			}// for index			
			// System.out.println("Insert reps done");
			Sequence nSeq = new Sequence(Alphabet.DNA4(), seqByte, currentInx,
					seq.getName());
			nSeq.writeFasta(outFasta);			

			logOS.print("Sequence " + seq.getName() + " : " + numSNPs + " SNPs "
					+ numIndels + " indels \n");

			//JapsaAnnotation.write(nSeq, null, outJsa);
		}// for
		outFasta.close();
		logOS.close();		
		//outJsa.close();
	}
	
	/**
	 * Generate a random seed if the input <=0
	 * @param seed
	 * @return
	 */
	public static int seed(int seed){
		if (seed <= 0)
			seed = new Random().nextInt();
		
		//make sure seed is not negative
		if (seed <0 )
			seed = - seed;
		
		return seed;
	}

	static class StructualVarition{
		static final int DELETION = 0;		
		static final int DUPLICATION = 1;

		String chr;		
		int svType;
		int start;
		int end;

		static StructualVarition parseSV(String str){
			StructualVarition sv = new StructualVarition(); 
			//example:	chr1:1000-1020:-,chr2:1300-1340:+
			String [] toks = str.split(":");
			sv.chr = toks[0];

			if (toks.length<2)
				sv.svType = StructualVarition.DELETION;
			else if (toks[2].charAt(0) == '-')
				sv.svType = StructualVarition.DELETION;
			else  if (toks[2].charAt(0) == '+')
				sv.svType = StructualVarition.DUPLICATION;
			else 
				return null;

			toks = toks[1].split("-");

			sv.start = (int) Double.parseDouble(toks[0]);
			sv.end = (int) Double.parseDouble(toks[1]);			

			return sv;
		}

	}

}
