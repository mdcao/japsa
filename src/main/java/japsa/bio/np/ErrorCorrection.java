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
 * 3 Apr 2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.bio.np;

import japsa.seq.Alphabet;
import japsa.seq.Alphabet.DNA;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Provides functionality for error correction and find the best allele/sequence.
 * Potentially will use HMM/FSM 
 * 
 * @author minhduc
 * 
 *
 */
public class ErrorCorrection {
    private static final Logger LOG = LoggerFactory.getLogger(ErrorCorrection.class);

	public static String prefix = "tmp";
	public static String msa = "kalign";

	public static double needle(Sequence seq1, Sequence seq2, String prefix) throws IOException, InterruptedException{
		String seq1File = prefix + seq1.getName() + ".fasta";
		String seq2File = prefix + seq2.getName() + ".fasta";
		seq1.writeFasta(seq1File);
		seq2.writeFasta(seq2File);
		String needleOut = prefix + "_alignment.needle";
		String cmd = "needle -gapopen 10 -gapextend 0.5 -asequence " 
			+ seq1File + " -bsequence " + seq2File + " -outfile " + needleOut;
		LOG.info("Running " + cmd);
		Process process = Runtime.getRuntime().exec(cmd);
		process.waitFor();					
		LOG.info("Run'ed " + cmd );

		BufferedReader scoreBf = new BufferedReader(new FileReader(needleOut));
		String scoreLine = null;					
		double score = 0;
		while ((scoreLine = scoreBf.readLine())!=null){
			String [] scoreToks = scoreLine.split(" ");					
			if (scoreToks.length == 3 && scoreToks[1].equals("Score:")){
				score = Double.parseDouble(scoreToks[2]);
				break;//while
			}					
		}//while
		scoreBf.close();
		return score;		
	}

    public static Sequence consensusSequence(ArrayList<Sequence> readList, String prefix, String msa) throws IOException, InterruptedException{
    	if (readList != null && readList.size() > 0)
    		return consensusSequence(readList, readList.size(), prefix, msa);
    	else
    		return null;
    	
    }
	public static Sequence consensusSequence(ArrayList<Sequence> readList, int max, String prefix, String msa) throws IOException, InterruptedException{
		//String faiFile = prefix + "_" + this.currentReadCount;
		Sequence consensus = null;
		if (readList != null && readList.size() > 0){
			//1.0 write alignment to faiFile
			if (readList.size() == 1){
				readList.get(0).setName("consensus");
				consensus = readList.get(0);
			}else{
				String faiFile = prefix + "_ai.fasta";//name of fasta files of reads mapped to the gene				
				String faoFile = prefix + "_ao.fasta";//name of fasta files of reads mapped to the gene
				{
					SequenceOutputStream faiSt = SequenceOutputStream.makeOutputStream(faiFile);
					int count = 0;
					for (Sequence seq:readList){
						LOG.info(seq.getName() + "  " + seq.length());
						seq.writeFasta(faiSt);
						count ++;
						if (count >= max)
						    break;//for
					}
					faiSt.close();
				}

				//2.0 Run multiple alignment
				{
					String cmd  = "";
					if (msa.startsWith("poa")){						
						cmd = "poa -read_fasta " + faiFile + " -clustal " + faoFile + " -hb blosum80.mat";
					}else if (msa.startsWith("muscle")){
						cmd = "muscle -in " + faiFile + " -out " + faoFile + " -maxiters 5 -quiet";				
					}else if (msa.startsWith("clustal")) {
						cmd = "clustalo --force -i " + faiFile + " -o " + faoFile;
					}else if (msa.startsWith("kalign")){
						cmd = "kalign -gpo 60 -gpe 10 -tgpe 0 -bonus 0 -q -i " + faiFile	+ " -o " + faoFile;
					}else if (msa.startsWith("msaprobs")){
						cmd = "msaprobs -o " + faoFile + " " + faiFile;
					}else if (msa.startsWith("mafft")){
						cmd = "mafft_wrapper.sh  " + faiFile + " " + faoFile;
					}else{
						LOG.error("Unknown msa function " + msa);
						return null;
					}

					LOG.info("Running " + cmd);
					Process process = Runtime.getRuntime().exec(cmd);
					process.waitFor();
					LOG.info("Done " + cmd);
				}


				if ("poa".equals(msa)){
					SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA(), readList.get(0).length());
					BufferedReader bf =  FastaReader.openFile(faoFile);
					String line = bf.readLine();
					while ( (line = bf.readLine()) != null){
						if (line.startsWith("CONSENS0")){
							for (int i = 10;i < line.length();i++){
								char c = line.charAt(i);
								int base = DNA.DNA().char2int(c);
								if (base >= 0 && base < 4)
									sb.append((byte)base);
							}//for							
						}//if
					}//while
					sb.setName("consensus");
					LOG.info(sb.getName() + "  " + sb.length());
					return sb.toSequence();
				}

				//3.0 Read in multiple alignment
				ArrayList<Sequence> seqList = new ArrayList<Sequence>();
				{
					SequenceReader msaReader = FastaReader.getReader(faoFile);
					Sequence nSeq = null;
					while ((nSeq = msaReader.nextSequence(Alphabet.DNA())) != null) {
						seqList.add(nSeq);						
					}
					msaReader.close();
				}

				//4.0 get consensus and write to facFile				
				{
					int [] coef = new int[seqList.size()];			
					for (int y = 0; y < seqList.size(); y++){
						coef[y] = 1;
						//if (seqList.get(y).getName().contains("twodi"))
						//	coef[y] = 2;				
					}

					//TODO: combine error profiles?? 
					int [] counts = new int[6];

					SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA(), seqList.get(0).length());
					for (int x = 0; x < seqList.get(0).length();x++){		
						Arrays.fill(counts, 0);
						for (int y = 0; y < seqList.size(); y++){					
							byte base = seqList.get(y).getBase(x);
							if (base >= 6) 
								counts[4] += coef[y];//N
							else
								counts[base] += coef[y];
						}//for y
						int maxIdx = 0;
						for (int y = 1; y < counts.length; y++){
							if (counts[y] > counts[maxIdx])
								maxIdx = y;					
						}//for y
						if (maxIdx < Alphabet.DNA.GAP){//not a gap
							sb.append((byte)maxIdx);
						}//if
					}//for x
					sb.setName("consensus");
					LOG.info(sb.getName() + "  " + sb.length());
					consensus = sb.toSequence();
				}
			}			

		}
		return consensus;
	}
}
