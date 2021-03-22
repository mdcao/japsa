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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.seq.Alphabet;
import japsa.seq.Alphabet.DNA;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

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

    public static Sequence consensusSequence(List<Sequence> readList, String prefix, String msa) throws IOException, InterruptedException{
    	if (readList != null && readList.size() > 0)
    		return consensusSequence(readList, readList.size(), prefix, msa);
    	else
    		return null;
    	
    }
    
    public static void writeAlignmentToFaiFile(List<Sequence> readList,String faiFile, String faoFile , int max)throws IOException, InterruptedException{
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
    }
    
    public static void runMultipleAlignment(String faiFile, String faoFile) throws IOException, InterruptedException{
    	{
			String[] cmd;
    	    File temp = File.createTempFile("tempfile", ".tmp"); 
			if (msa.startsWith("poa")){						
				//create a temporary matrix file similar to blosum80.mat	    		
	    	    PrintWriter printer = new PrintWriter(new BufferedWriter(new FileWriter(temp)));
	    	    printer.println("GAP-PENALTIES=12 6 6");
	    	    printer.println("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  ?  a  g  t  c  u  ]  n");
	    	    printer.println("A  7 -3 -3 -3 -1 -2 -2  0 -3 -3 -3 -1 -2 -4 -1  2  0 -5 -4 -1 -3 -2 -1 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"R -3  9 -1 -3 -6  1 -1 -4  0 -5 -4  3 -3 -5 -3 -2 -2 -5 -4 -4 -2  0 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"N -3 -1  9  2 -5  0 -1 -1  1 -6 -6  0 -4 -6 -4  1  0 -7 -4 -5  5 -1 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"D -3 -3  2 10 -7 -1  2 -3 -2 -7 -7 -2 -6 -6 -3 -1 -2 -8 -6 -6  6  1 -3 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"C -1 -6 -5 -7 13 -5 -7 -6 -7 -2 -3 -6 -3 -4 -6 -2 -2 -5 -5 -2 -6 -7 -4 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"Q -2  1  0 -1 -5  9  3 -4  1 -5 -4  2 -1 -5 -3 -1 -1 -4 -3 -4 -1  5 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"E -2 -1 -1  2 -7  3  8 -4  0 -6 -6  1 -4 -6 -2 -1 -2 -6 -5 -4  1  6 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"G  0 -4 -1 -3 -6 -4 -4  9 -4 -7 -7 -3 -5 -6 -5 -1 -3 -6 -6 -6 -2 -4 -3 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"H -3  0  1 -2 -7  1  0 -4 12 -6 -5 -1 -4 -2 -4 -2 -3 -4  3 -5 -1  0 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"I -3 -5 -6 -7 -2 -5 -6 -7 -6  7  2 -5  2 -1 -5 -4 -2 -5 -3  4 -6 -6 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"L -3 -4 -6 -7 -3 -4 -6 -7 -5  2  6 -4  3  0 -5 -4 -3 -4 -2  1 -7 -5 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"K -1  3  0 -2 -6  2  1 -3 -1 -5 -4  8 -3 -5 -2 -1 -1 -6 -4 -4 -1  1 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"M -2 -3 -4 -6 -3 -1 -4 -5 -4  2  3 -3  9  0 -4 -3 -1 -3 -3  1 -5 -3 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"F -4 -5 -6 -6 -4 -5 -6 -6 -2 -1  0 -5  0 10 -6 -4 -4  0  4 -2 -6 -6 -3 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"P -1 -3 -4 -3 -6 -3 -2 -5 -4 -5 -5 -2 -4 -6 12 -2 -3 -7 -6 -4 -4 -2 -3 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"S  2 -2  1 -1 -2 -1 -1 -1 -2 -4 -4 -1 -3 -4 -2  7  2 -6 -3 -3  0 -1 -1 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"T  0 -2  0 -2 -2 -1 -2 -3 -3 -2 -3 -1 -1 -4 -3  2  8 -5 -3  0 -1 -2 -1 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"W -5 -5 -7 -8 -5 -4 -6 -6 -4 -5 -4 -6 -3  0 -7 -6 -5 16  3 -5 -8 -5 -5 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"Y -4 -4 -4 -6 -5 -3 -5 -6  3 -3 -2 -4 -3  4 -6 -3 -3  3 11 -3 -5 -4 -3 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"V -1 -4 -5 -6 -2 -4 -4 -6 -5  4  1 -4  1 -2 -4 -3  0 -5 -3  7 -6 -4 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"B -3 -2  5  6 -6 -1  1 -2 -1 -6 -7 -1 -5 -6 -4  0 -1 -8 -5 -6  6  0 -3 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"Z -2  0 -1  1 -7  5  6 -4  0 -6 -5  1 -3 -6 -2 -1 -2 -5 -4 -4  0  6 -1 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"X -1 -2 -2 -3 -4 -2 -2 -3 -2 -2 -2 -2 -2 -3 -3 -1 -1 -5 -3 -2 -3 -1 -2 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"? -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"a -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9  4 -2 -2 -2 -2 -9  0\n"
	    	    			+	"g -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -2  4 -2 -2 -2 -9  0\n"
	    	    			+	"t -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -2 -2  4 -2  4 -9  0\n"
	    	    			+	"c -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -2 -2 -2  4 -2 -9  0\n"
	    	    			+	"u -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -2 -2  4 -2  4 -9  0\n"
	    	    			+	"] -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9\n"
	    	    			+	"n -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9  0  0  0  0  0 -9  0");
	    	    
	    	    printer.close();			    	
				//invoke
				cmd = new String[]{"poa", "-read_fasta", faiFile, "-clustal", faoFile, "-hb", temp.getAbsolutePath()};
				
			}else if (msa.startsWith("spoa")){
				cmd = new String[]{"spoa", faiFile};	
			}else if (msa.startsWith("abpoa")){
				cmd = new String[]{"abpoa", faiFile, "-o", faoFile};
			}else if (msa.startsWith("muscle")){
				cmd = new String[]{"muscle", "-in", faiFile, "-out", faoFile, "-maxiters", "5", "-quiet"};				
			}else if (msa.startsWith("clustal")) {
				cmd = new String[]{"clustalo", "--force", "-i", faiFile, "-o", faoFile};
			}else if (msa.startsWith("kalign3")){
				cmd = new String[]{"kalign", "-i", faiFile, "-o", faoFile};
			}else if (msa.startsWith("kalign")){
				cmd = new String[]{"kalign", "-gpo", "60", "-gpe", "10", "-tgpe", "0", "-bonus", "0", "-q", "-i", faiFile, "-o", faoFile};
			}else if (msa.startsWith("msaprobs")){
				cmd = new String[]{"msaprobs", "-o", faoFile, faiFile};
			}else if (msa.startsWith("mafft")){
				cmd = new String[]{"mafft_wrapper.sh", faiFile, faoFile};
			}else{
				LOG.error("Unknown msa function " + msa);
				throw new InterruptedException("Unknown msa function " + msa);
			}

			LOG.info("Running " + Arrays.toString(cmd));
			ProcessBuilder builder = new ProcessBuilder(cmd).redirectError(new File("kalign_err.txt"));
			if (msa.startsWith("spoa")){
				builder.redirectOutput(new File(faoFile));			
			}
			Process process = builder.start();
			process.waitFor();
			LOG.info("Done " + Arrays.toString(cmd));
			temp.deleteOnExit();
		}
    }
    
    public static Sequence 	readPOAOutput(String faoFile, int seql) throws IOException{
    	SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA(), seql);
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
    public static Sequence 	readSPOAOutput(String faoFile, int seql) throws IOException{
    	SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA(), seql);
		BufferedReader bf =  FastaReader.openFile(faoFile);
		String line = bf.readLine();
		while ( (line = bf.readLine()) != null){
			if (line.startsWith("Consensus")){
				continue;						
			}//if
			sb.append(new Sequence(Alphabet.DNA(), line, "consensus"));
		}//while
		sb.setName("consensus");
		LOG.info(sb.getName() + "  " + sb.length());
		return sb.toSequence();
    }
    public static ArrayList<Sequence> readMultipleAlignment(String faoFile) throws IOException{
    	ArrayList<Sequence> seqList = new ArrayList<Sequence>();
		{
			SequenceReader msaReader = FastaReader.getReader(faoFile);
			Sequence nSeq = null;
			while ((nSeq = msaReader.nextSequence(Alphabet.DNA())) != null) {
				seqList.add(nSeq);						
			}
			msaReader.close();
		}
		return seqList;

    }

    public static Sequence getConsensus(ArrayList<Sequence> seqList){
    	
    		Sequence consensus = null;
			int [] coef = new int[seqList.size()];			
			for (int y = 0; y < seqList.size(); y++){
				coef[y] = 1;
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
		return consensus;
    }
    
	public static Sequence consensusSequence(List<Sequence> readList, int max, String prefix, String msa) throws IOException, InterruptedException{
		//String faiFile = prefix + "_" + this.currentReadCount;
		Sequence consensus = null;
		if (readList != null && readList.size() > 0){
			if (readList.size() == 1){
				readList.get(0).setName("consensus");
				consensus = readList.get(0);
			}else{
				
				String faiFile = prefix + "_ai.fasta";//name of fasta files of reads mapped to the gene				
				String faoFile = prefix + "_ao.fasta";//name of fasta files of reads mapped to the gene
				//1.0 write alignment to faiFile

				writeAlignmentToFaiFile(readList, faiFile, faoFile, max);
				
				//2.0 Run multiple alignment
				try{
					runMultipleAlignment(faiFile, faoFile);
				}catch(InterruptedException exc){
					return null;
				}
				
				if ("poa".equals(msa))
					return readPOAOutput(faoFile, readList.get(0).length());	
				else if("spoa".equals(msa)){
					return readSPOAOutput(faoFile, readList.get(0).length());
				}else if("abpoa".equals(msa)) {
					SequenceReader reader = SequenceReader.getReader(faoFile);
					return reader.nextSequence(DNA.DNA());
					
				}
				

				//3.0 Read in multiple alignment
				ArrayList<Sequence> seqList = readMultipleAlignment(faoFile);
				//4.0 get consensus and write to facFile	
				
				consensus = getConsensus(seqList);
				
				Files.deleteIfExists(Paths.get(faiFile));
				Files.deleteIfExists(Paths.get(faoFile));
			}
		}
		return consensus;
	}
}
