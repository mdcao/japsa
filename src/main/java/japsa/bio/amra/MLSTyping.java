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
 * 28/05/2014 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsa.bio.amra;



import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

/**
 * @author minhduc
 *
 */
public class MLSTyping{

	int numGenes = 7;	
	public ArrayList<MLSType> profiles;
	private String [] geneNames;
	private ArrayList<Sequence> [] geneSeqs;
	private int alleleNo2AlleleIndex[][];
	//Link the profile Number to the geneIndex. This is to make indexing transparant
	//e.g., inernalLink[3][4] point to the index in geneSeqs[3] of allele 4

	String baseDir = ".";
	BitSet[] bitSets;
	/**
	 * A redesign
	 * @param mlstBase: path to mlst scheme, containing profile.dat file
	 * @throws IOException
	 */

	@SuppressWarnings("unchecked")
	public MLSTyping(String mlstBase) throws IOException{
		baseDir = mlstBase;
		BufferedReader br = new BufferedReader(new FileReader(mlstBase + "/profile.dat"));
		String line = br.readLine();
		String [] toks = line.trim().split("\t");

		bitSets = new BitSet[numGenes];
		geneNames = new String[numGenes];
		profiles = new ArrayList<MLSType>();

		for (int i =0; i < numGenes;i ++){
			geneNames[i] = toks[i+1];
			bitSets[i] = new BitSet();
		}

		while ((line = br.readLine()) != null) {			
			toks = line.trim().split("\t");
			MLSType profile = new MLSType("ST" + toks[0], numGenes);

			for (int i = 0; i < geneNames.length; i++){
				int alleleNo = Integer.parseInt(toks[i+1]);				
				profile.geneAlleles[i] = alleleNo;
				bitSets[i].set(alleleNo);
			}
			profiles.add(profile);
		}
		br.close();

		geneSeqs = new ArrayList [numGenes];
		alleleNo2AlleleIndex = new int[numGenes][];		

		for (int i = 0; i < numGenes; i++){			
			geneSeqs[i] = SequenceReader.readAll(mlstBase + "/" + geneNames[i] + ".fas", Alphabet.DNA());

			int x = geneSeqs[i].size() - 1;
			String alleleName = geneSeqs[i].get(x).getName();
			
			int index_ = alleleName.lastIndexOf('_');
			
			int alleleNo = Integer.parseInt(alleleName.substring(1 + index_));
			alleleNo2AlleleIndex[i] = new int[alleleNo + 1];
			Arrays.fill(alleleNo2AlleleIndex[i], -1);

			for (; x >=0;x--){
				alleleName = geneSeqs[i].get(x).getName();
				index_ = alleleName.lastIndexOf('_');				
				alleleNo = Integer.parseInt(alleleName.substring(1 + index_));
				alleleNo2AlleleIndex[i][alleleNo] = x;
			}//for
		}
	}

	/**
	 * Return list of alleles for the gene with number geneNo
	 * @param geneNo
	 * @return
	 */
	public ArrayList<Sequence> alleles(int geneNo){
		return geneSeqs[geneNo];
	}

	/**
	 * Determine if a particular is used (path of some ST)
	 * @param geneNo
	 * @param alleleNo
	 * @return
	 */
	public boolean useAlleleNo(int geneNo, int alleleNo){
		return bitSets[geneNo].get(alleleNo);
	}

	
	public int alleleNo2AlleleIndex(int geneNo, int alleleNo){
		return alleleNo2AlleleIndex[geneNo][alleleNo];
	}
	public String getBase(){
		return baseDir;
	}

	public ArrayList<MLSType> getProfiles(){
		return profiles;
	}
	
	public String getGeneName(int geneNo){
		return geneNames[geneNo];
	}
	
	public int alleleIndex(MLSType type, int geneNo){
		return alleleNo2AlleleIndex[geneNo][type.getAllele(geneNo)];
	}

	/**
	 * Return the best ST
	 * @param seqs
	 * @param allelesDir
	 * @param table
	 * @param blastn
	 * @return
	 * @throws InterruptedException
	 * @throws IOException
	 */
	public static String bestMlst(ArrayList<Sequence> seqs, String mlstDir) throws InterruptedException, IOException{
		MLSTyping mlstScheme = new MLSTyping(mlstDir);
		String line;		
		HashMap<String, String> stMap = new HashMap<String, String>();
		for (MLSType profile:mlstScheme.profiles){
			String ST = profile.getST();

			String stKey = "";
			for (int i = 0; i < mlstScheme.numGenes; i++){
				stKey += mlstScheme.geneNames[i] + "_" + profile.geneAlleles[i] + "|";
			}
			stMap.put(stKey, ST);
		}		

		//Read in alleles
		String key = "";
		int typeScore = 0;

		for (String gene:mlstScheme.geneNames){			
			ProcessBuilder pb = new ProcessBuilder("blastn", "-subject", "-",
				"-query", mlstDir + "/" + gene + ".fas", "-outfmt", "6 qseqid qlen nident gaps mismatch");

			String allele = "";
			int score = 1000000;//less is more

			Process process = pb.start();

			SequenceOutputStream out = new SequenceOutputStream(process.getOutputStream());
			for (Sequence seq:seqs){
				seq.writeFasta(out);
			}
			out.close();

			BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));			

			while ((line = br.readLine()) != null) {
				String [] toks = line.trim().split("\t");
				int myScore = Integer.parseInt(toks[1]) - Integer.parseInt(toks[2]) + 2 *  Integer.parseInt(toks[3]);
				if (myScore < score){
					allele = toks[0];
					score = myScore;					
				}
			}
			br.close();

			process.waitFor();
			key += allele + "|";
			typeScore += score;
		}

		String st = stMap.get(key);
		if (st == null)
			st = "NEW";

		return (key + "\t" + st + "\t" + typeScore);
	}


	public static MLSTyping topMlst(ArrayList<Sequence> seqs, String mlstDir) throws InterruptedException, IOException{
		MLSTyping mlstScheme = new MLSTyping(mlstDir);
		int [][] scoreMatrix = new int[mlstScheme.numGenes][];
		
		
		//can run in parallele with i
		for (int i = 0; i < mlstScheme.numGenes; i++){
			String gene = mlstScheme.geneNames[i];
			scoreMatrix[i] = new int[mlstScheme.geneSeqs[i].size()];
			Arrays.fill(scoreMatrix[i], 100000);			
			ProcessBuilder pb = new ProcessBuilder("blastn", "-subject", "-",
				"-query", mlstDir + "/" + gene + ".fas", "-outfmt", "6 qseqid qlen nident gaps mismatch");

			Process process = pb.start();
			SequenceOutputStream out = new SequenceOutputStream(process.getOutputStream());
			for (Sequence seq:seqs){
				seq.writeFasta(out);
			}
			out.close();

			BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));			

			String line;
			while ((line = br.readLine()) != null) {
				String [] toks = line.trim().split("\t");
				//This is the magic line to compute score, need to make it modular later
				int myScore = Integer.parseInt(toks[1]) - Integer.parseInt(toks[2]) + 2 *  Integer.parseInt(toks[3]);
				int alleleIndex = mlstScheme.alleleNo2AlleleIndex[i][Integer.parseInt(toks[0].split("_")[1])];
				if (scoreMatrix[i][alleleIndex] > myScore)
					scoreMatrix[i][alleleIndex] = myScore;
			}
			br.close();
			process.waitFor();
		}

		for (MLSType type:mlstScheme.profiles){
			type.typeScore = 0;
			for (int i = 0; i < mlstScheme.numGenes;i++){
				type.typeScore += scoreMatrix[i][mlstScheme.alleleNo2AlleleIndex[i][type.geneAlleles[i]]];
			}			
		}
		Collections.sort(mlstScheme.profiles);
		return mlstScheme;
	}


	public static class MLSType implements Comparable<MLSType>{		
		public double typeScore;		
		private String st;		
		private int [] geneAlleles;

		MLSType(String st, int numGenes){
			this.st = st;
			typeScore = 100000;			
			geneAlleles = new int[numGenes];
		}

		@Override
		public int compareTo(MLSType o) {
			return Double.compare(typeScore, o.typeScore);
		}

		public String getST(){
			return st;
		}
		
		public double getScore(){
			return typeScore;
		}
		
		public int getAllele(int geneNo){
			return geneAlleles[geneNo];
		}		
	}
}
