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

package japsa.bio.bac;



import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

/**
 * @author minhduc
 *
 */
public class MLSTyping{
	
	int numGenes = 7;	
	ArrayList<MLSTProfile> profiles;
	String [] genes;
		
	/**
	 * A redesign
	 * @param table
	 * @throws IOException
	 */
	public MLSTyping(String table) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(table));
		String line = br.readLine();
		String [] toks = line.trim().split("\t");		
		
		genes = new String[numGenes];
		profiles = new ArrayList<MLSTProfile>();

		for (int i =0; i < genes.length;i ++)
			genes[i] = toks[i+1];
		
		while ((line = br.readLine()) != null) {			
			toks = line.trim().split("\t");
			MLSTProfile profile = new MLSTProfile(toks[0], numGenes);
			
			for (int i = 0; i < genes.length; i++){
				profile.indivisialScore[i]  = Integer.parseInt(toks[i+1]);
			}
			
			profiles.add(profile);
		}
		br.close();
	}

	public static String bestMlst(ArrayList<Sequence> seqs, String allelesDir, String table, String blastn) throws InterruptedException, IOException{
		BufferedReader br = new BufferedReader(new FileReader(table));
		String line = br.readLine();
		String [] toks = line.trim().split("\t");		
		String [] genes = new String[7];		 

		for (int i =0; i < genes.length;i ++)
			genes[i] = toks[i+1];		

		HashMap<String, String> stMap = new HashMap<String, String>();

		while ((line = br.readLine()) != null) {
			toks = line.trim().split("\t");

			String ST = "ST_" + toks[0];
			String stKey = "";
			for (int i = 0; i < genes.length; i++){
				stKey += genes[i] + "_" + toks[i+1] + "_";
			}
			stMap.put(stKey, ST);			
		}
		br.close();

		//Read in alleles
		String key = "";
		int typeScore = 0;

		for (String gene:genes){			
			ProcessBuilder pb = new ProcessBuilder(blastn, "-subject", "-",
				"-query", allelesDir + "/" + gene + ".fas", "-outfmt", "6 qseqid qlen nident gaps mismatch");

			String allele = "";
			int score = 1000000;//less is more

			Process process = pb.start();

			SequenceOutputStream out = new SequenceOutputStream(process.getOutputStream());
			for (Sequence seq:seqs){
				seq.writeFasta(out);
			}
			out.close();

			br = new BufferedReader(new InputStreamReader(process.getInputStream()));			

			while ((line = br.readLine()) != null) {
				toks = line.trim().split("\t");
				int myScore = Integer.parseInt(toks[1]) - Integer.parseInt(toks[2]) + 2 *  Integer.parseInt(toks[3]);
				if (myScore < score){
					allele = toks[0];
					score = myScore;					
				}
			}
			br.close();

			process.waitFor();
			key += allele + "_";
			typeScore += score;
		}

		String st = stMap.get(key);
		if (st == null)
			st = "NEW";

		return (key + "\t" + st + "\t" + typeScore);
	}

	public static class MLSTProfile implements Comparable<MLSTProfile>{		
		int typeScore;		
		String type;
		int [] indivisialScore;

		MLSTProfile(String st, int numGenes){
			type = st;
			typeScore = 0;
			indivisialScore = new int[numGenes];
			for (int i =0; i < indivisialScore.length; i++){
				indivisialScore[i] = 1000000;	
				typeScore += indivisialScore[i]; 
			}			
		}
		@Override
		public int compareTo(MLSTProfile o) {
			return Integer.compare(typeScore, o.typeScore);
		}
	}
}
