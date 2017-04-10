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

/*****************************************************************************
 *                           Revision History                                
 * 7 Sep 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsa.tools.amra;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;

import japsa.bio.amra.ResistanceGeneDB;
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
	scriptName = "jsa.amra.genome2res",
	scriptDesc = "Finding resistance genes/classes in a genome"
	)
public class Genomes2ResistanceGeneCmd extends CommandLine {

	public Genomes2ResistanceGeneCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("input", null, "Name of the genome file",true);
		addString("output", null, "Name of the output file",true);		
		addString("resDB", null, "Name of the resistance gene database",true);

		addDouble("identity", 0.85, "Minimum identity");
		addDouble("coverage", 0.85, "Minimum coverage of gene");

		addStdHelp();
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLine cmdLine = new Genomes2ResistanceGeneCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String input =  cmdLine.getStringVal("input");
		String output = cmdLine.getStringVal("output");
		String resDBPath = cmdLine.getStringVal("resDB");

		double identity = cmdLine.getDoubleVal("identity");
		double coverage = cmdLine.getDoubleVal("coverage");
		
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output); 

		//Note: The database would contain files 
		//  resDB.fas -- gene sequences 
		//  resDB.map -- map genes (corresponds to resDB.fas) to resistance classes
		ResistanceGeneDB resDB = new ResistanceGeneDB(resDBPath);		
		ArrayList<Sequence> seqs = SequenceReader.readAll(input, Alphabet.DNA());			
		HashSet<String> res = blastn(seqs, resDB, coverage, identity);		
		for (String s:res){
			//System.out.println(s);
			String geneRes = resDB.getRes(s);
			if (geneRes != null && geneRes.length() > 0){
				//Logging.info(line + " : " + geneRes);
				String [] toks = geneRes.split(",");
				for (String tok:toks){
					sos.print(tok + "\t" + resDB.getClass(s) + "\n");
				}					
			}
		}	
		sos.close();
	}

	public static HashSet<String> blastn(ArrayList<Sequence> seqs, ResistanceGeneDB resDB, double minCov, double minID) throws IOException, InterruptedException{
		String blastn = "blastn";

		HashSet<String> res = new HashSet<String>();

		ProcessBuilder pb = new ProcessBuilder(blastn, 
			"-subject", 
			"-",
			"-query", 
			resDB.getSequenceFile(), 
			"-outfmt", 
			"7 qseqid qlen qstart qend sseqid slen sstart send length frames pident nident gaps mismatch score bitscore");
		/////    0      1     2     3    4     5     6     7     8      9      10     11    12    13       14     15
		
		//6 qseqid qlen length pident nident gaps mismatch
		//    0      1    2      3      4      5    6
		Process process = pb.start();

		//Pass on the genome to blastn
		SequenceOutputStream out = new SequenceOutputStream(process.getOutputStream());
		for (Sequence seq:seqs){
			seq.writeFasta(out);
		}
		out.close();

		//Read the output of blastn
		BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));
		String line;

		while ((line = br.readLine()) != null) {
			if (line.startsWith("#"))
				continue;

			String [] toks = line.trim().split("\t");
			int length = Integer.parseInt(toks[8]);
			int qlen   = Integer.parseInt(toks[1]);
			if (minCov * qlen > length)
				continue;

			if (Double.parseDouble(toks[10]) < minID * 100)
				continue;
			//pass			

			res.add(toks[0]);			
		}
		br.close();
		process.waitFor();//Do i need this???

		return res;
	}
}
