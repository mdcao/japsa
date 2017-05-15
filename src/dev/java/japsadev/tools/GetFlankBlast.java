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
 * 25/09/2015 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsadev.tools;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.dev.flankFromBlastn", 
	scriptDesc = "Get sequencing with flanking sequences from blastn results"
	)
public class GetFlankBlast extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(GetFlankBlast.class);

	// Parse result from /sw/blast/current/bin/blastn -db /DataOnline/Data/Bacterial_Genome/Eskape/ftp.ncbi.nlm.nih.gov/blast/db/refseq_genomic -query F0.fasta -num_threads 16 -out F0_refseq.blastn -outfmt 
	// '7 qseqid qlen qstart qend sseqid slen sstart send length frames pident nident gaps mismatch score bitscore'
	 
	//CommandLine cmdLine;
	public GetFlankBlast(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("input", null, "Name of the genome file",true);
		addString("output", "-", "Name of the output file");				
		addString("source", null, "Source sequences",true);

		addStdHelp();
	}

	public static void main(String [] args) throws IOException, InterruptedException{
		GetFlankBlast cmdLine = new GetFlankBlast ();
		args = cmdLine.stdParseLine(args);

		String input = cmdLine.getStringVal("input");
		String output = cmdLine.getStringVal("output");
		String source = cmdLine.getStringVal("source");		

		double ratio = .85;


		ArrayList<Location> locs = new ArrayList<Location>(); 
		HashSet<String> seqSet = new HashSet<String>(); 

		{
			BufferedReader br = new BufferedReader(new FileReader(input));		
			String line = "";
			String currentGene = "";
			double score = 0;		

			int myStart = 0, myEnd = 1;
			String mySeq = "";
			boolean myRev = false;


			LOG.info("Read information");
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#"))
					continue;			
				String [] toks = line.trim().split("\t");
				//Start a new gene
				if (!currentGene.equals(toks[0])){
					if (score > 0){
						Location loc = new Location();
						loc.geneName = currentGene;
						loc.seqName = mySeq;
						loc.start = myStart;
						loc.end = myEnd;
						loc.rev = myRev;
						locs.add(loc);
						seqSet.add(mySeq);
						System.out.println(loc.geneName + " " + loc.seqName + ":" + loc.start + "-" + loc.end + (loc.rev?" -rev":""));
					}
					//starting a new gene				
					score = 0;
					currentGene = toks[0];
				}

				int myScore = Integer.parseInt(toks[14]);
				if (myScore <= score)
					continue;

				double length   = Double.parseDouble(toks[8]);
				double qlen     = Double.parseDouble(toks[1]);
				double slen     = Double.parseDouble(toks[5]);
				double identity = Double.parseDouble(toks[10]);
				if (identity < ratio * 100)
					continue;
				if (length/qlen < ratio)
					continue;			

				int sstart = Integer.parseInt(toks[6]);
				int send = Integer.parseInt(toks[7]);

				if (sstart <= 100 || send <= 100)
					continue;

				if (slen - sstart < 100 || slen - send < 100)
					continue;

				score = myScore;
				if (sstart < send){
					mySeq  = toks[4];
					myStart = sstart - 100;
					myEnd = send + 100;
					myRev = false;
					//location = toks[8] + ":" + (sstart - 100) + "-" + (send + 100);					
				}else{
					mySeq  = toks[4];
					myStart = send - 100;
					myEnd = sstart + 100;
					myRev = true;			
					//location = "-reverse " + toks[8] + ":" + (send - 100) + "-" + (sstart + 100);
				}
			}	
			if (score > 0){
				Location loc = new Location();
				loc.geneName = currentGene;
				loc.seqName = mySeq;
				loc.start = myStart;
				loc.end = myEnd;
				loc.rev = myRev;
				locs.add(loc);
				seqSet.add(mySeq);
				System.out.println(loc.geneName + " " + loc.seqName + ":" + loc.start + "-" + loc.end + (loc.rev?" -rev":""));
			}
			br.close();
		}

		LOG.info("Read source");
		HashMap<String, Sequence> seqMap = new HashMap<String, Sequence>();
		SequenceReader reader = SequenceReader.getReader(source);
		Sequence seq;
		while ( (seq = reader.nextSequence(Alphabet.DNA())) != null){
			if (seqSet.contains(seq.getName())){
				if (seqMap.put(seq.getName(), seq) != null){
					LOG.warn("Sequence " + seq.getName() + " duplicated");
				}	
			}			
		}
		reader.close();
		LOG.info("Extract");
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(output);
		//////////////////////////////////////////////////////
		//	out.print(currentGene + " " + location + "\n");
		for (Location loc:locs){			
			seq = seqMap.get(loc.seqName).subSequence(loc.start - 1, loc.end);
			if (loc.rev){
				seq = Alphabet.DNA.complement(seq);
				seq.setDesc(loc.seqName + ":" + loc.start + "-" + loc.end + " -");				
			}else{
				seq.setDesc(loc.seqName + ":" + loc.start + "-" + loc.end);
			}

			seq.setName(loc.geneName);			
			seq.writeFasta(out);
		}

		out.close();
	}

	static class Location{
		String geneName;
		String seqName;
		int start;
		int end;
		boolean rev;
	}
}

