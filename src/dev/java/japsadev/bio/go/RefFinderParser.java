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
 * 20 Jul 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsadev.bio.go;

import japsa.bio.gene.GeneDatabase;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.Logging;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * @author minhduc
 *
 */
public class RefFinderParser {	
	public RefFinderParser(String fileName) throws IOException{

	}

	

	/**
	 * @param args
	 * @throws IOException 
	 * 
	 */
	public static void main(String[] args) throws IOException {
		
		//Read in all sequences
		ArrayList<Sequence> seqs = SequenceReader.readAll(args[0],Alphabet.DNA());

		HashMap<String, Sequence> allele2Sequence = new HashMap<String, Sequence>();

		for (Sequence seq:seqs){
			allele2Sequence.put(seq.getName(), seq);
		}

//	/	HashMap<String, String> allele2Groups = new HashMap<String, String>();

		BufferedReader bf = SequenceReader.openFile("togene");
		String line = "";

		HashMap<String, String>  gene2Groups = new HashMap<String, String>();
		HashMap<String, String>  gene2dbGeneID = new HashMap<String, String>();

		GeneDatabase geneDB = new GeneDatabase();


		while ( (line = bf.readLine())!=null){
			String [] toks = line.trim().split("\t");			
			String alleleID = toks[0];
			String agroup = toks[1];
			String geneID = toks[4] + "_" + toks[2];

			if (gene2Groups.containsKey(geneID)){
				String ggroup = gene2Groups.get(geneID);
				
				if (agroup == null){
					Logging.error(alleleID + " xx  " + geneID );
				}

				if (!ggroup.equals(agroup))
					Logging.error( alleleID + " <>  " + geneID );

				GeneDatabase.GeneFamily  dbFamily = geneDB.getFamily(gene2dbGeneID.get(geneID));
				
				if (dbFamily == null)
					Logging.error("Problem finding " + geneID );
				else{
					dbFamily.addSequence(allele2Sequence.get(alleleID));					
				}


			}else{
				gene2Groups.put(geneID, agroup);
				Sequence seq = allele2Sequence.get(alleleID);				 
				String famID = geneDB.addNewFamily(seq);
				gene2dbGeneID.put(geneID, famID);		
				String desc = "geneID=" + geneID + ";dg="+agroup;
				
				geneDB.getFamily(famID).setDesc(desc);
			}
		}

		geneDB.write2File("F.fasta", false);
		geneDB.write2File("A.fasta", true);
		
		Logging.info(gene2Groups.size() + " group ");

		/********************************************************
		for (String relID:odo.relTypes.keySet()){
			System.out.println(relID + "\t" + odo.relTypes.get(relID));
		}

		for (String termID:odo.terms.keySet()){
			GoTerm term = odo.terms.get(termID);
			for (GoRelationship rel:term.relationship){
				if (rel.relType.ID.equals("is_a"))
					System.out.println(termID + "(" + term.name +")\t" + rel.relTerm.ID + "(" + rel.relTerm.name+")");
			}
		}
		/********************************************************/
	}
}
