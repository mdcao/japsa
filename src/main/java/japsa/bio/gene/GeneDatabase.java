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
 * 19 Mar 2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.bio.gene;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.seq.Alphabet.DNA;
import japsa.seq.SequenceOutputStream;
import japsa.util.Logging;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;


public class GeneDatabase implements Iterable<GeneDatabase.GeneFamily>{
	ArrayList<GeneFamily> geneFamilies;
	String dbID = "JSA";	
	
	public GeneDatabase(){
		geneFamilies = new ArrayList<GeneFamily>();
	}	
	
	public String addNewFamily(Sequence seq){
		GeneDatabase.GeneFamily newFam = new GeneDatabase.GeneFamily(size());
		String geneID = newFam.addSequence(seq);
		geneFamilies.add(newFam);
		return geneID;
	}
	
	
	/**
	 * Write database to a fasta file
	 * @param fileName
	 * @throws IOException
	 */
	public void write2File(String fileName, boolean includeAlleles) throws IOException{
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(fileName);
		for (GeneDatabase.GeneFamily family:geneFamilies){
			family.represetationSequence().writeFasta(sos);
			if (includeAlleles){
				for (Sequence seq:family){
					seq.writeFasta(sos);
				}
			}
			
		}		
		sos.close();
	}
	
	
	/**
	 * Read an instance of a database from a file. Assuming the consistency
	 * 
	 * @param fileName
	 * @return
	 * @throws IOException 
	 */
	public static GeneDatabase readDB(String fileName) throws IOException{
		GeneDatabase db = new GeneDatabase();
		SequenceReader reader = SequenceReader.getReader(fileName);
		Alphabet.DNA dna = Alphabet.DNA();
		Sequence seq;		
		while ((seq = reader.nextSequence(dna))!=null){
			String [] toks = seq.getName().split("\\|");
			if (toks.length == 2){
				//rep
				GeneDatabase.GeneFamily newFam = new GeneDatabase.GeneFamily(db.size());
				db.geneFamilies.add(newFam);				
			}else if (toks.length == 3){
				//a gene
				GeneDatabase.GeneFamily fam = db.getFamily(Integer.parseInt(toks[1]));
				fam.geneAlleles.add(seq);
			}else{
				Logging.error("Unknown sequence " + seq.getName());
			}
		}		
		reader.close();
		return db;
	}
	
	/**
	 * Get the family based on the string ID
	 * @param sID
	 * @return
	 */
	public GeneDatabase.GeneFamily getFamily(String sID){
		String [] toks = sID.split("\\|");
		if (toks.length < 2)
			return null;
		
		if (!toks[0].equals(dbID))
			return null;
		
		int fID = -1;
		try{
			fID = Integer.parseInt(toks[1]);
		}catch(Exception e){
			//doing nothing
		}
		
		return geneFamilies.get(fID);
	}
	
	/**
	 * Get the family based on the number ID
	 * @param famID
	 * @return
	 */
	
	public GeneDatabase.GeneFamily getFamily(int famID){
		if (famID >= geneFamilies.size() || famID < 0){
			return null;
		}
		
		return geneFamilies.get(famID);
	}
	
	
	/**
	 * Return the number of families in the database
	 * @return
	 */
	public int size(){
		return geneFamilies.size();
	}
	

	/* (non-Javadoc)
	 * @see java.lang.Iterable#iterator()
	 */
	@Override
	public Iterator<GeneFamily> iterator() {
		return geneFamilies.iterator();

	}
	/////////////////////////////////////////////////////////////////////////
	
	public static class GeneFamily implements Iterable<Sequence>{
		private final int fID;
		private ArrayList<Sequence> geneAlleles;//known instance of this family
		Sequence rep = null;//The representation of this gene family		

		public GeneFamily(int id){
			fID = id;
			geneAlleles = new ArrayList<Sequence>();
		}


		public String familyID(){
			return "JSA|"+fID;
		}

		public Sequence represetationSequence(){
			return rep;
		}

		private void updateRep(Sequence seq){
			if (rep == null || seq.length() > rep.length()){
				rep = seq.clone();
				rep.setName(familyID());
			}
		}

		/**
		 * Add a new sequence to the new family. This will create a new allele
		 * if it is not already in the dababase
		 * @param seq
		 * @return
		 */
		public String addSequence(Sequence seq){
			//This allele already in the database
			for (int i = 0; i < geneAlleles.size(); i++){
				Sequence eSeq = geneAlleles.get(i);
				if (eSeq.match(seq) == 0)
					return eSeq.getName();
				if (eSeq.match(DNA.complement(seq)) == 0)
					return eSeq.getName();
			}

			Sequence nSeq = seq.clone();
			nSeq.setDesc(nSeq.getName() + " " + nSeq.getDesc());			
			nSeq.setName(familyID() + "|" + (geneAlleles.size()));			
			geneAlleles.add(nSeq);

			updateRep(seq);			
			return nSeq.getName();			
		}

		/* (non-Javadoc)
		 * @see java.lang.Iterable#iterator()
		 */
		@Override
		public Iterator<Sequence> iterator() {
			return geneAlleles.iterator();
		}
	}

}