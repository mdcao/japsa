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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;


public class GeneDatabase extends ArrayList<GeneDatabase.GeneFamily>{
    private static final Logger LOG = LoggerFactory.getLogger(GeneDatabase.class);
	String dbID = "JSA";	
	
	public GeneDatabase(){
	    super();
	}	
	
	public String addNewFamily(Sequence seq){
		GeneDatabase.GeneFamily newFam = new GeneDatabase.GeneFamily(size());
		String geneID = newFam.addSequence(seq);
		add(newFam);
		return geneID;
	}
	
	
	/**
	 * Write database to a fasta file
	 * @param fileName
	 * @throws IOException
	 */
	public void write2File(String fileName, boolean includeAlleles) throws IOException{
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(fileName);
		for (GeneDatabase.GeneFamily family:this){
			Sequence rep = family.represetationSequence();			
			rep.writeFasta(sos);
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
			String [] toks = seq.getName().split("_");
			if (toks.length == 2){
				//rep
				GeneDatabase.GeneFamily newFam = new GeneDatabase.GeneFamily(db.size());
				db.add(newFam);
			}else if (toks.length == 3){
				//a gene
				GeneDatabase.GeneFamily fam = db.getFamily(Integer.parseInt(toks[1]));
				fam.add(seq);
			}else{
				LOG.error("Unknown sequence " + seq.getName());
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
		String [] toks = sID.split("_");
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
		
		return get(fID);
	}
	
	/**
	 * Get the family based on the number ID
	 * @param famID
	 * @return
	 */
	
	public GeneDatabase.GeneFamily getFamily(int famID){
	    //TODO: Check if really need to validate famID
		if (famID >= size() || famID < 0){
			return null;
		}
		
		return get(famID);
	}

	/////////////////////////////////////////////////////////////////////////
	
	public static class GeneFamily extends ArrayList<Sequence>{
		private final int fID;

		int repIndex = -1;		
		String desc = "";
		
		public GeneFamily(int id){
		    super();
			fID = id;
		}

		/**
		 * @return the desc
		 */
		public String getDesc() {
			return desc;
		}


		/**
		 * @param desc the desc to set
		 */
		public void setDesc(String desc) {			
			this.desc = desc;
		}


		public String familyID(){
			return "JSA_"+fID;
		}

		public Sequence represetationSequence(){
			Sequence rep = get(repIndex).clone();
			rep.setDesc(desc + ";index=" +repIndex + ";size=" + size());
			rep.setName(familyID());
			return rep;
		}

		private void updateRep(int newIndex){
			if (repIndex < 0 || get(newIndex).length() > get(repIndex).length()){
				repIndex = newIndex;
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
            //TODO: Need to see if this feature is really needed
			for (int i = 0; i < size(); i++){
				Sequence eSeq = get(i);
				if (eSeq.match(seq) == 0)
					return eSeq.getName();
				if (eSeq.match(DNA.complement(seq)) == 0)
					return eSeq.getName();
			}

			Sequence nSeq = seq.clone();
			nSeq.setDesc(nSeq.getName() + " " + nSeq.getDesc());			
			nSeq.setName(familyID() + "_" + (size()));
			add(nSeq);
			updateRep(size() - 1);
			return nSeq.getName();			
		}
	}
}