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
 * 07/09/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.np;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.HTSUtilities;
import japsa.util.IntArray;
import japsa.util.Logging;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @author minhduc
 *
 */
public class GeneStrainTyping {
	//TODO: make the below private
	HashSet<String> addedGenes = new HashSet<String>(); 
	PresenceAbsence lcTyping;

	ArrayList<Sequence> geneList;
	HashMap<String, Sequence> geneMap;

	HashMap<String, ArrayList<Sequence>> alignmentMap;

	HashSet<String> targetGenes;

	public int readNumber = 100;
	public SequenceOutputStream datOS = null;


	int currentReadCount = 0;
	long currentBaseCount = 0;
	int currentReadAligned = 0;

	public IntArray hoursArray = null;
	public IntArray readCountArray = null;
	int arrayIndex = 0;

	long startTime;
	long firstReadTime = 0;	
	public int timestamp = 5000;


	public GeneStrainTyping(){		
		startTime = System.currentTimeMillis();		
	}



	/**
	 * Analysis of drug resistance profile
	 * @param mcoordFile
	 * @throws IOException
	 */


	/**
	 * Read genes from gene file to a list + map: for random access
	 * TODO: make private
	 * @param geneFile
	 * @throws IOException
	 */
	public void readGenes(String geneFile) throws IOException{
		geneList = SequenceReader.readAll(geneFile, Alphabet.DNA());		 
		geneMap = new HashMap<String, Sequence>();

		for (Sequence gene:geneList){
			geneMap.put(gene.getName(), gene);		
		}
	}


	/**
	 * Read the gene profile of all known strains (for strain typing)
	 * TODO:make private 
	 * @param profileFile
	 * @throws IOException
	 */
	public void readKnowProfiles(String profileFile) throws IOException{
		String line;
		BufferedReader reader = new BufferedReader (new FileReader(profileFile));		
		ArrayList<GeneProfile> myProfileList = new ArrayList<GeneProfile>(); 

		String currentStrainID = "";
		GeneProfile profile = null;

		while ((line = reader.readLine()) != null){
			if (line.startsWith("#"))
				continue;
			String [] toks = line.trim().split("\t");
			String strainID = toks[0];
			String geneFamID = toks[1];
			

			if (strainID.equals(currentStrainID)){
				profile.addGene(geneFamID);
			}else{
				profile = new GeneProfile(strainID);
				currentStrainID = strainID;
				profile.addGene(geneFamID);
				myProfileList.add(profile);				
			}
		}		
		reader.close();
		double threshold = 0.999;
		//Checking
		HashSet<Integer> removeList = new HashSet<Integer>();   
		for (int i = 0; i < myProfileList.size();i++){
			if (removeList.contains(i))
				continue;

			GeneProfile aProfile = myProfileList.get(i);
			for (int j = i + 1; j < myProfileList.size();j++){
				if (removeList.contains(j))
					continue;

				GeneProfile bProfile = myProfileList.get(j);		
				double distance = distance(aProfile.genes, bProfile.genes);

				if (distance > threshold){
					Logging.warn("CHECK  " + aProfile.strainID + " similar to " + bProfile.strainID + " " + distance);
					removeList.add(j);
					continue;
				}

				//if (aProfile.genes.equals(bProfile.genes)){
				//	Logging.warn("CHECK  " + aProfile.strainID + " equals " + bProfile.strainID);
				//	removeList.add(j);
				//	continue;
				//}				
				//if (aProfile.genes.containsAll(bProfile.genes)){
				//	Logging.warn("CHECK1 " + aProfile.strainID + "(" + aProfile.genes.size() + ") contains all of " + bProfile.strainID + " (" + bProfile.genes.size() + ")");
				//}
				//if (bProfile.genes.containsAll(aProfile.genes)){
				//	Logging.warn("CHECK2 " + bProfile.strainID + "(" + bProfile.genes.size() + ") contains all of " + aProfile.strainID + " (" + aProfile.genes.size() + ")");
				//}
			}
		}

		ArrayList<GeneProfile> profileList = new ArrayList<GeneProfile>();
		for (int i = 0; i< myProfileList.size();i++){
			if (!removeList.contains(i))
				profileList.add(myProfileList.get(i));
		}
		Logging.info("There are " + myProfileList.size() +" strains");
		lcTyping = new PresenceAbsence(profileList);
	}




	double[] posterior = new double[0];// = lcTyping.calcPosterior();
	double[][] samp = null;// =lcTyping.calcPosterior(1000);
	double[][] ranges = null;// = lcTyping.getRanges(samp, 0.99);


	/**
	 * Compute the alignment score between a gene and a list of (errornous) reads
	 * that were aligned to the gene. The algorithm is:
	 *  - Get all read sequences (which were previously trimmed)
	 *  - Call a MSA method to alignment them
	 *  - Make a consensus sequence of those reads
	 *  - Align the consensus sequence to the gene sequence using needle
	 *  - Return the needle alignment score
	 * @param gene : the gene sequence
	 * @param readList: an array list of read sequences aligned to this gene

	 * @throws IOException 
	 * @throws InterruptedException 
	 */


	private double alignmentScore2(Sequence gene, ArrayList<Sequence> readList){
		double score = 0;
		if (readList != null){
			score = readList.size();
		}
		return score;
	}


	HashSet<String> mentionedStrain = new HashSet<String>(); 
	double threshold = 0;

	private ArrayList<LCTypingResult> makePresenceTyping(int top) throws IOException, InterruptedException{

		int step = currentReadCount;
		if (hoursArray != null) 
			step = hoursArray.get(arrayIndex);

		//HashSet<String> myGenes = new HashSet<String>();	
		boolean compute = false;
		for (Sequence gene:geneList){			
			ArrayList<Sequence> alignmentList =  alignmentMap.get(gene.getName());
			//This method use the simple scoreing as it involes ten thousands of genes
			if (alignmentScore2(gene, alignmentList) > threshold){
				//myGenes.add(gene.getName());

				if (!addedGenes.contains(gene.getName())){
					lcTyping.likelihood(100,gene.getName());
					addedGenes.add(gene.getName());
					compute = true;//only need to compute if new evidence is observed
				}				
			}
		}
		Logging.info(step + ": Found " + addedGenes.size() + "  " + compute);

		if (compute){
			posterior = lcTyping.calcPosterior();
			samp =lcTyping.calcPosterior(1000);
			ranges = lcTyping.getRanges(samp, 0.99);
		}
		ArrayList<LCTypingResult> lcT = new ArrayList<LCTypingResult>(); 
		for(int i=0; i<posterior.length; i++){
			LCTypingResult lts = new LCTypingResult();
			lts.strainID = lcTyping.spl[i].species;
			lts.postProb = posterior[i];
			lts.l = ranges[i][0];
			lts.h = ranges[i][1];
			lcT.add(lts);
		}
		Collections.sort(lcT);

		for (int i = 0; i < top && i < lcT.size();i++){
			LCTypingResult lr  = lcT.get(i);

			if (lr.postProb < 0.010)
				break;
			datOS.print(step + "\t" + currentReadCount + "\t" + currentBaseCount + "\t" + lr.strainID + "\t" + lr.postProb +"\t" + (lr.postProb - lr.l) + "\t" + (lr.h -lr.postProb)  +"\t"+addedGenes.size());
			datOS.println();			
		}
		datOS.flush();

		return lcT;
	}


	/**
	 * @param bamFile
	 * @param geneFile
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void typing(String bamFile, int top) throws IOException, InterruptedException{		
		alignmentMap = new HashMap<String, ArrayList<Sequence>> ();

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if ("-".equals(bamFile))
			samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			samReader = SamReaderFactory.makeDefault().open(new File(bamFile));

		SAMRecordIterator samIter = samReader.iterator();				

		String readName = "";
		//A dummy sequence
		Sequence readSequence = new Sequence(Alphabet.DNA(),1,"");
		while (samIter.hasNext()){
			SAMRecord record = samIter.next();
			if (firstReadTime <=0)
				firstReadTime = System.currentTimeMillis();

			//if (this.twoDOnly && !record.getReadName().contains("twodim")){
			//	continue;
			//}

			if (!record.getReadName().equals(readName)){
				readName = record.getReadName();

				currentReadCount ++;	
				currentBaseCount += record.getReadLength();


				if (hoursArray != null){
					if (arrayIndex < hoursArray.size() && currentReadCount >= this.readCountArray.get(arrayIndex)){						
						makePresenceTyping(top);
						arrayIndex ++;
					}
				}else{				
					if (currentReadCount % readNumber == 0){
						makePresenceTyping(top);
					}
				}

				//Get the read
				if (!record.getReadUnmappedFlag()){
					readSequence = new Sequence(Alphabet.DNA(), record.getReadString(), readName);
					if (record.getReadNegativeStrandFlag()){
						readSequence = Alphabet.DNA.complement(readSequence);
						readSequence.setName(readName);
					}
				}
			}

			if (record.getReadUnmappedFlag())
				continue;			
			//assert: the read sequence is stored in readSequence with the right direction

			currentReadAligned ++;
			String	geneID = record.getReferenceName();
			if (!geneMap.containsKey(geneID))
				continue;

			int refLength =  geneMap.get(geneID).length();


			ArrayList<Sequence> alignmentList = alignmentMap.get(geneID);
			if (alignmentList == null){
				alignmentList = new ArrayList<Sequence>();
				alignmentMap.put(geneID, alignmentList);
			}
			//put the sequence into alignment list

			Sequence readSeq = HTSUtilities.spanningSequence(record, readSequence, refLength, 0);

			if (readSeq == null){
				Logging.warn("Read sequence is NULL sequence ");
			}else{
				alignmentList.add(readSeq);
			}
		}//while	
		samIter.close();
		samReader.close();

		makePresenceTyping(top);

	}


	public static class LCTypingResult implements Comparable<LCTypingResult>{
		String strainID;
		double postProb, l, h;
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(LCTypingResult o) {
			return Double.compare(o.postProb, postProb);
		}
	}


	/**
	 * Set up: -- read gff files, extract gene sequences, generate profile for each strain
	 * 
	 * @param file
	 * @param out
	 * @param profile
	 * @throws IOException
	 */

	public static class GeneProfile implements Comparable<GeneProfile>{
		String strainID;
		double score = 0.0;
		double f1 = 0, precision, recall;
		HashSet<String> genes;

		public GeneProfile(String id){
			strainID = id;
			genes = new HashSet<String>();
		}

		public void addGene(String geneID){
			genes.add(geneID);
		}

		public HashSet<String>  getGeneList(){
			return genes;
		}

		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(GeneProfile o) {
			return Double.compare(o.score, score);
		}

		public String strainID(){
			return strainID;
		}
	}


	private static double distance (HashSet<String> s1,HashSet<String> s2){		
		int count= 0;
				
		for (String st:s1){
			if (s2.contains(st))
				count ++;
		}
		return count *2.0 / (s1.size() + s2.size());
	}




}