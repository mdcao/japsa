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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @author minhduc
 *
 */
public class RealtimeStrainTyping {
    private static final Logger LOG = LoggerFactory.getLogger(RealtimeStrainTyping.class);
	RealtimeStrainTyper typer;
	private double minQual = 0;
	private boolean twoDOnly = false;

	//TODO: make the below private
	ArrayList<Sequence> geneList;
	HashMap<String, Sequence> geneMap;
	HashMap<String, ArrayList<Sequence>> alignmentMap;

	int  currentReadCount = 0;
	long currentBaseCount = 0;
	int  currentReadAligned = 0;


	public RealtimeStrainTyping(int minRead, int minTime, String geneDB, String output) throws IOException{
		typer = new RealtimeStrainTyper(this, geneDB, output);
		typer.setReadPeriod(minRead);
		typer.setTimePeriod(minTime * 1000);

		readGenes(geneDB + File.separatorChar + "geneFam.fasta" );
	}

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

	public void setMinQual(double qual) {
		minQual = qual;
	}

	/**
	 * @param twoOnly the twoOnly to set
	 */
	public void setTwoOnly(boolean twoOnly) {
		this.twoDOnly = twoOnly;
	}

	/**
	 * @param bamFile
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void typing(String bamFile) throws IOException, InterruptedException{
		LOG.info("Strain typing ready at " + new Date());

		alignmentMap = new HashMap<String, ArrayList<Sequence>> ();

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if ("-".equals(bamFile))
			samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			samReader = SamReaderFactory.makeDefault().open(new File(bamFile));

		SAMRecordIterator samIter = samReader.iterator();

		Thread typerThread = new Thread(typer);
		typerThread.start();

		String readName = "";
		//A dummy sequence
		Sequence readSequence = new Sequence(Alphabet.DNA(),1,"");
		while (samIter.hasNext()){
			SAMRecord record = samIter.next();

			if (this.twoDOnly && !record.getReadName().contains("twodim")){
				continue;
			}

			if (!record.getReadName().equals(readName)){
				readName = record.getReadName();

				currentReadCount ++;	
				currentBaseCount += record.getReadLength();

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

			if(record.getMappingQuality() < minQual)
				continue;
			//assert: the read sequence is stored in readSequence with the right direction

			currentReadAligned ++;
			String	geneID = record.getReferenceName();
			if (!geneMap.containsKey(geneID))
				continue;

			synchronized(this){
				int refLength =  geneMap.get(geneID).length();

				ArrayList<Sequence> alignmentList = alignmentMap.get(geneID);
				if (alignmentList == null){
					alignmentList = new ArrayList<Sequence>();
					alignmentMap.put(geneID, alignmentList);
				}
				//put the sequence into alignment list
				//extract the portion of read that align
				Sequence readSeq = HTSUtilities.spanningSequence(record, readSequence, refLength, 0);

				if (readSeq == null){
					LOG.warn("Read sequence is NULL sequence ");
				}else{
					alignmentList.add(readSeq);
				}
			}
		}//while	
		samIter.close();
		samReader.close();
		typer.stopWaiting();
	}

	public static class LCTypingResult implements Comparable<LCTypingResult>{
		String strainID;
		double postProb, l, h;

		@Override
		public int compareTo(LCTypingResult o) {
			return Double.compare(o.postProb, postProb);
		}
	}


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

	protected static double distance (HashSet<String> s1,HashSet<String> s2){		
		int notIn = 0;
		int intersect = 0;

		for (String st:s1){
			if (s2.contains(st))
				intersect ++;
			else notIn ++;
		}
		return (intersect + 0.0) / (notIn + s2.size());
	}

	public static class RealtimeStrainTyper extends RealtimeAnalysis{
		//Set of genes that seen from the sample
		HashSet<String> seenGenes = new HashSet<String>();
		double threshold = 0;
		public SequenceOutputStream datOS = null;
		RealtimeStrainTyping typing;
		PresenceAbsence lcTyping;

		double[] posterior = new double[0];// = lcTyping.calcPosterior();
		double[][] samp = null;// =lcTyping.calcPosterior(1000);
		double[][] ranges = null;// = lcTyping.getRanges(samp, 0.99);

		RealtimeStrainTyper(RealtimeStrainTyping typing, String geneDB, String output) throws IOException{
			super();
			this.typing = typing;

			datOS = SequenceOutputStream.makeOutputStream(output);			
			datOS.print("time\tstep\treads\tbases\tstrain\tprob\tlow\thigh\tgenes\n");
			readKnowProfiles(geneDB + File.separatorChar + "profileFam.dat");
		}

		private void readKnowProfiles(String profileFile) throws IOException{
			String line;
			BufferedReader reader = new BufferedReader (new FileReader(profileFile));		
			ArrayList<RealtimeStrainTyping.GeneProfile> myProfileList = new ArrayList<RealtimeStrainTyping.GeneProfile>(); 

			String currentStrainID = "";
			RealtimeStrainTyping.GeneProfile profile = null;

			while ((line = reader.readLine()) != null){
				if (line.startsWith("#"))
					continue;
				String [] toks = line.trim().split("\t");
				String strainID = toks[0];
				String geneFamID = toks[1];


				if (strainID.equals(currentStrainID)){
					profile.addGene(geneFamID);
				}else{
					profile = new RealtimeStrainTyping.GeneProfile(strainID);
					currentStrainID = strainID;
					profile.addGene(geneFamID);
					myProfileList.add(profile);				
				}
			}
			reader.close();

			/*****************************************************************
			//Checking			
			double [][] mtx = new double[myProfileList.size()][myProfileList.size()];

			for (int i = 0; i < myProfileList.size();i++){
				mtx[i][i] = 0;
				RealtimeStrainTyping.GeneProfile aProfile = myProfileList.get(i);
				for (int j = i + 1; j < myProfileList.size();j++){					
					RealtimeStrainTyping.GeneProfile bProfile = myProfileList.get(j);					 
					mtx[i][j] = mtx[j][i] = 1 - distance(aProfile.genes, bProfile.genes);
				}
			}

			PrintStream out = new PrintStream(new FileOutputStream("infile"));
			out.println(" " + myProfileList.size());
			for (int s = 0; s < myProfileList.size(); s++) {
				String name = myProfileList.get(s).strainID;
				int start = name.length() - 10;
				if (start < 0) 
					start = 0;
				
				out.printf("%-12s ", name.substring(start));
				for (int x = 0; x < myProfileList.size(); x++) {
					out.printf(" %10f ", mtx[s][x]);
				}
				out.println();
			}
			out.close();

			/*****************************************************************/
			LOG.info("There are " + myProfileList.size() +" strains");
			lcTyping = new PresenceAbsence(myProfileList);
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#analysis()
		 */
		@Override
		protected void analysis() {
			// TODO Auto-generated method stub
			try {
				makePresenceTyping(10);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#getCurrentRead()
		 */
		@Override
		protected int getCurrentRead() {			
			return typing.currentReadCount;
		}

		private ArrayList<LCTypingResult> makePresenceTyping(int top) throws IOException, InterruptedException{
			LOG.info("Perform an analysis at " + new Date());
			long step = (lastTime - startTime)/1000;//convert to second

			//int step = typing.currentReadCount;
			//HashSet<String> myGenes = new HashSet<String>();	
			boolean compute = false;

			synchronized(typing){
				for (Sequence gene:typing.geneList){			
					ArrayList<Sequence> alignmentList =  typing.alignmentMap.get(gene.getName());
					//This method use the simple scoreing as it involes ten thousands of genes
					if (alignmentScore2(gene, alignmentList) > threshold){
						//myGenes.add(gene.getName());
						if (!seenGenes.contains(gene.getName())){
							lcTyping.likelihood(100,gene.getName());
							seenGenes.add(gene.getName());
							compute = true;//only need to compute if new evidence is observed
						}				
					}
				}
			}
			LOG.info(timeNow + ": Found " + seenGenes.size() + "  " + compute);

			if (compute){
				posterior = lcTyping.calcPosterior();
				samp =lcTyping.calcPosterior(1000);
				ranges = lcTyping.getRanges(samp, 0.95);
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
				datOS.print(timeNow + "\t" + step  + "\t" + lastReadNumber + "\t" + typing.currentBaseCount + "\t" + lr.strainID + "\t" + lr.postProb +"\t" + (lr.postProb - lr.l) + "\t" + (lr.h -lr.postProb)  +"\t"+seenGenes.size());
				datOS.println();			
			}
			datOS.flush();

			LOG.info("End an analysis at " + new Date());
			return lcT;
		}
		/**
		 * A very simple score where would return the number of reads
		 * @param gene
		 * @param readList
		 * @return
		 */
		private double alignmentScore2(Sequence gene, ArrayList<Sequence> readList){			
			return (readList == null)?0:readList.size();
			//double score = 0;
			//if (readList != null){
			//	score = readList.size();
			//}
			//return score;
		}
		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#close()
		 */
		@Override
		protected void close() {		
			try {
				datOS.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}			
		}
	}
}