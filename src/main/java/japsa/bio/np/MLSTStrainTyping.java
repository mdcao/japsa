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

import japsa.bio.alignment.ProbFSM.Emission;
import japsa.bio.alignment.ProbFSM.ProbThreeSM;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.HTSUtilities;
import japsa.util.IntArray;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;
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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;


/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.np.mlstStrainTyping", scriptDesc = "Strain typing using MLST system")
public class MLSTStrainTyping {

	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws Exception 
	 * @throws OutOfMemoryError 
	 */
	public static void main(String[] args) throws IOException, InterruptedException{
		/*********************** Setting up script ****************************/
		Deployable annotation = MLSTStrainTyping.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options]", annotation.scriptDesc());		
		/**********************************************************************/		

				
		cmdLine.addString("mlst", "MLST_profiles.txt",  "MLST file");
		cmdLine.addString("bamFile", null,  "The bam file");
		cmdLine.addString("geneFile", null,  "The gene file");		

		cmdLine.addInt("top", 10,  "The number of top strains");		
		cmdLine.addString("msa", "kalign",
				"Name of the msa method, support poa, kalign, muscle and clustalo");
		cmdLine.addString("tmp", "tmp/t",  "Temporary folder");
		cmdLine.addString("hours", null,  "The file containging hours against yields, if set will output acording to tiime");

		cmdLine.addInt("read", 500,  "Number of reads before a typing, NA if timestamp is set");
		cmdLine.addBoolean("twodonly", false,  "Use only two dimentional reads");		

		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/

		String mlst = cmdLine.getStringVal("mlst");
		String bamFile = cmdLine.getStringVal("bamFile");
		String geneFile = cmdLine.getStringVal("geneFile");
		String msa = cmdLine.getStringVal("msa");
		String tmp = cmdLine.getStringVal("tmp");		
		String hours = cmdLine.getStringVal("hours");

		int top = cmdLine.getIntVal("top");		
		int read = cmdLine.getIntVal("read");		

		boolean twodonly = cmdLine.getBooleanVal("twodonly");


		MLSTStrainTyping paTyping = new MLSTStrainTyping();		
		paTyping.msa = msa;
		paTyping.prefix = tmp;			

		paTyping.twoDOnly = twodonly;
		paTyping.readNumber = read;
		if (hours !=null){
			BufferedReader bf = SequenceReader.openFile(hours);
			String line = bf.readLine();//first line
			paTyping.hoursArray = new IntArray();
			paTyping.readCountArray = new IntArray();

			while ((line = bf.readLine())!= null){
				String [] tokens = line.split("\\s");
				int hrs = Integer.parseInt(tokens[0]);
				int readCount = Integer.parseInt(tokens[2]);

				paTyping.hoursArray.add(hrs);
				paTyping.readCountArray.add(readCount);	
			}
		}


		if (paTyping.readNumber < 1)
			paTyping.readNumber = 1;

		paTyping.readGenes(geneFile);
		paTyping.readMLSTProfiles(mlst);
		paTyping.typing(bamFile,  top);	

	}
	/////////////////////////////////////////////////////////////////////////////

	ArrayList<GeneProfile> profileList;
	ArrayList<Sequence> geneList;
	HashMap<String, Sequence> geneMap;
	HashMap<String, ArrayList<Sequence>> alignmentMap;


	int currentReadCount = 0;
	long currentBaseCount = 0;

	String prefix = "tmp";	
	String msa = "kalign";		
	boolean twoDOnly = false;
	int readNumber = 100;
	SequenceOutputStream datOS = null;


	IntArray hoursArray = null;
	IntArray readCountArray = null;
	int arrayIndex = 0;
	
	public MLSTStrainTyping(){
	
	}


	private void readMLSTProfiles(String mlstFile) throws IOException{
		BufferedReader bf = SequenceReader.openFile(mlstFile);
		String line = bf.readLine();
		String [] toks = line.trim().split("\t");
		String [] geneNames = Arrays.copyOfRange(toks, 1, toks.length);
		//mlstProfiles = new ArrayList<int []>();

		profileList = new ArrayList<GeneProfile>();
		//do a check
		for (int i = 0; i < geneNames.length;i++){
			if (!geneNames[i].equals(geneList.get(i).getName())){
				bf.close();
				throw new RuntimeException("Gene file does not match MLST profile");
			}
		}
		while ((line = bf.readLine())!=null){
			toks = line.split("\t");

			GeneProfile profile = new GeneProfile(toks[0]);			
			for (int i = 0; i < geneNames.length;i++){
				profile.addGene(geneNames[i] + "_" + toks[1 + i]);
			}
			profileList.add(profile);
		}//
		bf.close();
	}
	/**
	 * Read genes from gene file to a list + map: for random access
	 * @param geneFile
	 * @throws IOException
	 */
	private void readGenes(String geneFile) throws IOException{
		geneList = SequenceReader.readAll(geneFile, Alphabet.DNA());		 
		geneMap = new HashMap<String, Sequence>();

		for (Sequence gene:geneList){
			geneMap.put(gene.getName(), gene);		
		}
	}


	private void makeMLSTTyping(int top) throws IOException, InterruptedException{
		if (msa.equals("hmm"))
			makeMLSTTyping2(top);
		else
			makeMLSTTyping1(top);

	}
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


	private void makeMLSTTyping2(int top) throws IOException, InterruptedException{		
		//Get consensus read sequence of reads mapped to each gene

		HashMap<String, Double> scoreMap = new HashMap<String, Double>(); 
		for (GeneProfile profile: profileList){
			profile.score = 0;
			for (int i = 0; i < geneList.size(); i++){
				String geneID = geneList.get(i).getName();				
				Logging.info("Compute score for " + profile.genes.get(i));				
				Double scoreObj = scoreMap.get(geneID + "_" + profile.genes.get(i));

				if (scoreObj != null){
					Logging.info(" Found catched " + scoreObj);
					profile.score += scoreObj;					
				}else{
					Logging.info(" FSM alignment :");
					ArrayList<Sequence> alignmentList =  alignmentMap.get(geneID);
					if (alignmentList != null){
						Sequence fa = FastaReader.getReader("geneAlleles/out_" + profile.genes.get(i) + ".fasta").nextSequence(Alphabet.DNA16());
						double thisScore = 0;

						for (Sequence readSeq:alignmentList){
							ProbThreeSM tsm = new ProbThreeSM(fa);
							double cost = 1000000000;

							for (int c = 0; c < 10; c++){
								tsm.resetCount();
								Emission retState = tsm.alignGenerative(readSeq);
								if (cost  <= retState.myCost)
									break;

								//if (cost - retState.myCost < 0.5){
								//	cost = retState.myCost;
								//	break;
								//}
								cost = retState.myCost;
								int emitCount = tsm.updateCount(retState);
								Logging.info("Iter " + c + " : " + emitCount + " states and " + cost + " bits " + readSeq.length() + "bp " + readSeq.getName() + " by " + fa.getName());
								tsm.reEstimate();	
							}
							Logging.info(" Saving " + (readSeq.length() * 2 - cost) + " on " + readSeq.getName() + " by " + fa.getName());						
							if (cost < readSeq.length() * 2){
								thisScore += (readSeq.length() * 2 - cost);
							}
						}
						Logging.info(" Computed " + geneID + "_" + profile.genes.get(i) + " : " + thisScore);
						scoreMap.put(geneID + "_" + profile.genes.get(i), thisScore);
						profile.score += thisScore;					//tsm.showProb();
					}//if not null
				}//else (null)
			}//for i
		}//for profile		

		Collections.sort(profileList);

		if (top > profileList.size())
			top = profileList.size();

		System.out.println("============================================== " + currentReadCount);
		for (int i = 0; i < top;i++){
			GeneProfile profile = profileList.get(i);
			System.out.println(profile.strainID + "\t" + profile.score);					
		}
	}

	private void makeMLSTTyping1(int top) throws IOException, InterruptedException{		
		//Get consensus read sequence of reads mapped to each gene
		HashMap<String, Sequence> consensusMap = new HashMap<String, Sequence>();

		for (int i = 0; i < geneList.size();i++){	
			String geneID = geneList.get(i).getName();
			ArrayList<Sequence> alignmentList =  alignmentMap.get(geneID);
			Sequence consensus = ErrorCorrection.consensusSequence(alignmentList, prefix + "_" + geneID + "_" + this.currentReadCount, msa);			
			if (consensus != null){				
				consensusMap.put(geneID, consensus);
			}
		}		

		HashMap<String, Double> scoreMap = new HashMap<String, Double>(); 
		for (GeneProfile profile: profileList){
			profile.score = 0;
			for (int i = 0; i < geneList.size(); i++){
				String geneID = geneList.get(i).getName();
				Sequence consensus = consensusMap.get(geneID);
				Logging.info("Compute score for " + profile.genes.get(i));
				if (consensus != null){
					Double scoreObj = scoreMap.get(geneID + "_" + profile.genes.get(i));
					if (scoreObj != null){
						profile.score += scoreObj;
						Logging.info(" Found catched " + profile.score);
					}else{
						
						Logging.info(" FSM alignment :");
						Sequence fa = FastaReader.getReader("geneAlleles/out_" + profile.genes.get(i) + ".fasta").nextSequence(Alphabet.DNA16());
						ProbThreeSM tsm = new ProbThreeSM(fa);
						
						double cost = 100000000;						
						for (int c = 0; c < 10; c++){
							tsm.resetCount();
							Emission retState = tsm.alignGenerative(consensus);
							if (cost  <= retState.myCost)
								break;
							
							cost = retState.myCost;
							int emitCount = tsm.updateCount(retState);
							Logging.info("Iter " + c + " : " + emitCount + " states and " + cost + " bits " + consensus.length() + "bp " + consensus.getName() + " by " + fa.getName());
							tsm.reEstimate();	
						}

						
						double thisScore = consensus.length() * 2 - cost; 
						profile.score += thisScore;

						scoreMap.put(geneID + "_" + profile.genes.get(i), thisScore);
						Logging.info(" Computed " + geneID + "_" + profile.genes.get(i) + " : " + thisScore);
						//tsm.showProb();
					}						

				}
			}
		}		

		Collections.sort(profileList);

		if (top > profileList.size())
			top = profileList.size();

		System.out.println("============================================== " + currentReadCount);
		for (int i = 0; i < top;i++){
			GeneProfile profile = profileList.get(i);
			System.out.println(profile.strainID + "\t" + profile.score);						
		}
	}

	private void makeMLSTTyping3(int top) throws IOException, InterruptedException{
		String [] facFiles = new String[geneList.size()];

		//Get consensus read sequence of reads mapped to each gene
		for (int i = 0; i < geneList.size();i++){	
			String geneID = geneList.get(i).getName();
			ArrayList<Sequence> alignmentList =  alignmentMap.get(geneID);
			Sequence consensus = ErrorCorrection.consensusSequence(alignmentList, prefix + "_" + geneID + "_" + this.currentReadCount, msa);


			if (consensus != null){				
				facFiles[i] = prefix + geneID + "_" + this.currentReadCount + "_consensus.fasta";//name of fasta files of reads mapped to the gene				
				consensus.writeFasta(facFiles[i]);
			}else
				facFiles[i] = null;
		}		

		for (GeneProfile profile: profileList){
			profile.score = 0;
			for (int i = 0; i < geneList.size(); i++){
				if (facFiles[i] != null){

					String faAFile = "geneAlleles/out_" + profile.genes.get(i) + ".fasta";
					String needleOut = prefix + profile.genes.get(i) + "_" + this.currentReadCount + "_consensus.needle";;

					File f = new File(needleOut);
					if (!f.exists()){//only run if file not exists (this allele has not done before)
						String cmd = "needle -gapopen 10 -gapextend 0.5 -asequence " 
								+ faAFile + " -bsequence " + facFiles[i] + " -outfile " + needleOut;					

						Logging.info("Running " + cmd);
						Process process = Runtime.getRuntime().exec(cmd);
						process.waitFor();					

						Logging.info("Run'ed " + cmd );
					}
					BufferedReader scoreBf = new BufferedReader(new FileReader(needleOut));
					String scoreLine = null;					

					while ((scoreLine = scoreBf.readLine())!=null){
						String [] scoreToks = scoreLine.split(" ");					
						if (scoreToks.length == 3 && scoreToks[1].equals("Score:")){
							profile.score += Double.parseDouble(scoreToks[2]);
							break;//while
						}					
					}//while
					scoreBf.close();					
				}
			}
		}		

		Collections.sort(profileList);

		if (top > profileList.size())
			top = profileList.size();

		System.out.println("============================================== " + currentReadCount);
		for (int i = 0; i < top;i++){
			GeneProfile profile = profileList.get(i);
			System.out.println(profile.strainID + "\t" + profile.score);						
		}
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
			if (this.twoDOnly && !record.getReadName().contains("twodim")){
				continue;
			}

			if (!record.getReadName().equals(readName)){
				readName = record.getReadName();

				currentReadCount ++;	
				currentBaseCount += record.getReadLength();

				if (hoursArray != null){
					if (arrayIndex < hoursArray.size() && currentReadCount >= this.readCountArray.get(arrayIndex)){
						makeMLSTTyping(top);
						arrayIndex ++;
					}
				}else{				
					if (currentReadCount % readNumber == 0){
						makeMLSTTyping(top);
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

			Sequence readSeq = HTSUtilities.spanningSequence(record, readSequence, refLength,0);

			if (readSeq == null){
				Logging.warn("Read sequence is NULL sequence ");
			}else{
				alignmentList.add(readSeq);
			}
		}//while	
		samIter.close();
		samReader.close();
		makeMLSTTyping(top);
	}

	private static class GeneProfile implements Comparable<GeneProfile>{
		String strainID;
		double score = 0;
		ArrayList<String> genes;

		public GeneProfile(String id){
			strainID = id;
			genes = new ArrayList<String>();
		}

		public void addGene(String geneID){
			genes.add(geneID);
		}

		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(GeneProfile o) {
			double comp = score - o.score;

			if (comp < 0)
				return 1;
			else if (comp > 0)
				return -1;
			else 
				return 0;
		}		
	}

	public static String compare(HashSet<String> s1,HashSet<String> s2){
		String ret = "";
		int count1 = 0,count2=0,count=0;

		for (String st:s1){
			if (s2.contains(st))
				count ++;
			else{ count1++;
			ret = ret + ";" +st;
			}
		}
		ret = ret + "#";

		for (String st:s2){
			if (!s1.contains(st)){			
				count2++;			
				ret = ret + ";" +st;
			}
		}		
		return "Common" + count + "#"+count1+"#"+count2+"#"+ret;
	}
	
}