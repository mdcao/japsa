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
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.np.resistGenes", scriptDesc = "Antibiotic resistance genes identification")
public class ResistanceGene {

	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws Exception 
	 * @throws OutOfMemoryError 
	 */
	public static void main(String[] args) throws IOException, InterruptedException{
		/*********************** Setting up script ****************************/
		Deployable annotation = ResistanceGene.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options]", annotation.scriptDesc());		
		/**********************************************************************/		

		cmdLine.addString("output", "output.dat",  "Output file");
		cmdLine.addString("bamFile", null,  "The bam file");
		//cmdLine.addString("geneFile", null,  "The gene file");
		//cmdLine.addString("figure", null,  "Figure file");
		cmdLine.addString("mcoordFile", null,  "Alignment with dnadiff between the genome (illumina) and genes database");

		cmdLine.addDouble("scoreThreshold", 1.5,  "The alignment score threshold");
		cmdLine.addString("msa", "kalign",
				"Name of the msa method, support poa, kalign, muscle and clustalo");
		cmdLine.addString("global", "needle",
				"Name of the global method, support needle and hmm");
		cmdLine.addString("tmp", "tmp/t",  "Temporary folder");
		cmdLine.addString("hours", null,  "The file containging hours against yields, if set will output acording to time");
		cmdLine.addInt("read", 500,  "Number of reads before a typing, NA if timestamp is set");		
		cmdLine.addInt("timestamp", 0,  "Number of seconds between internval");
		cmdLine.addDouble("il", 0.9,   "Threshold for Illumina");
		cmdLine.addBoolean("twodonly", false,  "Use only two dimentional reads");

		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/

		String output = cmdLine.getStringVal("output");
		String bamFile = cmdLine.getStringVal("bam");
		String mcoordFile = cmdLine.getStringVal("mcoordFile");
		String msa = cmdLine.getStringVal("msa");
		String global = cmdLine.getStringVal("global");

		//String figure = cmdLine.getStringVal("figure");
		String tmp = cmdLine.getStringVal("tmp");
		String hours = cmdLine.getStringVal("hours");

		double scoreThreshold = cmdLine.getDoubleVal("scoreThreshold");				
		int read = cmdLine.getIntVal("read");
		int timestamp = cmdLine.getIntVal("timestamp");

		//double np = cmdLine.getDoubleVal("np");
		double il = cmdLine.getDoubleVal("il");

		boolean twodonly = cmdLine.getBooleanVal("twodonly");

		ResistanceGene paTyping = new ResistanceGene();		
		paTyping.msa = msa;
		paTyping.global = global;

		paTyping.prefix = tmp;
		paTyping.scoreThreshold = scoreThreshold;
		paTyping.twoDOnly = twodonly;
		paTyping.readNumber = read;
		DateFormat df = new SimpleDateFormat("dd/MM/yy HH:mm:ss");
		Logging.info("START : " + df.format(Calendar.getInstance().getTime()));

		if (paTyping.readNumber < 1)
			paTyping.readNumber = 1;

		paTyping.datOS = SequenceOutputStream.makeOutputStream(output);

		paTyping.getGeneClassInformation(mcoordFile);
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



		//paTyping.readGenes(geneFile);
		//paTyping.getAROMap(mcoordFile);		

		//paTyping.npThreshold = np;
		paTyping.ilThreshold = il;		



		//paTyping.datOS.print("step\treads\tbases");
		//for (int i = 0; i < paTyping.functionNames.size();i++){
		//	paTyping.datOS.print('\t');
		//	paTyping.datOS.print(paTyping.functionNames.get(i));
		//}
		//paTyping.datOS.println();
		paTyping.typing(bamFile);

		paTyping.datOS.close();


	}

	//double npThreshold = 0.9;
	double ilThreshold = 0.9;

	/////////////////////////////////////////////////////////////////////////////


	HashMap<String, HashSet<String>> aroMap;//map a gene to function/annotations
	HashMap<String, ArrayList<Sequence>> alignmentMap;
	HashMap<String, ArrayList<SAMRecord>> samMap;

	HashSet<String> targetGenes;
	HashSet<String> targetGroup;
	HashSet<String> targetPrimaryGroup;


	HashMap<String, Double> maxScores = new HashMap<String, Double>(); 

	int currentReadCount = 0;
	long currentBaseCount = 0;

	String prefix = "tmp";	
	String msa = "kalign";
	String global = "needle";

	double scoreThreshold = 2;

	//HashMap<String,XYSeries> plotData = new HashMap<String,XYSeries>();
	//String figureFile = "out.png";
	boolean twoDOnly = false;

	String functionFile = "Function.txt";
	String gene2FunctionFile = "Gene2Function.txt";

	ArrayList<String> functions = new ArrayList<String>();	
	ArrayList<String> functionNames = new ArrayList<String>();
	int readNumber = 100;
	SequenceOutputStream datOS = null;

	IntArray hoursArray = null;
	IntArray readCountArray = null;
	int arrayIndex = 0;
	
	long firstReadTime = 0;

	public ResistanceGene(){
	}

	HashMap<String, ArrayList<Sequence>> familyMap;

	HashMap<String, String> gene2Group;
	HashMap<String, String> gene2PrimaryGroup;
	HashSet<String> groupSet, primaryGroupSet, geneSet;

	ArrayList<Sequence> geneList;
	HashMap<String, Sequence> geneMap;
	

	private void getGeneClassInformation(String mcoordFile) throws IOException{
		ArrayList<Sequence> genes = FastaReader.readAll("geneAlleles90.fasta", Alphabet.DNA());
		HashMap<String, Sequence> myMap = new HashMap<String, Sequence>();

		gene2Group = new HashMap<String, String>();	
		gene2PrimaryGroup = new HashMap<String, String>();		

		familyMap  = new HashMap<String, ArrayList<Sequence>>(); 
		geneMap    = new HashMap<String, Sequence>();
		geneList   = new ArrayList<Sequence> ();

		for (Sequence g:genes){
			myMap.put(g.getName(), g);
		}
		
		
		groupSet = new HashSet<String>();
		primaryGroupSet = new HashSet<String>();
		geneSet = new HashSet<String>();
		

		String fn = "GeneClassification.csv";

		BufferedReader bf = SequenceReader.openFile(fn);
		String line = bf.readLine();
		while ((line = bf.readLine())!= null){
			if (line.startsWith("#"))
				continue;

			String [] toks = line.split("\t");
			if (!toks[0].equals("0"))
				continue;							

			String famID = "JSA_"+toks[5];
			String geneID = "JSA_"+toks[5] + "_" + toks[6];			
			Sequence gene = myMap.remove(geneID);

			if (gene != null){
				ArrayList<Sequence> fam = familyMap.get(famID);
				if (fam == null){//fam not yet entered
					fam = new ArrayList<Sequence>();
					familyMap.put(famID, fam);			
					Sequence famSeq = myMap.get(famID);
					if(famSeq == null){
						Logging.exit(" Gene Family " + famID + " not found",1);
					}
					geneList.add(famSeq);
					geneMap.put(famID, famSeq);
					gene2Group.put(famID, toks[3]);
					gene2PrimaryGroup.put(famID, toks[2]);

				}
				fam.add(gene);
				gene2Group.put(geneID, toks[3]);
				gene2PrimaryGroup.put(geneID, toks[2]);
			}
		}
		bf.close();
					
		
		groupSet.addAll(gene2Group.values());
		primaryGroupSet.addAll(gene2PrimaryGroup.values());

		//Read the target genes (get from Illumina)
		targetGenes = new HashSet<String>();	
		targetGroup = new HashSet<String>();
		targetPrimaryGroup = new HashSet<String>();
		bf = new BufferedReader (new FileReader(mcoordFile));


		while ((line = bf.readLine()) != null){
			String [] toks = line.trim().split("\t");			
			String [] tt = toks[12].split("_");
			String tGeneID = tt[0] + "_" + tt[1];

			if (Double.parseDouble(toks[10]) >= this.ilThreshold * 100 
					&& Double.parseDouble(toks[6]) >= this.ilThreshold * 100 
					&& geneMap.containsKey(tGeneID)){
				targetGenes.add(tGeneID);
				targetGroup.add(gene2Group.get(tGeneID));
				targetPrimaryGroup.add(gene2PrimaryGroup.get(tGeneID));				
			}
		}		

		datOS.print("#Target : " + targetGenes.size() + "\n");
		for (String tGeneID: targetGenes){
			datOS.print("#TG " + tGeneID + " " + gene2Group.get(tGeneID) + " " + gene2PrimaryGroup.get(tGeneID));
			datOS.println();
		}

		for (String group:targetGroup){
			datOS.print("#TC " + group);
			datOS.println();
		}
		for (String group:targetPrimaryGroup){
			datOS.print("#TP " + group);
			datOS.println();
		}

		datOS.flush();
		bf.close();
	}


	int [] targetProfile = null;
	public double  checkNeedle(String consensusFile, Sequence gene) throws IOException, InterruptedException{
		//Needle the gene
		String geneID = gene.getName();
		String faAFile = "geneAlleles/out_" + geneID + ".fasta";
		String needleOut = prefix + geneID + "_" + this.currentReadCount + "_consensus.needle";

		String cmd = "needle -gapopen 10 -gapextend 0.5 -asequence " 
				+ faAFile + " -bsequence " + consensusFile + " -outfile " + needleOut;
		Logging.info("Running " + cmd);
		Process process = Runtime.getRuntime().exec(cmd);
		process.waitFor();		
		Logging.info("Run'ed " + cmd );

		BufferedReader scoreBf = new BufferedReader(new FileReader(needleOut));
		String scoreLine = null;
		double score = 0;
		while ((scoreLine = scoreBf.readLine())!=null){
			String [] scoreToks = scoreLine.split(" ");					
			if (scoreToks.length == 3 && scoreToks[1].equals("Score:")){
				score += Double.parseDouble(scoreToks[2]);
				break;//while
			}					
		}//while
		scoreBf.close();
		return score / gene.length();
	}

	public double checkHMM(Sequence consensus, Sequence gene){
		
		if (gene.length() > 2700 || consensus.length() > 4000 || gene.length() * consensus.length() > 6000000){
			Logging.info("SKIP " + gene.getName() + " " + gene.length() + " vs " + consensus.length());			
			return 0;
		}

		ProbThreeSM tsmF = new ProbThreeSM(gene);
		double cost = 100000000;						
		for (int c = 0; c < 10; c++){
			tsmF.resetCount();
			Emission retState = tsmF.alignGenerative(consensus);
			if (cost  <= retState.myCost)
				break;//for c

			cost = retState.myCost;
			int emitCount = tsmF.updateCount(retState);
			Logging.info("Iter " + c + " : " + emitCount + " states and " + cost + " bits " + consensus.length() + "bp " + consensus.getName() + " by " + gene.getName());
			tsmF.reEstimate();	
		}				
		return (consensus.length() * 2 - cost) / gene.length();
	}

	private void addPreditedGene(String geneID){
		pGenes.add(geneID);
		pGroup.add(gene2Group.get(geneID));
		pPrimaryGroup.add(gene2PrimaryGroup.get(geneID));		
	}
	HashSet<String> pGenes = new HashSet<String>();
	HashSet<String> pGroup = new HashSet<String>();
	HashSet<String> pPrimaryGroup = new HashSet<String>();
	private void antiBioticsProfile() throws IOException, InterruptedException{


		int step = currentReadCount;
		if (hoursArray != null) 
			step = hoursArray.get(arrayIndex);

		Logging.info("STEP " + step);
		//Get list of genes from my
		for (Sequence gene:geneList){
			String geneID = gene.getName();
			//TODO: to remove
			if (pPrimaryGroup.contains(gene2PrimaryGroup.get(geneID)))
				continue;
			
			if (pGenes.contains(geneID))
				continue;

			ArrayList<Sequence> alignmentList =  alignmentMap.get(geneID);
			Sequence 
				consensus = ErrorCorrection.consensusSequence(alignmentList, prefix + "_" + geneID + "_" + this.currentReadCount, msa);

			if (consensus == null){
				//Not consider this gene at all
				continue;//gene
			}

			if (global.equals("hmm")){
				//	double cutoff = 0.6;				
				{	
					double score = checkHMM(consensus, gene);
					
					//if (score == 0 && gene2PrimaryGroup.get(geneID).equals("Daptomycin-Rifampin")){
					//	score = 0.8;
					//	Logging.info("XXXX " + geneID);
					//}
					
					Logging.info("REF: " + geneID + " " + score + " " + gene.length() + " " + consensus.length()+ " " + gene2Group.get(geneID)+ " " + gene2PrimaryGroup.get(geneID));

					Double oldScore = maxScores.get(gene2PrimaryGroup.get(geneID));
					if (oldScore == null || oldScore < score){
						maxScores.put(gene2PrimaryGroup.get(geneID), score);
					}					

					Logging.info("SGF: " + score + " " + geneID + " " + targetGenes.contains(geneID));
					Logging.info("SCF: " + score + " " + gene2Group.get(geneID) + " " + targetGroup.contains(gene2Group.get(geneID)));
					Logging.info("SPF: " + score + " " + gene2PrimaryGroup.get(geneID) + " " + targetPrimaryGroup.contains(gene2PrimaryGroup.get(geneID)));
					
					if (score >= scoreThreshold){
						addPreditedGene(geneID);
						Logging.info("ADDF " + geneID + " " + gene2Group.get(geneID)+ " " + gene2PrimaryGroup.get(geneID) + " " + step + " " + geneID);						
						continue;//for gene
					}
				}

				/*************************************************/		
				ArrayList<Sequence> alleleList = familyMap.get(geneID);
				for (Sequence allele:alleleList){
					double score = checkHMM(consensus, allele);
					Logging.info("REA: " + allele.getName() + " " + score  + " " + allele.length() + " " + consensus.length()+ " " + gene2Group.get(geneID)+ " " + gene2PrimaryGroup.get(geneID));					
					Double oldScore = maxScores.get(gene2PrimaryGroup.get(geneID));
					if (oldScore == null || oldScore < score){
						maxScores.put(gene2PrimaryGroup.get(geneID), score);
					}

					
					Logging.info("SGA: " + score + " " + geneID + " " + targetGenes.contains(geneID));
					Logging.info("SCA: " + score + " " + gene2Group.get(geneID) + " " + targetGroup.contains(gene2Group.get(geneID)));
					Logging.info("SPA: " + score + " " + gene2PrimaryGroup.get(geneID) + " " + targetPrimaryGroup.contains(gene2PrimaryGroup.get(geneID)));
					
					
					if (score >= scoreThreshold){
						addPreditedGene(geneID);
						Logging.info("ADDA " + geneID + " " + gene2Group.get(geneID)+ " " + gene2PrimaryGroup.get(geneID) + " " + step + " " + allele.getName());
						break;//for allele
					}
				}

			}else{
				String consensusFile = prefix + "consensus" + geneID + "_" + this.currentReadCount + ".fasta"; 
				consensus.writeFasta(consensusFile);				
				{	
					double score = checkNeedle(consensusFile, gene);
					Logging.info("REF: " + geneID + " " + score + " " + gene.length() + " " + consensus.length()+ " " + gene2Group.get(geneID)+ " " + gene2PrimaryGroup.get(geneID));					
					Double oldScore = maxScores.get(gene2PrimaryGroup.get(geneID));
					if (oldScore == null || oldScore < score){
						maxScores.put(gene2PrimaryGroup.get(geneID), score);
					}
					
					
					Logging.info("SGF: " + score + " " + geneID + " " + targetGenes.contains(geneID));
					Logging.info("SCF: " + score + " " + gene2Group.get(geneID) + " " + targetGroup.contains(gene2Group.get(geneID)));
					Logging.info("SPF: " + score + " " + gene2PrimaryGroup.get(geneID) + " " + targetPrimaryGroup.contains(gene2PrimaryGroup.get(geneID)));
					if (score >= scoreThreshold){					
						Logging.info("ADDF " + geneID + " " + gene2Group.get(geneID) + " " + gene2PrimaryGroup.get(geneID) + " " + step + " " + geneID);
						addPreditedGene(geneID);
						continue;//for gene
					}
				}

				//Needle every allele
				ArrayList<Sequence> alleleList = familyMap.get(geneID);			
				for (Sequence allele:alleleList){					
					double score = checkNeedle(consensusFile, allele);					
					Logging.info("REA: " + geneID + " " + score + " " + allele.length() + " " + consensus.length() + " " + gene2Group.get(geneID)+ " " + gene2PrimaryGroup.get(geneID) );
					
					Double oldScore = maxScores.get(gene2PrimaryGroup.get(geneID));
					if (oldScore == null || oldScore < score){
						maxScores.put(gene2PrimaryGroup.get(geneID), score);
					}									
					
					Logging.info("SGA: " + score + " " + geneID + " " + targetGenes.contains(geneID));
					Logging.info("SCA: " + score + " " + gene2Group.get(geneID) + " " + targetGroup.contains(gene2Group.get(geneID)));
					Logging.info("SPA: " + score + " " + gene2PrimaryGroup.get(geneID) + " " + targetPrimaryGroup.contains(gene2PrimaryGroup.get(geneID)));
					
					
					if (score >= scoreThreshold){
						Logging.info("ADDA " + geneID + " " + gene2Group.get(geneID)+ " " + gene2PrimaryGroup.get(geneID) + " " + step  + " " + allele.getName());
						addPreditedGene(geneID);			
						break;//for allele
					}				
				}//for allele
			}
		}

		Logging.info("===Found " + pGenes.size() + " vs " + geneList.size() + "  " + alignmentMap.size() + " " + targetGenes.size());

		/********************************************/
		int TP = 0, FN = 0;
		for (String geneID: targetGenes){
			if (pGenes.contains(geneID))
				TP ++;
			else{
				FN ++;
				//System.out.println("FN " + geneID);

			}
		}//for geneID

		int FP = pGenes.size() - TP;

		double precision = (TP * 1.0) / (TP + FP);
		double recall = (TP * 1.0) / (TP + FN);			
		double f1 = 2 * precision * recall/ (precision + recall);
		/********************************************/
		datOS.print("O\t" + step + "\t" + currentReadCount + "\t" + currentBaseCount + "\t" +  f1 + "\t" + recall + "\t" + precision + "\t" + pGenes.size() + "\t" + targetGenes.size());
		datOS.println();

		datOS.print("G\t" + step + "\t" + currentReadCount + "\t" + currentBaseCount + "\t" + stats(targetGenes, pGenes, familyMap.keySet() ) + "\t" + pGenes.size() + "\t" + targetGenes.size());
		datOS.println();

		datOS.print("C\t" + step + "\t" + currentReadCount + "\t" + currentBaseCount + "\t" + stats(targetGroup, pGroup, groupSet) + "\t" + pGroup.size() + "\t" + targetGroup.size());
		datOS.println();

		datOS.print("P\t" + step + "\t" + currentReadCount + "\t" + currentBaseCount + "\t" + stats(targetPrimaryGroup, pPrimaryGroup, primaryGroupSet) + "\t" + pPrimaryGroup.size() + "\t" + targetPrimaryGroup.size());		
		datOS.println();

		//System.out.printf("%d %6.4f %6.4f %6.4f %d %d, %d\n",step,f1,recall, precision, myGenes.size(), targetGenes.size(), currentReadCount);		

		for (String x:pGenes){
			datOS.print("#PG\t" + step + "\t" + x + "\t" + targetGenes.contains(x));
			datOS.println();
		}

		for (String x:pGroup){
			datOS.print("#PC\t" + step + "\t" + x + "\t" + targetGroup.contains(x));
			datOS.println();
		}
		for (String x:pPrimaryGroup){
			datOS.print("#PP\t" + step + "\t" + x + "\t" + targetPrimaryGroup.contains(x));
			datOS.println();
		}


		datOS.println();
		datOS.flush();
		/****************************************************************
		int [] myProfile = new int[functions.size()];
		for (String geneID: myGenes){			
			HashSet<String> fList = aroMap.get(geneID);			
			if (fList != null){				
				for (String geneFunction: fList){
					int index = functions.indexOf(geneFunction);		
					if (index >= 0)
						myProfile[index] ++;
				}
			}		 
		}

		//antibioticProfiles.add(myProfile);
		datOS.print(step + "\t" + currentReadCount + "\t" + currentBaseCount);
		for (int i =0 ;i < myProfile.length;i++){
			datOS.print("\t" + myProfile[i]);
		}
		datOS.println();
		datOS.flush();
		/****************************************************************/
	}
	/**
	 * @param bamFile
	 * @param geneFile
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void typing(String bamFile) throws IOException, InterruptedException{		
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

			
			if (this.twoDOnly && !record.getReadName().contains("twodim")){
				continue;
			}

			if (!record.getReadName().equals(readName)){
				readName = record.getReadName();

				currentReadCount ++;	
				currentBaseCount += record.getReadLength();

				if (hoursArray != null){
					if (arrayIndex < hoursArray.size() && currentReadCount >= this.readCountArray.get(arrayIndex)){
						antiBioticsProfile();
						arrayIndex ++;
					}
				}else{				
					if (currentReadCount % readNumber == 0){
						antiBioticsProfile();
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

			Sequence readSeq = HTSUtilities.spanningSequence(record, readSequence, refLength, 20);

			if (readSeq == null){
				Logging.warn("Read sequence is NULL sequence ");
			}else{
				alignmentList.add(readSeq);
			}
		}//while	
		samIter.close();
		samReader.close();
		antiBioticsProfile();
		//makeStackedChart();
		
		for (String pg:maxScores.keySet()){			
			double score = maxScores.get(pg);
			Logging.info("SCORE: " + pg + " " + score + " " + targetPrimaryGroup.contains(pg));
		}

		DateFormat df = new SimpleDateFormat("dd/MM/yy HH:mm:ss");
		Logging.info("END : " + df.format(Calendar.getInstance().getTime()));
	}

	/************************************************************************
	ArrayList<int[]> antibioticProfiles = new ArrayList<int[]>(); 
	private void makeStackedChart() throws IOException{
		if (antibioticProfiles.size() <=0)
			return;

		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int r = 0; r < functions.size(); r++) {
			String rowKey = functions.get(r);
			for (int c = 0; c < antibioticProfiles.size(); c++) {
				Double columnKey = new Double(100 * c);
				dataset.addValue(antibioticProfiles.get(c)[r], rowKey, columnKey);
			}
		}


		final JFreeChart chart = ChartFactory.createStackedAreaChart(
				"Antibiotic gene profile",      // chart title
				"Read number",                // domain axis label
				"Gene number",                   // range axis label
				dataset,                   // data
				PlotOrientation.VERTICAL,  // orientation
				true,                      // include legend
				true,
				false
				);
		chart.setBackgroundPaint(Color.white);


		final CategoryPlot plot = (CategoryPlot) chart.getPlot();


		//plot.setForegroundAlpha(0.5f);
		plot.setBackgroundPaint(Color.WHITE);		
		//plot.setDomainGridlinePaint(Color.white);
		//plot.setRangeGridlinePaint(Color.white);

		final CategoryAxis domainAxis = plot.getDomainAxis();
		domainAxis.setLowerMargin(0.0);
		domainAxis.setUpperMargin(0.0);

		// change the auto tick unit selection to integer units only...
		final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
		//rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());		


		File chartFile = new File(figureFile);
		ChartUtilities. saveChartAsPNG(chartFile , chart , 640 , 480 );
	}
/
	 * @param <K>************************************************************************/
	private <K> String stats(HashSet<K> a, HashSet<K> p, Collection<K> universe){
		int TP = 0, FP = 0, FN = 0;
		for (K s:p){
			if (a.contains(s))
				TP ++;//True positive
			else{
				FP ++;
			}			
		}
		FN = a.size() - TP;
		
		int TN = universe.size() - TP - FP - FN;
		
		double precision = (TP * 1.0) / (TP + FP);
		double recall = (TP * 1.0) / (TP + FN);
		double specificity = ( TN * 1.0 ) / (TN + FP);
		
		
		double f1 = 2 * precision * recall/ (precision + recall);
		double accuracy = (TP + TN + 0.0)/(universe.size());
		
		return f1 +"\t" + precision + "\t" + recall + "\t" + specificity + "\t" + accuracy + "\t" + TP + "\t" + TN + "\t" + FP + "\t" + FN;		
	}


}