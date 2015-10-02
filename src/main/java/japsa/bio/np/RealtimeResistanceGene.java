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
import japsa.bio.alignment.ProbFSM.ProbOneSM;
import japsa.bio.alignment.ProbFSM.ProbThreeSM;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.HTSUtilities;
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
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;


/**
 * 
 * @author minhduc
 *
 */

public class RealtimeResistanceGene {
	ResistanceGeneFinder resistFinder;		

	private HashMap<String, ArrayList<Sequence>> alignmentMap;
	int currentReadCount = 0;
	long currentBaseCount = 0;

	public String msa = "kalign";
	private String global = "hmm";

	public double scoreThreshold = 2;
	public boolean twoDOnly = false;

	public int numThead = 4;

	public RealtimeResistanceGene(int read, int time, String output, String resDB, String tmp) throws IOException{
		resistFinder = new ResistanceGeneFinder(this, output, resDB, tmp);
		resistFinder.setReadPeriod(read);
		resistFinder.setTimePeriod(time * 1000);
	}

	/**
	 * @param bamFile
	 * @param geneFile
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void typing(String bamFile) throws IOException, InterruptedException{
		//DateFormat df = new SimpleDateFormat("dd/MM/yy HH:mm:ss");
		//Logging.info("START : " + df.format(Calendar.getInstance().getTime()));		

		Logging.info("Resistance identification ready at " + new Date());
		alignmentMap = new HashMap<String, ArrayList<Sequence>> ();
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if ("-".equals(bamFile))
			samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			samReader = SamReaderFactory.makeDefault().open(new File(bamFile));

		SAMRecordIterator samIter = samReader.iterator();

		Thread t = new Thread(resistFinder, "SSS");
		t.start();

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
			//assert: the read sequence is stored in readSequence with the right direction
			String	geneID = record.getReferenceName();
			if (!resistFinder.geneMap.containsKey(geneID))
				continue;

			int refLength =  resistFinder.geneMap.get(geneID).length();

			int refStart = record.getAlignmentStart();
			int refEnd   = record.getAlignmentEnd();

			if(refStart > 99 || refEnd < refLength - 99)
				continue;			

			synchronized(this){
				ArrayList<Sequence> alignmentList = alignmentMap.get(geneID);
				if (alignmentList == null){
					alignmentList = new ArrayList<Sequence>();
					alignmentMap.put(geneID, alignmentList);				
				}
				//put the sequence into alignment list
				Sequence readSeq = HTSUtilities.readSequence(record, readSequence, 99, refLength-99);
				alignmentList.add(readSeq);
			}//synchronized(this)
		}//while
		resistFinder.stopWaiting();
		samIter.close();
		samReader.close();

		Logging.info("END : " + new Date());
	}	

	//TODO: way to improve performance:
	//1. 
	//3. options: gene or antibiotics class
	//4. 
	//5. Future improve: incrementally multiple alignment

	public static class ResistanceGeneFinder extends RealtimeAnalysis{
		public String prefix = "tmp";		
		RealtimeResistanceGene resistGene;

		HashMap<String, ArrayList<Sequence>> alignmentMapSnap = new HashMap<String, ArrayList<Sequence> >();
		HashMap<String, String> gene2GeneName;
		HashMap<String, String> gene2Group;			
		HashMap<String, Sequence> geneMap;
		ArrayList<String> geneList = new ArrayList<String>();
		//Set of genes confirmed to have found
		HashSet<String> predictedGenes = new HashSet<String>();

		SequenceOutputStream sos;

		public ResistanceGeneFinder(RealtimeResistanceGene resistGene, String output, String resDB, String tmp) throws IOException{
			this.resistGene = resistGene;			
			getGeneClassInformation(resDB + "/DB.fasta");

			File file = new File(resDB + "/geneList");
			if (file.exists()){
				readGeneInformation(resDB + "/geneList");
			}			

			Logging.info("geneList = " + geneList.size());
			Logging.info("geneMap = " + geneMap.size());
			Logging.info("gene2Group = " + gene2Group.size());
			Logging.info("gene2GeneName = " + gene2GeneName.size());

			sos = SequenceOutputStream.makeOutputStream(output);
			prefix = tmp;
		}	

		private void antiBioticAnalysis(){

			try {
				sos.print("##" + timeNow  + "\t" + (this.lastTime - this.startTime) + "\t" + this.lastReadNumber + "\n");

				//1. Make a snapshot of the current alignment
				synchronized(resistGene){
					//lastTime = System.currentTimeMillis();
					//lastReadNumber = resistGene.currentReadCount;
					for (String gene:resistGene.alignmentMap.keySet()){
						ArrayList<Sequence> readMap = resistGene.alignmentMap.get(gene);					 
						alignmentMapSnap.put(gene, (ArrayList<Sequence>) readMap.clone());
					}
				}//synchronized(resistGene)

				runIndex ++;
				//Now can make the call
				antiBioticsProfile();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
		/****
		 * 
		 * @throws IOException
		 * @throws InterruptedException
		 */

		int runIndex = 0;		
		private void antiBioticsProfile() throws IOException, InterruptedException{						 
			int jobNo = 0;			
			//Get list of genes from my
			ExecutorService executor = Executors.newFixedThreadPool((resistGene.numThead > 2)?resistGene.numThead-2:1);

			for (String geneID: geneList){				
				if (predictedGenes.contains(geneID))
					continue;

				ArrayList<Sequence> alignmentList =  alignmentMapSnap.get(geneID);
				Sequence consensus = 
					ErrorCorrection.consensusSequence(alignmentList, prefix + "_" + geneID + "_" + runIndex, resistGene.msa);

				if (consensus == null){
					continue;//gene
				}

				Sequence gene = geneMap.get(geneID);
				gene = gene.subSequence(99, gene.length() - 100);

				FSMThread thread = new FSMThread();		
				thread.resGeneFinder = this;
				thread.consensus =consensus;
				thread.gene = gene;
				thread.geneID = geneID;

				executor.execute(thread);
				jobNo ++;

				/*****************************************************************
				if (resistGene.global.equals("hmm")){
					double score = fsmAlignment(consensus, gene);
					Logging.info("SGF: " + score + " " + geneID + " " + alignmentList.size() + " " + (geneID) + " " + gene2Group.get(geneID));

					if (score >= resistGene.scoreThreshold){
						addPreditedGene(geneID);
						Logging.info("ADDF " + geneID);//
						//Logging.info("ADDF " + geneID + " " + resistGene.gene2Group.get(geneID)+ " " + resistGene.gene2PrimaryGroup.get(geneID) + " " + geneID);						
						continue;//for gene
					}					
				}else{
					/*****************************************************************
					String consensusFile = prefix + "consensus" + geneID + "_" + runIndex + ".fasta"; 
					consensus.writeFasta(consensusFile);				
					{	
						double score = checkNeedle(consensusFile, gene);
						Logging.info("SGF: " + score + " " + geneID + " " + alignmentList.size() + " " + (geneID) + " " + gene2Group.get(geneID));

						if (score >= resistGene.scoreThreshold){
							addPreditedGene(geneID);
							Logging.info("ADDF " + geneID);							
							continue;//for gene
						}
					}					
				}
				/*****************************************************************/
			}
			executor.shutdown(); 
			executor.awaitTermination(3, TimeUnit.DAYS);			
			Logging.info("===Found " + predictedGenes.size() + " vs " + geneMap.size() + "  " + alignmentMapSnap.size() + " with " + jobNo);
		}

		private void addPreditedGene(String geneID) throws IOException{
			predictedGenes.add(geneID);			
			sos.print(timeNow + "\t" + (this.lastTime - this.startTime)/1000 + "\t" + lastReadNumber + "\t" + resistGene.currentBaseCount + "\t" + geneID + "\t" + gene2GeneName.get(geneID) + "\t" + gene2Group.get(geneID) + "\n");			
			sos.flush();
		}


		private static double fsmAlignment(Sequence consensus, Sequence gene){
			//if (gene.length() > 2700 || consensus.length() > 4000 || gene.length() * consensus.length() > 6000000){
			//	Logging.info("SKIP " + gene.getName() + " " + gene.length() + " vs " + consensus.length());			
			//	return 0;
			//}
			//ProbThreeSM tsmF = new ProbThreeSM(gene);
			ProbOneSM tsmF = new ProbOneSM(gene);
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

		private double checkNeedle(String consensusFile, Sequence gene) throws IOException, InterruptedException{
			//Needle the gene
			String geneID = gene.getName();
			String faAFile = "geneAlleles/out_" + geneID + ".fasta";
			String needleOut = prefix + geneID + "_" + this.lastReadNumber + "_consensus.needle";

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

		private void addGeneInfo(HashMap<String, String> map, String key, String info){
			String s = map.get(key);
			if (s == null)
				map.put(key, info);
			else{
				if (!s.contains(info)){
					s = s + ", " + info;
					map.put(key, s);
				}
			}			
		}

		private void readGeneInformation(String geneInfoFile) throws IOException{
			BufferedReader bf = SequenceReader.openFile(geneInfoFile);
			String line = "";
			while((line = bf.readLine())!=null){
				String [] toks = line.trim().split(" ");
				if(toks.length <3)
					continue;

				addGeneInfo(gene2Group, toks[0], toks[2]);
				addGeneInfo(gene2GeneName, toks[0], toks[1]);

			}
			bf.close();
		}
		private void getGeneClassInformation(String geneFile) throws IOException{
			ArrayList<Sequence> drGeneList = FastaReader.readAll(geneFile, Alphabet.DNA());

			geneMap    = new HashMap<String, Sequence>();
			gene2Group = new HashMap<String, String>();
			gene2GeneName = new HashMap<String, String>();

			for (Sequence seq:drGeneList){
				geneMap.put(seq.getName(), seq);
				geneList.add(seq.getName());

				String desc = seq.getDesc();
				String [] toks = desc.split(";");
				for (String tok:toks){
					if (tok.startsWith("dg=")){
						addGeneInfo(gene2Group, seq.getName(), tok.substring(3));						
					}
					if (tok.startsWith("geneID=")){
						String proteinID = tok.substring(7);
						addGeneInfo(gene2GeneName, seq.getName(), proteinID);
					}
				}				
			}					
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#close()
		 */
		@Override
		protected void close() {
			try {
				sos.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#analysis()
		 */
		@Override
		protected void analysis() {
			antiBioticAnalysis();

		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#getCurrentRead()
		 */
		@Override
		protected int getCurrentRead() {
			return resistGene.currentReadCount;

		}
	}

	static class FSMThread implements Runnable{
		ResistanceGeneFinder resGeneFinder;
		Sequence consensus, gene;
		String geneID;

		/* (non-Javadoc)
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run() {
			double score = ResistanceGeneFinder.fsmAlignment(consensus, gene);
			Logging.info("SGF: " + score + " " + geneID + " " + " " + (geneID) + " " + resGeneFinder.gene2Group.get(geneID));

			if (score >= resGeneFinder.resistGene.scoreThreshold){
				synchronized(resGeneFinder){
					try {
						resGeneFinder.addPreditedGene(geneID);
						Logging.info("ADDF " + geneID);//
					} catch (IOException e) {						
						e.printStackTrace();
					}
				}
			}
		}

	}
}