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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;

import htsjdk.samtools.SAMRecord;
import japsa.bio.alignment.ProbFSM.Emission;
import japsa.bio.alignment.ProbFSM.ProbOneSM;
//import japsa.bio.alignment.ProbFSM.ProbThreeSM;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.tools.seq.CachedFastqWriter;
import japsa.tools.seq.CachedOutput;
import japsa.util.HTSUtilities;

/**
 * 
 * @author minhduc
 *
 */

public class RealtimeResistanceGene {
	private static final Logger LOG = LoggerFactory.getLogger(RealtimeResistanceGene.class);
	public static boolean writeAlignedOnly = false;
	public static boolean JSON = false;
	public static boolean OUTSEQ=true;
	private static HashMap<String, ArrayList<String>> allAlignedReads;
	public static boolean realtimeAnalysis = false;
	public static boolean runKAlign = false;
	public static boolean writeSep=true;
//	public final CachedOutput fqw ;
	Map<String, CachedOutput> fqw = new HashMap<String, CachedOutput>();
	//File fqDir = new File("fastqs");
	
	private ResistanceGeneFinder resistFinder;
	private HashMap<String, ArrayList<Sequence>> alignmentMap;
	private Integer currentReadCount = 0;
	private Long currentBaseCount = 0L;

	public String msa = "kalign";	

	public Double scoreThreshold = 2D;
	public Boolean twoDOnly = false;
	public Integer numThead = 4;

	public void getOutfiles(List<String> res){
		Iterator<CachedOutput> it =  this.fqw.values().iterator();
		while(it.hasNext()){
			CachedOutput nxt = it.next();
			nxt.getOutFile(res);
		}
	}
	
	public RealtimeResistanceGene(Integer read, Integer time,File outdir,  String outputFile, String resDB, String recordPrefix) throws IOException {
	//	File outF = new File(outputFile);
		
	this.outdir = outdir;
    File geneListFile = new File(resDB + "/geneList");
    File fastaFile = new File(resDB + "/DB.fasta");
    OutputStream outStream = SequenceOutputStream.makeOutputStream(outputFile);
    InputStream resDBInputStream = new FileInputStream(geneListFile);
    InputStream fastaInputStream = new FileInputStream(fastaFile);

		resistFinder = new ResistanceGeneFinder(this, outStream, resDBInputStream, fastaInputStream, recordPrefix);
		resistFinder.setReadPeriod(read);
		resistFinder.setTimePeriod(time * 1000);
	}

	
	
	public RealtimeResistanceGene(Integer read, Integer time, OutputStream outStream, 
			InputStream resDBInputStream, InputStream fastaInputStream, String recordPrefix) throws IOException {
	//    this.fqw = null;
		//fqw = new CachedFastqWriter(new File("."), "resistance.fq", true, writeAlignedOnly);	
this.outdir = new File("./");
		resistFinder = new ResistanceGeneFinder(this, outStream, resDBInputStream, fastaInputStream, recordPrefix);
	    resistFinder.setReadPeriod(read);
	    resistFinder.setTimePeriod(time * 1000);
	}

	public void setScoreThreshold(Double value) {
		this.scoreThreshold = value;
	}

	/**
   * Support legacy filename-based calls
	 * @param bamFile
	 * @throws IOException
	 * @throws InterruptedException
	
	public void typing(String bamFile) throws IOException, InterruptedException {
		InputStream bamInputStream;

		if ("-".equals(bamFile))
			bamInputStream = System.in;
		else
			bamInputStream = new FileInputStream(bamFile);

		typing(bamInputStream);
	} */

	
	public void typing(Iterator<SAMRecord> samIter) throws IOException, InterruptedException{
		LOG.info("Resistance identification ready at " + new Date());

		alignmentMap = new HashMap<String, ArrayList<Sequence>> ();
		
		if(OUTSEQ)
			allAlignedReads = new HashMap<String, ArrayList<String>> ();
		
		//SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
	//	SamReader samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamInputStream));
	//	SAMRecordIterator samIter = samReader.iterator();
		if(realtimeAnalysis){
		Thread t = new Thread(resistFinder, "SSS");
		t.start();
		}

		String readName = "";
		//A dummy sequence
		Sequence readSequence = new Sequence(Alphabet.DNA(),1,"");
		Map<String, List<SAMRecord>> records = new HashMap<String, List<SAMRecord>>();
		List<String>resist = new ArrayList<String>();
		while (samIter.hasNext()){
			SAMRecord record = samIter.next();

			if (this.twoDOnly && !record.getReadName().contains("twodim"))
				continue;

			if (!record.getReadName().equals(readName)) {
				writeAll(records, false, readName, resist);
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
			
			String	geneID = record.getReferenceName();
			int refLength =  resistFinder.geneMap.get(geneID).length();
			int refStart = record.getAlignmentStart();
			int refEnd   = record.getAlignmentEnd();
			String resclass = resistFinder.gene2Group.get(geneID);
			if (!resistFinder.geneMap.containsKey(geneID))
				continue;
			if(refStart > 99 || refEnd < refLength - 99)
				continue;			

			synchronized(this){
				if (runKAlign && alignmentMap.get(geneID) == null)
					alignmentMap.put(geneID, new ArrayList<Sequence>());

				//put the sequence into alignment list
				List<SAMRecord> recs = records.get(resclass);
				if(recs==null) records.put(resclass,recs = new ArrayList<SAMRecord>());
				if(writeSep) recs.add(record);
				Sequence readSeq = HTSUtilities.readSequence(record, readSequence, 99, refLength-99);
				if(runKAlign) alignmentMap.get(geneID).add(readSeq);
				
				if(OUTSEQ){
					resist.add(geneID);
					//System.err.println(resist.size());
					
				}
			}//synchronized(this)
		}//while
		this.writeAll(records, true,readName, resist);
	//	this.fqw.close();
//for(Iterator<CachedOutput> it = this.fqw.values().iterator(); it.hasNext();){
//it.next().close();
	//	}
		resistFinder.stopWaiting();
		if(!realtimeAnalysis){
			resistFinder.run();
		}
		//samIter.close();
		//samReader.close();

		LOG.info("END : " + new Date());
	}	
	final File outdir;
	
	/** also clears 
	 * @param resist 
	 * @param readName */
	private void writeAll(Map<String, List<SAMRecord>> records, boolean close, String readName, List<String> resist) {
		
		Collections.sort(resist);
		if(resist.size()>0){
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<resist.size(); i++){
				sb.append(resist.get(i));
				if(i<resist.size()-1) sb.append(";");
			}
			resist.clear();
			String sbs = sb.toString();
			ArrayList<String> li = allAlignedReads.get(sbs);
			if(li==null) allAlignedReads.put(sbs, li = new ArrayList<String>());
			li.add(readName);
		}
		
		Iterator<String> it = records.keySet().iterator();
		while(it.hasNext()){
			String resc = it.next();
			List<SAMRecord> rec = records.get(resc);
			if(rec.size()==0) continue;
			CachedOutput co = this.fqw.get(resc);
			if(co==null){
				co = new CachedFastqWriter(outdir, resc, false, writeAlignedOnly);	
				fqw.put(resc, (CachedFastqWriter) co);
			}
			
			for(int i=0; i<rec.size(); i++){
				co.write(rec.get(i), resc);
			}
			if(close) co.close();
			rec.clear();
		}
		
	}
	
	 static double fsmAlignment(Sequence consensus, Sequence gene){
			ProbOneSM tsmF = new ProbOneSM(gene);
			double cost = 100000000;						
			for (int c = 0; c < 10; c++){
				tsmF.resetCount();
				Emission retState = tsmF.alignGenerative(consensus);
				if (cost  <= retState.myCost)
					break;//for c

				cost = retState.myCost;
				int emitCount = tsmF.updateCount(retState);
				LOG.info("Iter " + c + " : " + emitCount + " states and " + cost + " bits " + consensus.length() + "bp " + consensus.getName() + " by " + gene.getName());
				tsmF.reEstimate();	
			}				
			return (consensus.length() * 2 - cost) / gene.length();
		}


	//TODO: way to improve performance:
	//1. 
	//3. options: gene or antibiotics class
	//4. 
	//5. Future improve: incrementally multiple alignment

	public  class ResistanceGeneFinder extends RealtimeAnalysis{
		private final Map<String, ArrayList<Sequence>> alignmentMapSnap = new HashMap<String, ArrayList<Sequence>>();
		private final Map<String, String> gene2GeneName = new HashMap<String,String>();
		private final Map<String, String> gene2Group = new HashMap<String,String>();
		private final Map<String, Sequence> geneMap = new HashMap<String,Sequence>();
		private final List<String> geneList = new ArrayList<String>();
		CachedOutput fwq;
		
		//Set of genes confirmed to have found
		private final HashSet<String> predictedGenes = new HashSet<String>();

		private final Gson gson = new GsonBuilder().serializeNulls().create();

		private String recordPrefix = "tmp";
		private RealtimeResistanceGene resistGene;
		private SequenceOutputStream sequenceOutputStream;

	    public ResistanceGeneFinder(RealtimeResistanceGene resistGene, String outputFile, String resDB, String recordPrefix) throws IOException {
	    	this(resistGene, new FileOutputStream(outputFile), resDB, recordPrefix);
	    }
	
	    public ResistanceGeneFinder(RealtimeResistanceGene resistGene, OutputStream outStream, String resDB, String recordPrefix) throws IOException {
	    	this(resistGene, outStream, new FileInputStream(new File(resDB+"/geneList")), new FileInputStream(new File(resDB+"/DB.fasta")), recordPrefix);
	    }
	
	    public ResistanceGeneFinder(RealtimeResistanceGene resistGene, OutputStream outStream, InputStream resDBInputStream, InputStream fastaInputStream, String recordPrefix) throws IOException {
			this.resistGene = resistGene;
			
			
			readGeneClassInformation(fastaInputStream); // step 1
			readGeneInformation(resDBInputStream);      // step 2
	
			LOG.info("geneList = " + geneList.size());
			LOG.info("geneMap = " + geneMap.size());
			LOG.info("gene2Group = " + gene2Group.size());
			LOG.info("gene2GeneName = " + gene2GeneName.size());
	
			this.recordPrefix = recordPrefix;
			this.sequenceOutputStream = new SequenceOutputStream(outStream);
		}

	    private void readGeneClassInformation(InputStream fastaInputStream) throws IOException{
	    	List<Sequence> drGeneList = new ArrayList<Sequence>();
	
	    	SequenceReader sequenceReader = FastaReader.getReader(fastaInputStream);
	
	    	while (true) {
	    		Sequence fastaSequence = sequenceReader.nextSequence(Alphabet.DNA());
	    		if (fastaSequence == null)
	    			break;
	    		//LOG.info("fastaSequence = "+fastaSequence);
	    		drGeneList.add(fastaSequence);
	    	}
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
	    private void readGeneInformation(InputStream geneInfoInputStream) throws IOException {
	    	BufferedReader bf = SequenceReader.openInputStream(geneInfoInputStream);
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
	    private void addGeneInfo(Map<String, String> map, String key, String info){
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

	    @SuppressWarnings("unchecked")
		private void antiBioticAnalysis(){

	    	try {
	    		if (!JSON) {
	    			sequenceOutputStream.print("##" + lastTime + "\t" + (this.lastTime - this.startTime) + "\t" + this.lastReadNumber + "\n");
	    		}
    			else {
		          JsonObject jo = new JsonObject();
		          jo.addProperty("timestamp", lastTime);
		          jo.addProperty("timeLast", this.lastTime);
		          jo.addProperty("timeStart", this.startTime);
		          jo.addProperty("timeWaited", (this.lastTime - this.startTime));
		          jo.addProperty("lastReadNumber", this.lastReadNumber);
		          sequenceOutputStream.print(gson.toJson(jo));
		          sequenceOutputStream.println();
    			}
	    		sequenceOutputStream.flush();

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
				e.printStackTrace();
			} catch (InterruptedException e) {
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

				List<Sequence> alignmentList = alignmentMapSnap.get(geneID);
				Sequence consensus = 
					ErrorCorrection.consensusSequence(alignmentList, recordPrefix + "_" + geneID + "_" + runIndex, resistGene.msa);

				if (consensus == null){
					continue;//gene
				}

				Sequence gene = geneMap.get(geneID);
				gene = gene.subSequence(99, gene.length() - 100);
				gene.setName(geneID);//for log only

				FSMThread thread = new FSMThread();		
				thread.resGeneFinder = this;
				thread.consensus =consensus;
				thread.gene = gene;
				thread.geneID = geneID;

				executor.execute(thread);
				jobNo ++;
			}
			executor.shutdown(); 
			executor.awaitTermination(3, TimeUnit.DAYS);			
			LOG.info("===Found " + predictedGenes.size() + " vs " + geneMap.size() + "  " + alignmentMapSnap.size() + " with " + jobNo);
		}

		private void addPreditedGene(String geneID) throws IOException {
			predictedGenes.add(geneID);
			if (!JSON) {
				sequenceOutputStream.print(lastTime + "\t" + (this.lastTime - this.startTime) / 1000 + "\t" + lastReadNumber + "\t" + resistGene.currentBaseCount + "\t" + geneID + "\t" + gene2GeneName.get(geneID) + "\t" + gene2Group.get(geneID) + "\n");
			}
			else {
		        JsonObject jo = new JsonObject();
		        jo.addProperty("timeStamp", lastTime);
		        jo.addProperty("timeLast", this.lastTime);
		        jo.addProperty("timeStart", this.startTime);
		        jo.addProperty("timeWaited", (this.lastTime - this.startTime));
		        jo.addProperty("lastReadNumber", this.lastReadNumber);
		        jo.addProperty("currentBaseCount", resistGene.currentBaseCount);
		        jo.addProperty("currentReadCount", resistGene.currentReadCount);
		        jo.addProperty("geneID", geneID);
		        jo.addProperty("geneName", gene2GeneName.get(geneID));
		        jo.addProperty("geneGroup", gene2Group.get(geneID));
		        sequenceOutputStream.print(gson.toJson(jo));
		        sequenceOutputStream.println();
			}
			sequenceOutputStream.flush();
		}

		


		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#close()
		 */
		@Override
		protected void close() {
			try {
				sequenceOutputStream.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			//print out
			if(OUTSEQ){
				try (BufferedWriter bw = new BufferedWriter(new FileWriter("genes2reads.map"))) {
					for(String sp:predictedGenes){
						ArrayList<String> readList = allAlignedReads.get(sp);
						if(readList.size()==0)
							continue;
						bw.write(">"+sp+"\n");
						for(String read:readList)
							bw.write(read+"\n");
					}			

				} catch (IOException e) {

					e.printStackTrace();

				}
			}
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#analysis()
		 */
		@Override
		protected void analysis() {
			//antiBioticAnalysis();

		}
		@Override
		protected void lastAnalysis() {
			if(RealtimeResistanceGene.runKAlign)		antiBioticAnalysis();
			else{
				try{
				Iterator<String> genes = allAlignedReads.keySet().iterator();
				
				while(genes.hasNext()){
					String gene = genes.next();
					List<String> reads = allAlignedReads.get(gene);
					String[] gnes = gene.split(";");
					StringBuffer nmes = new StringBuffer();
					StringBuffer groups = new StringBuffer();
					for(int k=0; k<gnes.length; k++){
						String xtra = k<gnes.length-1 ? ";" : "";
						nmes.append(gene2GeneName.get(gnes[k])+xtra);
						groups.append(gene2Group.get(gnes[k])+xtra);
					}
					this.sequenceOutputStream.print(gene+"\t"+nmes.toString()+"\t"+groups.toString()+"\t"+reads.size());//+"\t"+bases);
					this.sequenceOutputStream.println();
				}
				}catch(Exception exc){
					exc.printStackTrace();
				}
			}

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
			double score =fsmAlignment(consensus, gene);
			LOG.info("SGF: score=" + score + " geneID=" + geneID + " group=" + resGeneFinder.gene2Group.get(geneID));

			if (score >= resGeneFinder.resistGene.scoreThreshold){
				synchronized(resGeneFinder){
					try {
						resGeneFinder.addPreditedGene(geneID);
						LOG.info("ADDF " + geneID);//
					} catch (IOException e) {						
						e.printStackTrace();
					}
				}
			}
		}

	}

	public void close() {
	Iterator<CachedOutput> it = this.fqw.values().iterator();
	while(it.hasNext()){
		it.next().close();
	}
		
	}
}