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

import japsa.bio.alignment.ProbFSM;
import japsa.bio.alignment.ProbFSM.Emission;
import japsa.bio.alignment.ProbFSM.ProbThreeSM;
import japsa.bio.amra.MLSTyping;
import japsa.bio.amra.MLSTyping.MLSType;
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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * @author minhduc
 *
 */
public class RealtimeMLST{
    private static final Logger LOG = LoggerFactory.getLogger(RealtimeMLST.class);
	/////////////////////////////////////////////////////////////////////////////	
	private double minQual = 0;
	private boolean twoDOnly = false;
	//private int numThread = 16;

	ArrayList<Sequence> [] alignmentLists;

	ArrayList<Sequence> geneList;
	HashMap<String, Sequence> geneMap;	

	int currentReadCount = 0;
	long currentBaseCount = 0;

	RealtimeMLSTyper typer;	

	@SuppressWarnings("unchecked")
	public RealtimeMLST(String mlstDir, String output, int minRead, int minTime) throws IOException{
		typer = new RealtimeMLSTyper(this, mlstDir, output);		
		typer.setReadPeriod(minRead);
		typer.setTimePeriod(minTime * 1000);


		geneList = SequenceReader.readAll(mlstDir + "/bwaIndex/genes.fasta", Alphabet.DNA());		 
		geneMap = new HashMap<String, Sequence>();
		for (Sequence gene:geneList){
			geneMap.put(gene.getName(), gene);
		}		
		alignmentLists = new ArrayList[geneList.size()];
		for (int i = 0; i < geneList.size(); i++)
			alignmentLists[i] = new ArrayList<Sequence>();
	}


	/**
	 * @param minQual the minQual to set
	 */
	public void setMinQual(double minQual) {
		this.minQual = minQual;
	}

	/**
	 * @param twoDOnly the twoDOnly to set
	 */
	public void setTwoDOnly(boolean twoDOnly) {
		this.twoDOnly = twoDOnly;
	}
	/**
	 * @param bamFile
	 * @param top
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void typing(String bamFile, int top) throws IOException, InterruptedException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if ("-".equals(bamFile))
			samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			samReader = SamReaderFactory.makeDefault().open(new File(bamFile));
		SAMRecordIterator samIter = samReader.iterator();

		//ExecutorService executor = Executors.newFixedThreadPool((numThread>2) ? (numThread - 2) : 1);

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

			if (record.getMappingQuality() < minQual)
				continue;

			//assert: the read sequence is stored in readSequence with the right direction
			String	geneID = record.getReferenceName();
			if (!geneMap.containsKey(geneID))
				continue;

			int refLength =  geneMap.get(geneID).length();
			Sequence readSeq = HTSUtilities.spanningSequence(record, readSequence, refLength, 0);

			if (readSeq == null){
				LOG.warn("Read sequence is NULL sequence ");
			}else{
				//MLFSMThread mlFSM = new MLFSMThread(this.typer, record.getReferenceIndex(), readSeq);
				//executor.execute(mlFSM);

				synchronized(alignmentLists){
					alignmentLists[record.getReferenceIndex()].add(readSeq);					
				}
			}
		}//while
		samIter.close();
		samReader.close();

		//executor.shutdown(); 
		//executor.awaitTermination(3, TimeUnit.DAYS);
		typer.stopWaiting();//Tell typer to stop		
		//typer.makeMLSTTyping0(top);		
	}	

	public static class RealtimeMLSTyper extends RealtimeAnalysis{
		RealtimeMLST typing;
		public SequenceOutputStream countsOS;
		MLSTyping mlstScheme;
		int numGenes = 7;
		int top = 10;
		double [][] typerScoreMatrix; 		
		String prefix = "_tmp_" + System.currentTimeMillis() + "_";

		RealtimeMLSTyper(RealtimeMLST typing, String mlstDir, String output) throws IOException{
			this.typing = typing;
			countsOS = SequenceOutputStream.makeOutputStream(output);
			mlstScheme = new MLSTyping(mlstDir);			

			typerScoreMatrix = new double[numGenes][];
			for (int i = 0; i < numGenes; i++){				
				typerScoreMatrix[i] = new double[mlstScheme.alleles(i).size()];
				Arrays.fill(typerScoreMatrix[i], 0);
			}
		}


		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#close()
		 */
		@Override
		protected void close() {
			try{				
				countsOS.close();
			}catch (Exception e){
				e.printStackTrace();
			}
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#analysis()
		 * i is the source
		 */
		@Override
		protected void analysis(int i) {
			try {
				makeTypingConsensus();
				makeTypingMLwithFSM(top);
			} catch (Exception e) {				 
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


		private void makeTypingConsensus() throws IOException, InterruptedException{
			int count = typing.currentReadCount;

			ExecutorService executor = Executors.newFixedThreadPool(7);
			synchronized(typing.alignmentLists){
				for (int i = 0; i < typing.alignmentLists.length;i++){
					Sequence consensus = ErrorCorrection.consensusSequence(typing.alignmentLists[i], prefix + mlstScheme.getGeneName(i) + "kalign" + count, "kalign");
					if (consensus == null){
						LOG.warn("No read found for " + mlstScheme.getGeneName(i));
						continue;
					}

					MLFSMThread mlFSM = new MLFSMThread(this, i, consensus);
					executor.execute(mlFSM);
				}//for i
			}//synchronized
			executor.shutdown(); 
			executor.awaitTermination(3, TimeUnit.DAYS);
		}

		private void makeTypingMLwithFSM(int top) throws IOException{
			synchronized(typerScoreMatrix){
				for (MLSType type:mlstScheme.profiles){
					type.typeScore = 0;
					for (int i = 0; i < numGenes;i++){
						type.typeScore += typerScoreMatrix[i][mlstScheme.alleleIndex(type, i)];
					}
				}
			}

			Collections.sort(mlstScheme.profiles);
			Collections.reverse(mlstScheme.profiles);

			if (top > mlstScheme.profiles.size())
				top = mlstScheme.profiles.size();

			countsOS.print("============================================== " + typing.currentReadCount);
			countsOS.println();
			for (int i = 0; i < top;i++){
				MLSType profile = mlstScheme.profiles.get(i);
				countsOS.print(profile.getST() + "\t" + profile.getScore());
				countsOS.println();
			}
		}


		@Override
		protected int numSources() {
			// need to fix this to deal with multiple sources
			return 1;
			// TODO Auto-generated method stub
//			return 0;
		}
	}

	static class MLFSMThread implements Runnable{
		RealtimeMLSTyper typer;
		int geneIndex;
		Sequence readSeq;

		MLFSMThread(RealtimeMLSTyper typer, int index, Sequence read){
			this.typer = typer;
			this.geneIndex = index;
			this.readSeq = read;
		}

		/* (non-Javadoc)
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run() {
			LOG.info("Running thread " + geneIndex + " on " + readSeq.getName());
			ArrayList<Sequence> alleles = typer.mlstScheme.alleles(geneIndex);
			int numAlleles = alleles.size();
			double [] myScore = new double[numAlleles];
			for (int x = 0; x < numAlleles; x ++){
				Sequence seq = alleles.get(x);
				int alleleNo = Integer.parseInt(seq.getName().split("_")[1]);
				if (!typer.mlstScheme.useAlleleNo(geneIndex, alleleNo))
					continue;
				
				ProbFSM tsm = new ProbThreeSM(seq);
				double cost = 1000000;

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
					LOG.info("Iter " + c + " : " + emitCount + " states and " + cost + " bits " + readSeq.length() + "bp " + readSeq.getName() + " by " + seq.getName());
					tsm.reEstimate();	
				}//for (iteration)
				//LOG.info(" Saving " + (readSeq.length() * 2 - cost) + " on " + readSeq.getName() + " by " + seq.getName());
				//if (cost < readSeq.length() * 2){
				myScore[x] = (readSeq.length() * 2 - cost);
				LOG.info("Score for " + seq.getName() + " " + myScore[x]);
				//}else
				//	myScore[x] = 0;
			}//for
			synchronized (typer.typerScoreMatrix){
				for (int x = 0; x < numAlleles; x ++){
					typer.typerScoreMatrix[geneIndex][x] = myScore[x];
				}
			}
			LOG.info("Done thread " + geneIndex + " on " + readSeq.getName());
		}//run		
	}
}