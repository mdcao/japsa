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
import japsa.bio.bac.MLSTyping;
import japsa.bio.bac.MLSTyping.MLSType;
import japsa.seq.Alphabet;
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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

/**
 * @author minhduc
 *
 */
public class RealtimeMLST{	
	/////////////////////////////////////////////////////////////////////////////	
	private double minQual = 0;
	private boolean twoDOnly = false;

	ArrayList<Sequence> geneList;
	HashMap<String, Sequence> geneMap;
	HashMap<String, ArrayList<Sequence>> alignmentMap;

	int currentReadCount = 0;
	long currentBaseCount = 0;

	public String prefix = "tmp";	
	public String msa = "kalign";	

	public int readNumber = 100;
	SequenceOutputStream datOS = null;



	MLSTyping mlstScheme;
	int numGenes = 7;
	public RealtimeMLST(String mlstDir) throws IOException{
		mlstScheme = new MLSTyping(mlstDir);
		geneList = SequenceReader.readAll(mlstDir + "/bwaIndex/genes.fas", Alphabet.DNA());		 
		geneMap = new HashMap<String, Sequence>();

		for (Sequence gene:geneList){
			geneMap.put(gene.getName(), gene);		
		}		
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



	private void makeMLSTTyping(int top) throws IOException, InterruptedException{
		makeMLSTTyping0(top);
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


	private void makeMLSTTyping0(int top) throws IOException, InterruptedException{		
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

		double [][] scoreMatrix = new double[numGenes][];
		for (int i = 0; i < numGenes; i++){
			String geneID = geneList.get(i).getName();
			scoreMatrix[i] = new double[mlstScheme.geneSeqs[i].size()];
			Arrays.fill(scoreMatrix[i], 0);
			Sequence consensus = consensusMap.get(geneID);
			if (consensus == null)
				continue;

			/************************** HMM ML *******************************/			
			for (int x = 0; x < mlstScheme.geneSeqs[i].size(); x ++){
				Sequence seq = mlstScheme.geneSeqs[i].get(x);
				ArrayList<Sequence> alignmentList =  alignmentMap.get(geneID);
				if (alignmentList != null){
					double thisScore = 0;
					for (Sequence readSeq:alignmentList){
						ProbOneSM tsm = new ProbOneSM(seq);
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
							Logging.info("Iter " + c + " : " + emitCount + " states and " + cost + " bits " + readSeq.length() + "bp " + readSeq.getName() + " by " + seq.getName());
							tsm.reEstimate();	
						}//for (iteration)
						Logging.info(" Saving " + (readSeq.length() * 2 - cost) + " on " + readSeq.getName() + " by " + seq.getName());						
						if (cost < readSeq.length() * 2){
							thisScore += (readSeq.length() * 2 - cost);
						}						
					}//for read sequence
					scoreMatrix[i][x] = thisScore;					
				}//if
			}
			/****************************************************************/

			/************************ Needle *******************************
			for (int x = 0; x < mlstScheme.geneSeqs[i].size(); x ++){
				Sequence seq = mlstScheme.geneSeqs[i].get(x);
				scoreMatrix[i][x] = ErrorCorrection.needle(consensus, seq, prefix + "_" + geneID + "_" + this.currentReadCount);
			}			
			/****************************************************************/

			/************************** HMM ********************************			
			for (int x = 0; x < mlstScheme.geneSeqs[i].size(); x ++){
				Sequence seq = mlstScheme.geneSeqs[i].get(x);
				ProbOneSM tsm = new ProbOneSM(seq);
				double cost = 100000000;						
				for (int c = 0; c < 10; c++){
					tsm.resetCount();
					Emission retState = tsm.alignGenerative(consensus);
					if (cost  <= retState.myCost)
						break;

					cost = retState.myCost;
					int emitCount = tsm.updateCount(retState);
					Logging.info("Iter " + c + " : " + emitCount + " states and " + cost + " bits " + consensus.length() + "bp " + consensus.getName() + " by " + seq.getName());
					tsm.reEstimate();
				}				 
				scoreMatrix[i][x] = consensus.length() * 2 - cost;
			}
			/****************************************************************/
		}
		///////////////////////////////////////////////////////////

		for (MLSType type:mlstScheme.profiles){
			type.typeScore = 0;
			for (int i = 0; i < numGenes;i++){
				type.typeScore += scoreMatrix[i][type.getAllele(i) - 1];
			}
		}
		Collections.sort(mlstScheme.profiles);
		Collections.reverse(mlstScheme.profiles);

		if (top > mlstScheme.profiles.size())
			top = mlstScheme.profiles.size();

		System.out.println("============================================== " + currentReadCount);
		for (int i = 0; i < top;i++){
			MLSType profile = mlstScheme.profiles.get(i);
			System.out.println(profile.getST() + "\t" + profile.getScore());						
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


				if (currentReadCount % readNumber == 0){
					makeMLSTTyping(top);
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

			if (record.getMappingQuality() < minQual)
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
}