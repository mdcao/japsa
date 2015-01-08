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
 * 20/12/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.hts.scaffold;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.Logging;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class ScaffoldGraph{
	ScaffoldDeque [] scaffolds;
	int [] head;
	ArrayList<Contig> contigs;				
	int nScaffolds;
	public double estimatedCov = 0;
	double estimatedLength = 0;
	ArrayList<ContigBridge> bridgeList = new ArrayList<ContigBridge>();
	HashMap<String, ContigBridge> bridgeMap= new HashMap<String, ContigBridge>();


	public ScaffoldGraph(String sequenceFile) throws IOException{
		//1. read in contigs
		SequenceReader reader = SequenceReader.getReader(sequenceFile);
		Sequence seq;
		contigs = new ArrayList<Contig>(); 

		int index = 0;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			Contig ctg = new Contig(index, seq);		

			String name = seq.getName();
			double mycov = 100;
			String [] toks = name.split("_");
			for (int i = 0; i < toks.length - 1;i++){
				if ("cov".equals(toks[i])){
					mycov = Double.parseDouble(toks[i+1]);
					break;
				}
			}
			estimatedCov += mycov * seq.length();
			estimatedLength += seq.length();
			ctg.setCoverage(mycov);

			contigs.add(ctg);
			index ++;
		}
		reader.close();

		estimatedCov /= estimatedLength;
		System.out.println("Cov " + estimatedCov + " Length " + estimatedLength);

		//2. Initialise scaffold graph
		scaffolds = new ScaffoldDeque[contigs.size()];
		nScaffolds = contigs.size();
		head = new int[contigs.size()];			

		for (int i = 0; i < contigs.size();i++){				
			scaffolds[i] = new ScaffoldDeque(contigs.get(i));
			head[i] = i;//pointer contig -> scaffold
			//point to the head of the scaffold
		}//for
	}//constructor


	/**
	 * Make connections between any two uniquely (non-repeat) contigs
	 * 
	 * @param bamFile
	 * @param minCov
	 * @param maxCov
	 * @param threshold
	 * @param qual
	 * @throws IOException
	 */
	public void makeConnections(String bamFile, double minCov, double maxCov,  int threshold, int qual) throws IOException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));	

		SAMRecordIterator iter = reader.iterator();

		int readID = -1;
		ReadFilling readFilling = null;
		ArrayList<AlignmentRecord> samList = null;// = new ArrayList<AlignmentRecord>();		

		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			if (rec.getReadUnmappedFlag())
				continue;
			if (rec.getMappingQuality() < qual)
				continue;

			AlignmentRecord myRec = new AlignmentRecord(rec, contigs.get(rec.getReferenceIndex()));
			
			//not the first occurance				
			if (readID == myRec.readID) {
				if (myRec.useful){				
					for (AlignmentRecord s : samList) {
						if (s.useful)
							this.addBridge(readFilling, s, myRec, minCov, maxCov, threshold);
					}
				}
			} else {
				//samList.clear();
				samList = new ArrayList<AlignmentRecord>();
				readID = myRec.readID;	
				readFilling = new ReadFilling(new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + readID/3 + "_" + readID % 3), samList);				
			}
			samList.add(myRec);

		}// while
		iter.close();

		//outOS.close();
		reader.close();		

		Logging.info("Sort list of bridges");		
		Collections.sort(bridgeList);		
	}



	/*********************************************************************************/
	private void addBridge(ReadFilling readSequence, AlignmentRecord a, AlignmentRecord b, double minCov, double maxCov, double minScore){
		if (a.contig.index > b.contig.index){
			AlignmentRecord t = a;a=b;b=t;
		}

		double rate = 1.0 * (Math.abs(a.refEnd - a.refStart) + Math.abs(b.refEnd - b.refStart))
				/
				(Math.abs(a.readEnd - a.readStart) + Math.abs(b.readEnd - b.readStart));		

		//See if this is reliable
		double score = Math.min(a.score, b.score);
		int alignP = (int) ((b.readStart - a.readStart) * rate);
		int alignD = (a.strand == b.strand)?1:-1;


		int overlap = Math.min(					
				Math.max(a.readEnd, a.readStart) - Math.min(b.readStart,b.readEnd),
				Math.max(b.readEnd, b.readStart) - Math.min(a.readStart,a.readEnd));

		if (    (overlap > minScore/2)       
				|| a.contig.getCoverage() > maxCov
				|| b.contig.getCoverage() > maxCov
				|| a.contig.getCoverage() < minCov
				|| b.contig.getCoverage() < minCov
				|| a.contig.length() < 2 * minScore
				|| b.contig.length() < 2* minScore
				|| score < minScore 
				){		
			//System.out.println("IGNORE"
			//		+ " " + a.refIndex 
			//		+ " " + b.refIndex
			//		+ " " + readID
			//		+ " " + a.pos() + " " + b.pos()
			//		+ " " + score
			//		+ " " + (contigs.get(a.refIndex).getCoverage()/this.estimatedCov)
			//		+ " " + (contigs.get(b.refIndex).getCoverage()/this.estimatedCov)
			//		+ " " + alignP
			//		+ " " + (alignD==1)
			//		);
			return;
		}

		int gP = (alignP + (a.strand ? a.refStart:-a.refStart) - (b.strand?b.refStart:-b.refStart));
		if (!a.strand)
			gP = -gP;		

		ScaffoldVector trans = new ScaffoldVector(gP, alignD);		

		int bridgeID = 0;
		ContigBridge bridge;
		while (true){
			String hash = ContigBridge.makeHash(a.contig.index, b.contig.index, bridgeID);
			bridge = bridgeMap.get(hash);
			if (bridge == null){
				bridge = new ContigBridge(a.contig, b.contig, bridgeID);
				bridge.addConnection(readSequence, a, b, trans, score);
				//add				
				bridgeList.add(bridge);
				bridgeMap.put(hash, bridge);

				break;
			}else if (bridge.consistentWith(trans)){
				//add
				bridge.addConnection(readSequence, a, b, trans, score);
				break;
			}else{
				bridgeID ++;
				//continue;
			}
		}

	}
	/*********************************************************************************/


	public void connectBridges(){		
		System.out.println(bridgeList.size());
		for (ContigBridge bridge:bridgeList){
			System.out.println("CONNECT " + bridge.hashKey + " " + bridge.getScore() + 
					" " + bridge.getConnections().size() + 
					" (" + bridge.getTransVector().toString() + 
					") " + bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig));
			bridge.display();

			Contig contigF = bridge.secondContig;
			Contig contigT = bridge.firstContig;

			int headF = head[contigF.index];
			int headT = head[contigT.index];			

			ScaffoldVector trans = bridge.getTransVector();

			//not sure if this is necccesary (yes, it is)
			if (headF < headT){
				//swap
				trans = ScaffoldVector.reverse(trans);
				Contig tmp = contigF;
				contigF = contigT;
				contigT = tmp;

				headF = head[contigF.index];
				headT = head[contigT.index];
			}

			if (headT == headF){				
				ScaffoldVector v = ScaffoldVector.composition(trans,contigT.getVector());
				int t_dis = Math.abs(contigF.getRelPos() -  v.getMagnitute());
				int t_len = (scaffolds[headT].getLast().rightMost() - scaffolds[headT].getFirst().leftMost());

				System.out.println("NOT Connect " + contigF.index + " (" + headF +") and " + contigT.index + " (" + headT +") ==== " + contigF.getRelPos() + " - " + v.getMagnitute() + "(" + t_dis + ") vs " + (t_len) + "  " + (t_dis * 1.0/t_len) );
				continue;
				//TODO: This is to close the circular chromosome/plasmid				
			}


			int posF = scaffolds[headF].isEnd(contigF);
			int posT = scaffolds[headT].isEnd(contigT);

			if (posT == 0 ){
				System.out.println("Opps " + contigF.index + " vs " + contigT.index);
				continue;
			}
			//assert: posT != 0, but note that posF may be equal to 0
			//checking if 



			System.out.println("Before Connect " + contigF.index + " (" + headF +") and " + contigT.index 
					+ " (" + headT +") " 
					+ (scaffolds[headT].getLast().rightMost() - scaffolds[headT].getFirst().leftMost()) 
					+ " " + (scaffolds[headF].getLast().rightMost() - scaffolds[headF].getFirst().leftMost()) 
					+ " " + (scaffolds[headT].getLast().rightMost() - scaffolds[headT].getFirst().leftMost() + scaffolds[headF].getLast().rightMost() - scaffolds[headF].getFirst().leftMost()));
			ScaffoldVector rev = ScaffoldVector.reverse(contigF.getVector());
			//rev = headF -> currentF				

			for (Contig ctg:scaffolds[headF]){					
				ctg.composite(rev);
				ctg.composite(trans);
				ctg.composite(contigT.getVector());
				//scaffolds[headT].addContig(ctg);

				int contigID = ctg.getIndex();
				head[contigID] = headT;											

				int newS = ctg.getRelPos();
				int newE = ctg.getRelPos() + ctg.getRelDir() * ctg.length();
				if (newS > newE){
					int t=newS;newS = newE;newE = t;
				}
			}
			
			scaffolds[headT].combineScaffold(scaffolds[headF], bridge, posT, posF);
			System.out.println("After Connect " + contigF.index + " (" + headF +") and " + contigT.index + " (" + headT +") " + (scaffolds[headT].getLast().rightMost() - scaffolds[headT].getFirst().leftMost()));
			nScaffolds --;
			scaffolds[headT].view();
		}
	}

	public void viewStatus(){
		System.out.println(nScaffolds);
		for (int i = 0; i < scaffolds.length;i++){
			if (head[i] == i){
				System.out.println("Scaffold " + i + " length " + (scaffolds[i].getLast().rightMost() - scaffolds[i].getFirst().leftMost()));
				scaffolds[i].viewSequence();
			}
		}
	}

	public void printScaffoldSequence(SequenceOutputStream out) throws IOException{			
		for (int i = 0; i < scaffolds.length;i++){
			if (head[i] == i){
				//TODO
				//Sequence s = scaffolds[i].scaffoldSequence(vectors, seqs);
				//s.setName("scaffold" + i);
				//s.writeFasta(out);
			}
		}
	}
}