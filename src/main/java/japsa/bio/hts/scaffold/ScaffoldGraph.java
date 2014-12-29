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
import japsa.bio.hts.scaffold.GapCloser.AlignmentRecord;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;


public class ScaffoldGraph{
	Scaffold [] scaffolds;
	int [] head;
	ArrayList<Contig> contigs;				
	int nScaffolds;
	public double estimatedCov = 0;
	double estimatedLength = 0;

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
		scaffolds = new Scaffold[contigs.size()];
		nScaffolds = contigs.size();
		head = new int[contigs.size()];			

		for (int i = 0; i < contigs.size();i++){				
			scaffolds[i] = new Scaffold(contigs.get(i));//TODO
			head[i] = i;//pointer contig -> scaffold
			//point to the head of the scaffold
		}//for
	}//constructor

	
	ArrayList<ContigBridge> bridgeList = new ArrayList<ContigBridge>();
	HashMap<String, ContigBridge> bridgeMap= new HashMap<String, ContigBridge>();	
	
	
	public void makeConnections(String bamFile, double cov, int threshold, int qual) throws IOException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));	

		SAMRecordIterator iter = reader.iterator();
		
		String readName = "";
		ArrayList<AlignmentRecord> samList = new ArrayList<AlignmentRecord>();

		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			if (rec.getReadUnmappedFlag())
				continue;
			if (rec.getMappingQuality() < qual)
				continue;

			AlignmentRecord myRec = new AlignmentRecord(rec, contigs.get(rec.getReferenceIndex()).length());
			//Logging.info(myRec.pos() + "#" + myRec.name + " " + myRec.refIndex);
			if (!myRec.useful)
				continue;

			if (rec.getReadName().equals(readName)) {
				for (AlignmentRecord s : samList) {
					this.addBridge(s, myRec, cov, threshold);
				}
			} else {
				samList.clear();				
			}
			samList.add(myRec);
			readName = rec.getReadName();
		}// while
		iter.close();

		//outOS.close();
		reader.close();
	}
	
	
	/*********************************************************************************/
	public void addBridge(AlignmentRecord a, AlignmentRecord b, double cov,double minScore){
		if (a.refIndex > b.refIndex){
			AlignmentRecord t = a;a=b;b=t;
		}

		String readName = a.name;
		double rate = 1.0 * (Math.abs(a.refEnd - a.refStart) + Math.abs(b.refEnd - b.refStart))
				/
				(Math.abs(a.readEnd - a.readStart) + Math.abs(b.readEnd - b.readStart));		
						
		//See if this is reliable
		double score = Math.min(a.refAlign, b.refAlign);
		int alignP = (int) ((b.readStart - a.readStart) * rate);
		int alignD = (a.strand == b.strand)?1:-1;
		

		int overlap = Math.min(					
				Math.max(a.readEnd, a.readStart) - Math.min(b.readStart,b.readEnd),
				Math.max(b.readEnd, b.readStart) - Math.min(a.readStart,a.readEnd));

		if (    (overlap > minScore/2)       
				||  contigs.get(a.refIndex).getCoverage() > cov * 1.4
				|| contigs.get(b.refIndex).getCoverage() > cov * 1.4
				|| contigs.get(a.refIndex).getCoverage() < cov / 1.6
				|| contigs.get(b.refIndex).getCoverage() < cov / 1.6
				|| a.refLength < 2 * minScore
				|| b.refLength < 2* minScore
				|| score < minScore 
				){				
			//System.out.println("IGNORE"
			//		+ " " + a.refIndex 
			//		+ " " + b.refIndex
			//		+ " " + readName
			//		+ " " + a.pos() + " " + b.pos()
			//		+ " " + score
			//		+ " " + (contigs.get(a.refIndex).getCoverage()/this.estimatedCov)
			//		+ " " + (contigs.get(b.refIndex).getCoverage()/this.estimatedCov)
			//		+ " " + alignP
			//		+ " " + (alignD==1)
			//		);

			return;
		}
		
		
		//if (readName.contains("channel_423_read_26_twodimentional")){
		//	int x = 1, y = 2;
		//	y += x;			
		//}
		
		int gP = (alignP + (a.strand ? a.refStart:-a.refStart) - (b.strand?b.refStart:-b.refStart));
		if (!a.strand)
			gP = -gP;		
		
		ScaffoldVector trans = new ScaffoldVector(gP, alignD);		
		
		int bridgeID = 0;
		ContigBridge bridge;
		while (true){
			String hash = ContigBridge.makeHash(a.refIndex, b.refIndex, bridgeID);
			bridge = bridgeMap.get(hash);
			if (bridge == null){
				bridge = new ContigBridge(contigs.get(a.refIndex), contigs.get(b.refIndex), bridgeID);
				bridge.addConnection(a, b, trans, score);
				//add				
				bridgeList.add(bridge);
				bridgeMap.put(hash, bridge);
				
				break;
			}else if (bridge.consistentWith(trans)){
				//add
				bridge.addConnection(a, b, trans, score);
				break;
			}else{
				bridgeID ++;
				//continue;
			}
		}
		
	}
	/*********************************************************************************/
	
	public void connectBridges(){
		Collections.sort(bridgeList);
		System.out.println(bridgeList.size());
		for (ContigBridge bridge:bridgeList){
			System.out.println("CONNECT " + bridge.hashKey + " " + bridge.getScore() + 
					" " + bridge.getConnections().size() + 
					" (" + bridge.getTransVector().toString() + 
					") " + bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig));
			bridge.display();			
			
			int currentF = bridge.secondContig.index;
			int currentT = bridge.firstContig.index;

			int headF = head[currentF];
			int headT = head[currentT];			
			
			
			ScaffoldVector trans = bridge.getTransVector();
			
			if (headF < headT){
				//swap
				trans = ScaffoldVector.reverse(trans);
				int tmp = currentF;
				currentF = currentT;
				currentT = tmp;

				headF = head[currentF];
				headT = head[currentT];
			}

			if (headT == headF){				
				ScaffoldVector v = ScaffoldVector.composition(trans,contigs.get(currentT).getVector());
				int t_dis = Math.abs(contigs.get(currentF).getRelPos() -  v.getMagnitute());
				int t_len = (scaffolds[headT].getEnd() - scaffolds[headT].getStart());

				System.out.println("NOT Connect " + currentF + " (" + headF +") and " + currentT + " (" + headT +") ==== " + contigs.get(currentF).getRelPos() + " - " + v.getMagnitute() + "(" + t_dis + ") vs " + (t_len) + "  " + (t_dis * 1.0/t_len) );

			}else{
				System.out.println("Before Connect " + currentF + " (" + headF +") and " + currentT + " (" + headT +") " + (scaffolds[headT].getEnd() - scaffolds[headT].getStart()) + " " + (scaffolds[headF].getEnd() - scaffolds[headF].getStart()) + " " + (scaffolds[headT].getEnd() - scaffolds[headT].getStart() + scaffolds[headF].getEnd() - scaffolds[headF].getStart()));
				ScaffoldVector revNew = ScaffoldVector.reverse(contigs.get(currentF).getVector());
				//rev = headF -> currentF
				
				for (Contig ctg:scaffolds[headF]){					
					ctg.composite(revNew);
					ctg.composite(trans);
					ctg.composite(contigs.get(currentT).getVector());
					scaffolds[headT].addContig(ctg);
					
					int contigID = ctg.getIndex();
					head[contigID] = headT;											

					int newS = ctg.getRelPos();
					int newE = ctg.getRelPos() + ctg.getRelDir() * ctg.length();
					if (newS > newE){
						int t=newS;newS = newE;newE = t;
					}

					if (newS < scaffolds[headT].getStart()){
						System.out.println("Extend " + headT + " start from " + scaffolds[headT].getStart() + " to " + newS);
						scaffolds[headT].setStart(newS); 
					}

					if (newE > scaffolds[headT].getEnd()){
						System.out.println("Extend " + headT + " end from " + scaffolds[headT].getEnd() + " to " + newE);
						scaffolds[headT].setEnd(newE); 
					}
				}
				System.out.println("After Connect " + currentF + " (" + headF +") and " + currentT + " (" + headT +") " + (scaffolds[headT].getEnd() - scaffolds[headT].getStart()));
				nScaffolds --;
				scaffolds[headT].viewNew();
			}			
		}
		
	}
	/*****************************************************************
	public void connect1(AlignmentRecord a, AlignmentRecord b, double cov,double minScore){
		
		String readName = a.name;
		double rate = (Math.abs(a.refEnd - a.refStart) + Math.abs(b.refEnd - b.refStart))
				/
				(Math.abs(a.readEnd - a.readStart) + Math.abs(b.readEnd - b.readStart));		
						
		//See if this is reliable
		double score = Math.min(a.refAlign, b.refAlign);

		int overlap = Math.min(					
				Math.max(a.readEnd, a.readStart) - Math.min(b.readStart,b.readEnd),
				Math.max(b.readEnd, b.readStart) - Math.min(a.readStart,a.readEnd));

		if ( overlap > minScore/2)
			return;
		//repeats
		if (       contigs.get(a.refIndex).getCoverage() > cov * 1.4
				|| contigs.get(b.refIndex).getCoverage() > cov * 1.4
				|| contigs.get(a.refIndex).getCoverage() < cov / 1.4
				|| contigs.get(b.refIndex).getCoverage() < cov / 1.4
				|| a.refLength < 2 * minScore
				|| b.refLength < 2* minScore
				|| score < minScore 
				){				
			System.out.println("IGNORE"
					+ " " + a.refIndex 
					+ " " + b.refIndex
					+ " " + readName
					+ " " + a.pos() + " " + b.pos()
					+ " " + score
					+ " " + contigs.get(a.refIndex).getCoverage()
					+ " " + contigs.get(b.refIndex).getCoverage()					
					);

			return;
		}		
		
		int currentF = b.refIndex;
		int currentT = a.refIndex;

		int headF = head[currentF];
		int headT = head[currentT];			
		
		int alignP = (int) ((b.readStart - a.readStart) * rate);
		int alignD = (a.strand == b.strand)?1:-1;
		
		int gP = (alignP + (a.strand ? a.refStart:-a.refStart) - (b.strand?b.refStart:-b.refStart));
		if (!a.strand)
			gP = -gP;		
		
//		System.out.println("XXX " + (alignD == dF) + " " + gP + " vs " + pF);
		
		ScaffoldVector trans = new ScaffoldVector(gP, alignD);
		
		if (headF < headT){
			//swap
			trans = ScaffoldVector.reverse(trans);
			int tmp = currentF;
			currentF = currentT;
			currentT = tmp;

			headF = head[currentF];
			headT = head[currentT];
		}

		if (headT == headF){				
			ScaffoldVector v = ScaffoldVector.composition(trans,contigs.get(currentT).getVector());
			int t_dis = Math.abs(contigs.get(currentF).getRelPos() -  v.getPos());
			int t_len = (scaffolds[headT].getEnd() - scaffolds[headT].getStart());

			System.out.println("NOT Connect " + currentF + " (" + headF +") and " + currentT + " (" + headT +") ==== " + contigs.get(currentF).getRelPos() + " - " + v.getPos() + "(" + t_dis + ") vs " + (t_len) + "  " + (t_dis * 1.0/t_len) );

		}else{
			System.out.println("Before Connect " + currentF + " (" + headF +") and " + currentT + " (" + headT +") " + (scaffolds[headT].getEnd() - scaffolds[headT].getStart()) + " " + (scaffolds[headF].getEnd() - scaffolds[headF].getStart()) + " " + (scaffolds[headT].getEnd() - scaffolds[headT].getStart() + scaffolds[headF].getEnd() - scaffolds[headF].getStart()));
			ScaffoldVector revNew = ScaffoldVector.reverse(contigs.get(currentF).getVector());
			//rev = headF -> currentF
			
			for (Contig ctg:scaffolds[headF]){					
				ctg.composite(revNew);
				ctg.composite(trans);
				ctg.composite(contigs.get(currentT).getVector());
				scaffolds[headT].addContig(ctg);
				
				int contigID = ctg.getIndex();
				head[contigID] = headT;											

				int newS = ctg.getRelPos();
				int newE = ctg.getRelPos() + ctg.getRelDir() * ctg.length();
				if (newS > newE){
					int t=newS;newS = newE;newE = t;
				}

				if (newS < scaffolds[headT].getStart()){
					System.out.println("Extend " + headT + " start from " + scaffolds[headT].getStart() + " to " + newS);
					scaffolds[headT].setStart(newS); 
				}

				if (newE > scaffolds[headT].getEnd()){
					System.out.println("Extend " + headT + " end from " + scaffolds[headT].getEnd() + " to " + newE);
					scaffolds[headT].setEnd(newE); 
				}
			}
			System.out.println("After Connect " + currentF + " (" + headF +") and " + currentT + " (" + headT +") " + (scaffolds[headT].getEnd() - scaffolds[headT].getStart()));
			nScaffolds --;
			scaffolds[headT].viewNew();
		}			
	}
	/*****************************************************************/
	public void viewStatus(){
		System.out.println(nScaffolds);
		for (int i = 0; i < scaffolds.length;i++){
			if (head[i] == i){
				System.out.println("Scaffold " + i + " length " + (scaffolds[i].getEnd() - scaffolds[i].getStart()));
				scaffolds[i].viewNew();
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