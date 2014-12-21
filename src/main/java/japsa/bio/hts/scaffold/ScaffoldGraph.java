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

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import work.GapCloser.AlignmentRecord;

import java.io.IOException;
import java.util.ArrayList;
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
	
	/*********************************************************************************/
	public void addBridge(AlignmentRecord a, AlignmentRecord b, double cov,double minScore){
		if (a.refIndex > b.refIndex){
			AlignmentRecord t = a;a=b;b=t;
		}
		
		double rate = 1.0 /0.95;
		String readName = a.name;		

		int sA = (int) (a.readStart * rate -  a.refLeft  * (a.strand?1:-1));
		int eA = (int) (a.readEnd * rate + a.refRight * (a.strand?1:-1)) - sA;

		int sB = (int) (b.readStart * rate - b.refLeft * (b.strand?1:-1)) - sA;
		int eB = (int) (b.readEnd * rate + b.refRight * (b.strand?1:-1)) - sA;

		sA = 0;
		if (eA < 0){
			eA = -eA;
			sB = -sB;
			eB = -eB;
		}
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
					+ " " + sA
					+ " " + eA
					+ " " + sB
					+ " " + eB
					+ " " + contigs.get(a.refIndex).getCoverage()
					+ " " + contigs.get(b.refIndex).getCoverage()					
					);
			return;
		}

		//int currentF = b.refIndex;
		//int currentT = a.refIndex;

		int pF = sB;
		int dF = (eB > sB) ? 1:-1;		
		ScaffoldVector trans = new ScaffoldVector(pF, dF);
		
		//if (currentF < currentT){
		//	//This is dead-code
		//	trans = ScaffoldVector.reverse(trans);
		//	int tmp = currentF;
		//	currentF = currentT;
		//	currentT = tmp;
		//}
		
		int bridgeID = 0;
		ContigBridge bridge;
		while (true){
			String hash = ContigBridge.makeHash(a.refIndex, b.refIndex, bridgeID);
			bridge = bridgeMap.get(hash);
			if (bridge == null){
				bridge = new ContigBridge(contigs.get(a.refIndex), contigs.get(b.refIndex), bridgeID);
				bridge.addConnection();
				
				
				//add				
				bridgeList.add(bridge);
				bridgeMap.put(hash, bridge);
				
				break;
			}else if (bridge.consistentWith(trans)){
				//add
				break;
			}else{
				bridgeID ++;
				//continue;
			}
		}
		
	}
	/*********************************************************************************/
	
	public void connect(AlignmentRecord a, AlignmentRecord b, double cov,double minScore){
		double rate = 1.0 /0.95;
		String readName = a.name;
		//if (!readName.contains("twod")){
		//	return;
		//}

		int sA = (int) (a.readStart * rate - a.refLeft * (a.strand?1:-1));
		int eA = (int) (a.readEnd * rate  + a.refRight * (a.strand?1:-1)) - sA;

		int sB = (int) (b.readStart * rate - b.refLeft * (b.strand?1:-1)) - sA;
		int eB = (int) (b.readEnd * rate  + b.refRight * (b.strand?1:-1)) - sA;

		sA = 0;
		if (eA < 0){
			eA = -eA;
			sB = -sB;
			eB = -eB;
		}
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
					+ " " + sA
					+ " " + eA
					+ " " + sB
					+ " " + eB
					+ " " + contigs.get(a.refIndex).getCoverage()
					+ " " + contigs.get(b.refIndex).getCoverage()					
					);

			return;
		}

		int currentF = b.refIndex;
		int currentT = a.refIndex;

		int headF = head[currentF];
		int headT = head[currentT];			

		int pF = sB;
		int dF = (eB > sB)?1:-1;
		
		ScaffoldVector trans = new ScaffoldVector(pF, dF);
		
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