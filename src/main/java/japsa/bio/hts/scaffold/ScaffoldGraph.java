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
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;

public class ScaffoldGraph{
	final static int maxRepeatLength=8000; //Koren S et al 2013
	//will need to fix this later -- MDC
	public static int marginThres = 1000;
	public static int minContigLength = 200;
	
	
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
	public void makeConnections(String bamFile, double minCov, int qual, SequenceOutputStream connectStr, SequenceOutputStream statStr) throws IOException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));	

		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		ReadFilling readFilling = null;
		ArrayList<AlignmentRecord> samList = null;// alignment record of the same read;		
		BitSet bitSet =null;
		int readScore = 0, readLength = 0;
		if(statStr != null)
			statStr.print("#readID\tlength\tcovered\tscore\n");
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			if (rec.getReadUnmappedFlag())
				continue;
			if (rec.getMappingQuality() < qual)
				continue;

			AlignmentRecord myRec = new AlignmentRecord(rec, contigs.get(rec.getReferenceIndex()));

			//////////////////////////////////////////////////////////////////
			// Just to save this
			if (connectStr != null && readID.equals(myRec.readID)) {								
				for (AlignmentRecord s : samList) {
					connectStr.print(myRec.contig.index + " " + s.contig.index + " " + readID + " " + myRec.useful + " " + s.useful + " " + myRec.pos() + " " + s.pos());
					connectStr.println();
					connectStr.print(s.contig.index + " " + myRec.contig.index + " " + readID + " " + s.useful + " " + myRec.useful + " " + s.pos() + " " + myRec.pos());
					connectStr.println();
				}
			}
			//////////////////////////////////////////////////////////////////
			// make bridge of contigs that align to the same (Nanopore) read. 
			// Note that SAM file MUST be sorted based on readID (samtools sort -n)
			
			//not the first occurrance				
			if (readID.equals(myRec.readID)) {				
				if (myRec.useful){				
					for (AlignmentRecord s : samList) {
						if (s.useful)
							this.addBridge(readFilling, s, myRec, minCov);
					}
				}
			} else {
				//samList.clear();			
				if(statStr != null && bitSet !=null)
					statStr.print(readID + "\t" + readLength + "\t" + bitSet.cardinality() + "\t" + readScore + "\n");
				bitSet = new BitSet(myRec.readLength);
				readScore = 0;
				readLength = 0;
				
				samList = new ArrayList<AlignmentRecord>();
				readID = myRec.readID;	
				readFilling = new ReadFilling(new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + readID), samList);	
			}
			bitSet.set(myRec.readAlignmentStart(), myRec.readAlignmentEnd());
			//statStr.print(readID + "\t" + myRec.readAlignmentStart() + "->" + myRec.readAlignmentEnd() + ": " + myRec.refStart + "->" + myRec.refEnd + "\n");

			readScore += myRec.score;
			readLength = myRec.readLength;
				
			samList.add(myRec);

		}// while
		iter.close();

		//outOS.close();
		reader.close();		

		Logging.info("Sort list of bridges");		
		Collections.sort(bridgeList);		
	}



	/*********************************************************************************/
	private void addBridge(ReadFilling readSequence, AlignmentRecord a, AlignmentRecord b, double minCov){
		if (a.contig.index > b.contig.index){
			AlignmentRecord t = a;a=b;b=t;
		}
		// Rate of aligned lengths: ref/read (illumina contig/nanopore read)
		double rate = 1.0 * (Math.abs(a.refEnd - a.refStart) + Math.abs(b.refEnd - b.refStart))
				/
				(Math.abs(a.readEnd - a.readStart) + Math.abs(b.readEnd - b.readStart));		

		//See if this is reliable
		double score = Math.min(a.score, b.score);
		int alignP = (int) ((b.readStart - a.readStart) * rate);
		int alignD = (a.strand == b.strand)?1:-1;

		//(rough) relative position from ref_b (contig of b) to ref_a (contig of a) in the assembled genome
		int gP = (alignP + (a.strand ? a.refStart:-a.refStart) - (b.strand?b.refStart:-b.refStart));
		if (!a.strand)
			gP = -gP;	
		// TODO: gP == contig length -> plasmid contig
		if (	a.contig.getIndex() == b.contig.getIndex() 
				&& alignD > 0
				&& (Math.abs(gP)*1.0 / a.contig.length()) < 1.1 
				&& (Math.abs(gP)*1.0 / a.contig.length()) > 0.9 
				&& a.readLength < 1.1* a.contig.length()
			)
		{
			System.out.printf("Potential CIRCULAR or TANDEM contig %s map to read %s(length=%d): (%d,%d)\n"
							, a.contig.getName(), a.readID, a.readLength, gP, alignD);
			a.contig.isCircular = true;							
		}		
		// overlap length on aligned read (<0 if not overlap)
		int overlap = Math.min(	a.readAlignmentEnd() - b.readAlignmentStart(), b.readAlignmentEnd() - a.readAlignmentStart());				

		if (	overlap > Math.min(	.5 * Math.min(a.readAlignmentEnd()-a.readAlignmentStart(), b.readAlignmentEnd()-b.readAlignmentStart()),
									minContigLength)      
				|| a.contig.getCoverage() < minCov	// filter out contigs with inappropriate cov
				|| b.contig.getCoverage() < minCov
				){		
			return;
		}

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
	// use override method in ScaffoldGraphDFS instead
	public void connectBridges(boolean mode){		
		if(!mode){
			System.err.println("Not applicable yet!");
			return;
		}
			
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

			//not sure if this is neccesary (yes, it is)
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
			
			// updating relative position of involved contigs based on the bridge of two main contigs. 
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
			// merging two dequeues
			scaffolds[headT].combineScaffold(scaffolds[headF], bridge, posT, posF);
			System.out.println("After Connect " + contigF.index + " (" + headF +") and " + contigT.index + " (" + headT +") " + (scaffolds[headT].getLast().rightMost() - scaffolds[headT].getFirst().leftMost()));
			nScaffolds --;
			scaffolds[headT].view();
		}
	}

	public void printSequences(SequenceOutputStream out) throws IOException{
		System.out.println(nScaffolds);
		for (int i = 0; i < scaffolds.length;i++){
			if ((head[i] == i 
				&& !isRepeat(scaffolds[i].element())
				&& scaffolds[i].element().length() > 1000)
				|| scaffolds[i].element().isCircular
				){
				System.out.println("Scaffold " + i + " length " + (scaffolds[i].getLast().rightMost() - scaffolds[i].getFirst().leftMost()));
				scaffolds[i].viewSequence(out);

			}
		}

/*		for (Contig contig:contigs){
			System.out.printf("Contig %s used %6.3f of  %6.3f (%6.3f) Left over %6.3f times or %6.3f%% \n",
					contig.getName(),
					contig.portionUsed,
					contig.coverage/this.estimatedCov,
					contig.coverage,
					contig.coverage/this.estimatedCov - contig.portionUsed,
					100 - contig.portionUsed*100.0/ (contig.coverage/this.estimatedCov)
					);

		}*/
		for (Contig contig:contigs){
				contig.display();
		}
				
	}	
	// To check if this contig is likely a repeat or a singleton. If FALSE: able to be used as a milestone.
	public boolean isRepeat(Contig ctg){
		if (ctg.length() > maxRepeatLength || ctg.coverage < 1.3 * estimatedCov) 
			return false;
		else if (ctg.coverage > 1.5 * estimatedCov || ctg.length() < 1000)
			return true;
		else{
			for(ContigBridge bridge:ctg.bridges){
				Contig other = bridge.firstContig.getIndex()==ctg.getIndex()?bridge.secondContig:bridge.firstContig;
				if(other.getIndex()==ctg.getIndex()) continue;
				int dist=bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig);
				if( dist<0 && dist>-ctg.length()*.25){
					if(other.length() > maxRepeatLength || other.getCoverage() < 1.5*estimatedCov)
						return true;
				}
			}
			//return false;
		}
		if(ctg.length() < 500 || ctg.coverage < .5 * estimatedCov) // maybe not repeat but crap 
			return true;
		else 
			return false;
	}


}