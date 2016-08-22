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
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.JapsaFeature;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.Logging;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class ScaffoldGraph{
	public static int maxRepeatLength=7500; //for ribosomal repeat cluster in bacteria (Koren S et al 2013), it's 9.1kb for yeast.
	public static int marginThres = 1000;
	public static int minContigLength = 300;
	public static int minSupportReads = 1;
	public static boolean verbose = false;
	public boolean annotation = false;
	
	public String prefix = "out";
	
	Scaffold [] scaffolds;
	ArrayList<Contig> contigs;				
	public static double estimatedCov = 0;
	private static double estimatedLength = 0;
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
			double mycov = 1.0;
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
		if(verbose) 
			System.out.println("Cov " + estimatedCov + " Length " + estimatedLength);

		//2. Initialise scaffold graph
		scaffolds = new Scaffold[contigs.size()];		

		for (int i = 0; i < contigs.size();i++){				
			scaffolds[i] = new Scaffold(contigs.get(i));
			//point to the head of the scaffold
			contigs.get(i).head = i;
		}//for
		
		
	}//constructor
	
	public Contig getContig(String name){
		Contig res = null;
		for(Contig ctg:contigs)
			if(ctg.getName().contains(name)){
				res = ctg;
				break;
			}
		return res;
	}
	
	
	public synchronized int getN50(){
		int [] lengths = new int[scaffolds.length];
		int count=0;
		double sum = 0;
		for (int i = 0; i < scaffolds.length;i++){
			if(scaffolds[i].isEmpty()) continue;
			int len = scaffolds[i].length();
			if ((contigs.get(i).head == i 
					&& !isRepeat(contigs.get(i))
					&& len > maxRepeatLength
					)
					|| scaffolds[i].closeBridge != null)
			{
				lengths[count] = len;
				sum+=len;
				count++;
			}
		}
			
		Arrays.sort(lengths);

		int index = lengths.length;
		double contains = 0;
		while (contains < sum/2){
			index --;
			contains += lengths[index];
		}

		return lengths[index];	
	}


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
	public void makeConnections(String bamFile, double minCov, int qual) throws IOException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		
		SamReader reader;
		if ("-".equals(bamFile))
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			reader = SamReaderFactory.makeDefault().open(new File(bamFile));	

		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		ReadFilling readFilling = null;
		ArrayList<AlignmentRecord> samList = null;// alignment record of the same read;	

		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			if (rec.getReadUnmappedFlag())
				continue;
			if (rec.getMappingQuality() < qual)
				continue;

			AlignmentRecord myRec = new AlignmentRecord(rec, contigs.get(rec.getReferenceIndex()));


			//////////////////////////////////////////////////////////////////
			// make bridge of contigs that align to the same (Nanopore) read. 
			// Note that SAM file MUST be sorted based on readID (samtools sort -n)
			
			//not the first occurrance				
			if (readID.equals(myRec.readID)) {				
				if (myRec.useful){				
					for (AlignmentRecord s : samList) {
						if (s.useful){
							this.addBridge(readFilling, s, myRec, minCov);
							//...update with synchronized
							
						}
					}
				}
			} else {
				
				samList = new ArrayList<AlignmentRecord>();
				readID = myRec.readID;	
				readFilling = new ReadFilling(new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + readID), samList);	
			}			
			samList.add(myRec);

		}// while
		iter.close();

		//outOS.close();
		reader.close();		

		Logging.info("Sort list of bridges");		
		//Collections.sort(bridgeList);		
	}



	/*********************************************************************************/
	protected void addBridge(ReadFilling readSequence, AlignmentRecord a, AlignmentRecord b, double minCov){
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
			if(verbose) 
				System.out.printf("Potential CIRCULAR or TANDEM contig %s map to read %s(length=%d): (%d,%d)\n"
							, a.contig.getName(), a.readID, a.readLength, gP, alignD);
			a.contig.cirProb ++;							
		}		
		else{
			a.contig.cirProb--;
			b.contig.cirProb--;
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

		int count = 0;
		ContigBridge bridge, bridge_rev;
		while (true){
			int brgID = count, revID = count;
			if(a.contig.getIndex()==b.contig.getIndex()){
				brgID = 2*count;
				revID = brgID+1;
			}
			String 	hash = ContigBridge.makeHash(a.contig.index, b.contig.index, brgID),
					hash_rev = ContigBridge.makeHash(b.contig.index, a.contig.index, revID);

			bridge = bridgeMap.get(hash);
			bridge_rev = bridgeMap.get(hash_rev);
			if (bridge == null){
				assert bridge_rev==null:hash_rev + " not null!";
				bridge = new ContigBridge(a.contig, b.contig, brgID);
				bridge_rev = new ContigBridge(b.contig, a.contig, revID);
				
				bridge.addConnection(readSequence, a, b, trans, score);
				bridge_rev.addConnection(readSequence, b, a, ScaffoldVector.reverse(trans), score);

				a.contig.bridges.add(bridge);
				b.contig.bridges.add(bridge_rev);
				
				bridgeMap.put(hash, bridge);
				bridgeMap.put(hash_rev, bridge_rev);
	
				break;
			}
			if ((a.contig.getIndex() != b.contig.getIndex()) && bridge.consistentWith(trans)){
				assert bridge_rev!=null:hash_rev + "is null!";
				bridge.addConnection(readSequence, a, b, trans, score);
				bridge_rev.addConnection(readSequence, b, a, ScaffoldVector.reverse(trans), score);
				break;
			}
			if(a.contig.getIndex() == b.contig.getIndex()){
				assert bridge_rev!=null:hash_rev + "is null";
				if(bridge.consistentWith(trans)){
					bridge.addConnection(readSequence, a, b, trans, score);
					bridge_rev.addConnection(readSequence, b, a, ScaffoldVector.reverse(trans), score);
					break;
				}
				if(bridge.consistentWith(ScaffoldVector.reverse(trans))){
					bridge_rev.addConnection(readSequence, b, a, trans, score);
					bridge.addConnection(readSequence, b, a, ScaffoldVector.reverse(trans), score);
					break;
				}		
			}
			count ++;
		}//while

	}
	/**********************************************************************************************/
	public ContigBridge getReversedBridge(ContigBridge bridge){
		String hash = ContigBridge.makeHash(bridge.secondContig.index, bridge.firstContig.index, bridge.orderIndex);
		return bridgeMap.get(hash);
	}
	/**********************************************************************************************/
	// check if it's possible to extend from *contig* with *bridge* to another extended-already contig (contigF)
	// use for markers and unique bridge only
	protected int extendDirection(Contig contig, ContigBridge bridge){
		Contig contigF = bridge.secondContig;
		ScaffoldVector trans = bridge.getTransVector(); //contigF -> contig
		int pointer = Integer.signum(trans.magnitude * trans.direction); //pointer < 0 => tail of contigF on bridge
		assert scaffolds[contigF.head].size() > 1 : contigF.head;
		
		int headF = contigF.head;
		int direction = 0; //direction of extension on scaffoldT (we need to return direction on scaffoldF)
		ScaffoldVector headT2contigF = ScaffoldVector.composition(trans, contig.getVector());
		int rEnd = contig.rightMost(), rEndF = contigF.rightMost(headT2contigF),
			lEnd =  contig.leftMost(), lEndF = contigF.leftMost(headT2contigF);
		if(rEndF > rEnd){
			direction = 1;
		}
		else if(lEndF < lEnd){
			direction = -1;
		}
		else 
			return 0;
		if(verbose)
			System.out.println("Examining extending direction from contig " + contig.getIndex() + " to " + bridge.hashKey);
		Scaffold scaffoldF = scaffolds[headF];
		Contig 	prevMarker = scaffoldF.nearestMarker(contigF, false), // previous marker of contigF on corresponding scaffold
				nextMarker = scaffoldF.nearestMarker(contigF, true); // next marker of contigF on corresponding scaffold		
		
		ScaffoldVector rev = ScaffoldVector.reverse(contigF.getVector()); //rev = contigF->headF	
		if(prevMarker != null){
			ScaffoldVector toPrev = ScaffoldVector.composition(prevMarker.getVector(),rev); //contigF->prevMarker
			if(scaffoldF.indexOf(prevMarker) > scaffoldF.indexOf(contigF) && scaffoldF.closeBridge != null)
				toPrev = ScaffoldVector.composition(ScaffoldVector.reverse(scaffoldF.circle), toPrev);
			ScaffoldVector headT2Prev = ScaffoldVector.composition(toPrev, headT2contigF);
			int rEndPrev = prevMarker.rightMost(headT2Prev),
				lEndPrev = prevMarker.leftMost(headT2Prev);
			if(verbose){
				System.out.printf("Extending from contigT %d to targeted contig (contigF) %d with previous contig (prevMarker) %d \n", contig.getIndex(), contigF.getIndex(), prevMarker.getIndex());
				System.out.println("...headT->contig, contigF and prevMarker: " + contig.getVector() + headT2contigF + headT2Prev);
			}
			if ((direction > 0?rEndPrev > rEndF: lEndPrev < lEndF)){
				if((pointer<0?contigF.nextScore:contigF.prevScore) < bridge.getScore()){
					if(verbose)
						System.out.printf("=> go from %d to %d to %d \n", contig.getIndex(), contigF.getIndex(), prevMarker.getIndex());
					return -1;
				}
				else{
					if(verbose)
						System.out.printf("Bridge score not strong enough: %.2f < %.2f (%.2f)\n", 
										bridge.getScore(), pointer<0?contigF.nextScore:contigF.prevScore,
												pointer<0?contigF.prevScore:contigF.nextScore);
				
					return 0;
				}
			}
		}
		if(nextMarker != null){
			ScaffoldVector toNext = ScaffoldVector.composition(nextMarker.getVector(),rev); //contigF->nextMarker
			if(scaffoldF.indexOf(nextMarker) < scaffoldF.indexOf(contigF) && scaffoldF.closeBridge != null)
				toNext = ScaffoldVector.composition(scaffoldF.circle, toNext);
			ScaffoldVector headT2Next = ScaffoldVector.composition(toNext, headT2contigF);
			int rEndNext = nextMarker.rightMost(headT2Next),
				lEndNext = nextMarker.leftMost(headT2Next);
			if(verbose){
				System.out.printf("Extending from contigT %d to targeted contig (contigF) %d with next contig (nextMarker) %d \n", contig.getIndex(), contigF.getIndex(), nextMarker.getIndex());
				System.out.println("...headT->contig, contigF and nextMarker: " + contig.getVector() + headT2contigF + headT2Next);
			}
			
			if ((direction > 0? rEndNext > rEndF : lEndNext < lEndF)){
				//if((rev.direction<0?contigF.nextScore:contigF.prevScore) < bridge.getScore()){
				if((pointer<0?contigF.nextScore:contigF.prevScore) < bridge.getScore()){	
					if(verbose)
						System.out.printf("=> go from %d to %d to %d \n", contig.getIndex(), contigF.getIndex(), nextMarker.getIndex());
					return 1;
				}
				else{
					if(verbose)
						System.out.printf("Bridge score not strong enough: %.2f < %.2f (%.2f)\n", 
										bridge.getScore(), pointer<0?contigF.nextScore:contigF.prevScore,
												pointer<0?contigF.prevScore:contigF.nextScore);
					return 0;
					
				}
			}else{
				if(verbose)
					System.out.printf("Direction conflict: %d, %d %d or %d %d \n", direction, rEndNext, rEndF, lEndNext, lEndF);
			}
		}
		return 0;
	}
	/*********************************************************************************/
	public synchronized boolean joinScaffold(Contig contig, ContigBridge bridge, boolean firstDir, int secondDir){		
		if(verbose) {
			System.out.println("PROCEED TO CONNECT " + bridge.hashKey + " with score " + bridge.getScore() + 
					", size " + bridge.getConnections().size() + 
					", vector (" + bridge.getTransVector().toString() + 
					"), distance " + bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig));
			bridge.display();
		}

		
		Contig contigF = bridge.secondContig, contigT = contig;
		ScaffoldVector trans = bridge.getTransVector();
					
		int headF = contigF.head,
			headT = contigT.head;
		Scaffold 	scaffoldF = scaffolds[headF],
					scaffoldT = scaffolds[headT];
		int	posT = scaffoldT.isEnd(contigT);
		if (posT == 0){
			if(verbose) 
				System.out.println("Impossible to jump from the middle of a scaffold " + headT + ": contig " + contigT.index);
			return false;
		}

		if(verbose) 
			System.out.println("Before joining " + contigF.index + " (" + headF +") to " + contigT.index 
				+ " (" + headT +") " 
				+ (scaffoldT.getLast().rightMost() - scaffoldT.getFirst().leftMost()) 
				+ " " + (scaffoldF.getLast().rightMost() - scaffoldF.getFirst().leftMost()) 
				+ " " + (scaffoldT.getLast().rightMost() - scaffoldT.getFirst().leftMost() + scaffoldF.getLast().rightMost() - scaffoldF.getFirst().leftMost()));
		//===================================================================================================
		int index = scaffoldF.indexOf(contigF),
			count = index;

		ScaffoldVector rev = ScaffoldVector.reverse(contigF.getVector()); //rev = contigF->headF	

		int addScf=-1; 

		if(secondDir == -1){
			if(headF==headT){
				//if(posT!=1)
				if(firstDir)	
					return false;
				else{
					Contig nextMarker = scaffoldF.nearestMarker(contigF, true);
					if(nextMarker!=null){
						Contig ctg = scaffoldF.remove(index+1);
						Scaffold newScf = new Scaffold(ctg);
						ContigBridge brg = scaffoldF.bridges.remove(index);
						while(true){
							if(scaffoldF.size()==index+1) break;
							ctg= scaffoldF.remove(index+1);
							brg = scaffoldF.bridges.remove(index);
							newScf.addRear(ctg,brg);
						}
						newScf.trim();
						changeHead(newScf, nextMarker);
						addScf=nextMarker.getIndex();
					}
					scaffoldF.setCloseBridge(getReversedBridge(bridge));
					changeHead(scaffoldF, contigF);
					//TODO shoud we separate scaffolds in scaffolds[] and others to avoid ambiguity
				}
			}else{
				Contig 	ctg = scaffoldF.remove(index);
				ContigBridge brg = getReversedBridge(bridge);
				//extend and connect
				while(true){
					ctg.composite(rev); // contigF->headF + headF->ctg = contigF->ctg
					ctg.composite(trans); // contigT->contigF + contigF->ctg = contigT->ctg
					ctg.composite(contigT.getVector()); //headT->contigT + contigT->ctg = headT->ctg : relative position of this ctg w.r.t headT
	
					ctg.head = headT;
					//if (posT == 1){
					if(!firstDir){
						scaffoldT.addFront(ctg,brg);
					}else{
						scaffoldT.addRear(ctg,getReversedBridge(brg));
					}	
					if(count<1) break;
					ctg = scaffoldF.remove(--count);
					brg = scaffoldF.bridges.remove(count);
					
				}
				if(scaffoldF.closeBridge!=null && !scaffoldF.isEmpty()){
					count = scaffoldF.size()-1;
					ctg = scaffoldF.removeLast();
					brg = scaffoldF.closeBridge;
					
					while(true){
						ctg.myVector = ScaffoldVector.composition(ScaffoldVector.reverse(scaffoldF.circle),ctg.myVector);
						ctg.composite(rev); // contigF->headF + headF->ctg = contigF->ctg
						ctg.composite(trans); // contigT->contigF + contigF->ctg = contigT->ctg
						ctg.composite(contigT.getVector()); //headT->contigT + contigT->ctg = headT->ctg : relative position of this ctg w.r.t headT
						//ctg.composite(ScaffoldVector.reverse(scaffoldF.circle)); //composite co tinh giao hoan k ma de day???
						ctg.head = headT;
						//if (posT == 1){ 
						if(!firstDir){
							scaffoldT.addFront(ctg,brg);
						}else{
							scaffoldT.addRear(ctg,getReversedBridge(brg));
						}	
						if(count<1) break;
						brg = scaffoldF.bridges.remove(count--);
						ctg = scaffoldF.remove(count);		
						
					}
				}
			
				//set the remaining
				scaffoldT.trim();
				if(scaffoldF.size() > 0){
					scaffoldF.trim();
					addScf=scaffoldF.getFirst().getIndex();
					changeHead(scaffoldF, scaffoldF.getFirst());
				}
				
			}
			scaffoldF = new Scaffold(contigs.get(headF));
		}
		else if(secondDir == 1){
			if(headF==headT){
					//if(posT!=-1)
					if(!firstDir)
						return false;
					else{
						Contig prevMarker = scaffoldF.nearestMarker(contigF, false);
						if(prevMarker!=null){
							Contig ctg = scaffoldF.remove(--count);
							Scaffold newScf = new Scaffold(ctg);
							ContigBridge brg = scaffoldF.bridges.remove(count);
							while(true){
								if(count<1) break;
								ctg= scaffoldF.remove(--count);
								brg = scaffoldF.bridges.remove(count);
								newScf.addFront(ctg,brg);
							}
							newScf.trim();
							changeHead(newScf, prevMarker);
							addScf=prevMarker.getIndex();
						}
						scaffoldF.setCloseBridge(bridge);
						changeHead(scaffoldF, contigF);
						
					}
			}else{
				Contig 	ctg = scaffoldF.remove(index);
				ContigBridge brg = bridge;
				//extend and connect
				while(true){
					ctg.composite(rev); // contigF->headF + headF->ctg = contigF->ctg
					ctg.composite(trans); // contigT->contigF + contigF->ctg = contigT->ctg
					ctg.composite(contigT.getVector()); //headT->contigT + contigT->ctg = headT->ctg : relative position of this ctg w.r.t headT
	
					ctg.head = headT;
					//if (posT == 1){ 
					if(!firstDir){
						scaffoldT.addFront(ctg,getReversedBridge(brg));
					}else{
						scaffoldT.addRear(ctg,brg);
					}				
					if(scaffoldF.size()==index) break;
					ctg = scaffoldF.remove(index);
					brg = scaffoldF.bridges.remove(index);
				}
				if(scaffoldF.closeBridge!=null && !scaffoldF.isEmpty()){
					ctg = scaffoldF.removeFirst();
					brg = scaffoldF.closeBridge;
					while(true){
						ctg.myVector = ScaffoldVector.composition(scaffoldF.circle,ctg.myVector);
						ctg.composite(rev); // contigF->headF + headF->ctg = contigF->ctg
						ctg.composite(trans); // contigT->contigF + contigF->ctg = contigT->ctg
						ctg.composite(contigT.getVector()); //headT->contigT + contigT->ctg = headT->ctg : relative position of this ctg w.r.t headT
						//ctg.composite(scaffoldF.circle);
						ctg.head = headT;
						//if (posT == 1){ 
						if(!firstDir){
						scaffoldT.addFront(ctg,getReversedBridge(brg));
						}else{
							scaffoldT.addRear(ctg,brg);
						}	
						if(scaffoldF.size()<1) break;
						brg = scaffoldF.bridges.removeFirst();
						ctg = scaffoldF.removeFirst();		
					}	
				}
				//set the remaining
				scaffoldT.trim();
				if(scaffoldF.size() > 0){
					scaffoldF.trim();
					addScf=scaffoldF.getLast().getIndex();
					changeHead(scaffoldF, scaffoldF.getLast());
				}
				
			}
			scaffoldF = new Scaffold(contigs.get(headF));
		}	
		else
			return false;
		
		//===================================================================================================
		if(verbose){ 
			System.out.println("After Joining: " + (addScf<0?1:2) + " scaffolds!");
			scaffolds[contigF.head].view();
			if(addScf >=0)
				scaffolds[addScf].view();
		}
		return true;
	}
	//change head of scaffold scf to newHead
	public void changeHead(Scaffold scf, Contig newHead){	
		if(isRepeat(newHead)){
			if(verbose)
				System.out.println("Cannot use repeat as a head! " + newHead.getName());
			return;
		}
		//Scaffold scf = scaffolds[scfIndex];
		int scfIndex = scf.scaffoldIndex;
		int headPos = scf.indexOf(newHead);
		if(headPos < 0){
			if(verbose)
				System.out.printf("Cannot find contig %d in scaffold %d\n" , newHead.getIndex(), scfIndex);
			return;
		}
		Scaffold newScf = new Scaffold(newHead.getIndex());
		ScaffoldVector rev = ScaffoldVector.reverse(newHead.getVector()); //rev = newHead->head	
		
		if(newHead.getRelDir() == 0){
			if(verbose)
				System.out.printf("Contig %d of scaffold %d got direction 0!\n" , newHead.getIndex(), scfIndex);
			return;
		}
		else if(newHead.getRelDir() > 0){
			while(!scf.isEmpty())
				newScf.add(scf.removeFirst());
			while(!scf.bridges.isEmpty())
				newScf.bridges.add(scf.bridges.removeFirst());
			if(scf.closeBridge != null){
				newScf.closeBridge = scf.closeBridge;
				newScf.circle = scf.circle;
			}
		}
		else{
			while(!scf.isEmpty())
				newScf.add(scf.removeLast());
			while(!scf.bridges.isEmpty())
				newScf.bridges.add(getReversedBridge(scf.bridges.removeLast()));
			if(scf.closeBridge != null){
				newScf.closeBridge = getReversedBridge(scf.closeBridge);
				newScf.circle = ScaffoldVector.reverse(scf.circle);
			}
		}

		for (Contig ctg:newScf){					
			ctg.composite(rev); // leftmost->head + head->ctg = leftmost->ctg
		}
		newScf.setHead(newHead.getIndex());
		scaffolds[newHead.getIndex()] = newScf;
		
	}
	public synchronized void printSequences() throws IOException{
		if(annotation){
			SequenceOutputStream aout = SequenceOutputStream.makeOutputStream(prefix+".anno.japsa");
			for (int i = 0; i < scaffolds.length;i++){
				if(scaffolds[i].isEmpty()) continue;
				int len = scaffolds[i].getLast().rightMost() - scaffolds[i].getFirst().leftMost();
				if ((contigs.get(i).head == i 
					&& !isRepeat(contigs.get(i))
					&& len > maxRepeatLength
					)
					|| scaffolds[i].closeBridge != null
					)
				{
					if(verbose) 
						System.out.println("Scaffold " + i + " estimated length " + len);
					scaffolds[i].viewAnnotation(aout);
				}
			}
			aout.close();
		} else{
			SequenceOutputStream 	fout = SequenceOutputStream.makeOutputStream(prefix+".fin.fasta"),
									jout = SequenceOutputStream.makeOutputStream(prefix+".fin.japsa");
			for (int i = 0; i < scaffolds.length;i++){
				if(scaffolds[i].isEmpty()) continue;
				int len = scaffolds[i].getLast().rightMost() - scaffolds[i].getFirst().leftMost();
				if ((contigs.get(i).head == i 
					&& !isRepeat(contigs.get(i))
					&& len > maxRepeatLength
					)
					|| scaffolds[i].closeBridge != null
					)
				{
					if(verbose) 
						System.out.println("Scaffold " + i + " estimated length " + len);
					
					scaffolds[i].viewSequence(fout, jout);
				}
			}
			fout.close();
			jout.close();
		}
	}	
	public synchronized void printRT(long tpoint) throws IOException{
		for (Contig contig:contigs){
			if(contig.oriRep.size() > 0){
				String fname = contig.getName() + ".rtout";
				File f = new File(fname);
				if(!f.exists())
					f.createNewFile();
				
				//BufferedWriter out = new BufferedWriter(new FileWriter(f.getPath(), true));
				FileWriter fw = new FileWriter(f,true);
				BufferedWriter bw = new BufferedWriter(fw);
				PrintWriter pw = new PrintWriter(bw);
				
				ArrayList<String> 	ctgList = new ArrayList<String>(),
									origList = new ArrayList<String>(), 
									resList = new ArrayList<String>(),
									genesList = new ArrayList<String>();
				
				for(Contig ctg:scaffolds[contig.head]){
					ctgList.add(ctg.getName());
					if(ctg.oriRep.size()>0)
						for(JapsaFeature ori:ctg.oriRep)
							origList.add(ori.getID());
					for (JapsaFeature feature:ctg.genes)
						genesList.add(feature.toString());
					for (JapsaFeature feature:ctg.resistanceGenes)
						resList.add(feature.toString());
				}
				float streamData=tpoint/1000000;
				pw.print(">");
				for(String ctg:ctgList)
					pw.printf("%s\t", ctg);
				
				pw.printf("\n>%.2fMpb\t%d genes\t", streamData, genesList.size());
				
				for(String ori:origList)
					pw.printf("+%s", ori);
				
				for(String genes:genesList)
					pw.print(" \n\t"+genes);
				pw.println("");
				
				for(String res:resList)
					pw.print(" \n\t"+res);
				pw.println("");
				
				pw.close();

			}
		}


	}	
	// To check if this contig is likely a repeat or a singleton. If FALSE: able to be used as a milestone.
	public static boolean isRepeat(Contig ctg){
		//for the case of AbySS when no coverage information of contigs is found
		if(estimatedCov == 1.0 && ctg.getCoverage() == 1.0){
			if(ctg.length() > maxRepeatLength)
				return false;
			else
				return true;
		}
		
		if (ctg.length() < minContigLength || ctg.getCoverage() < .3 * estimatedCov) return true;
		else if (ctg.length() > maxRepeatLength || ctg.getCoverage() < 1.3 * estimatedCov) 
			return false; 
		else if (ctg.getCoverage() > 1.5 * estimatedCov)
			return true;
		else{
			for(ContigBridge bridge:ctg.bridges){
				Contig other = bridge.firstContig.getIndex()==ctg.getIndex()?bridge.secondContig:bridge.firstContig;
				if(other.getIndex()==ctg.getIndex()) continue;
				int dist=bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig);
				if( dist<0 && dist>-ctg.length()*.25){
					if(other.length() > maxRepeatLength || other.getCoverage() < 1.3*estimatedCov)
						return true;
				}
			}
			
		}
		if(ctg.length() < 2*minContigLength) // further filter: maybe not repeat but insignificant contig 
			return true;
		else 
			return false;
	}


	public void connectBridges() {
		//TODO: to be override by children classes
	}


}