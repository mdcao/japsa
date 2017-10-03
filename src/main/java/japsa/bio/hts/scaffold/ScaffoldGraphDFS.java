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
 * 3 Jan 2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.bio.hts.scaffold;

import japsa.seq.JapsaFeature;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author minhduc&sonnguyen
 *
 */
public class ScaffoldGraphDFS extends ScaffoldGraph {
	/**
	 * @param sequenceFile
	 * @throws IOException
	 */
	public ScaffoldGraphDFS(String sequenceFile, String genesFile, String resistFile, String isFile, String oriFile) throws IOException, InterruptedException {
		super(sequenceFile);
		if(resistFile != null){
			readDb(resistFile, "Resistance genes", .9, 1.0);
			annotation = true;
		}
		if(isFile != null){
			readDb(isFile, "Insertion sites", .8, .9);
			annotation = true;
		}
		if(oriFile != null){
			readDb(oriFile, "Origin of replication", .8, .9);
			annotation = true;
		}
		if(genesFile != null){
			readGFF(genesFile);
			annotation = true;
		}
//		for (Contig contig:contigs){
//			for(JapsaFeature feature:contig.genes)
//				System.out.println(contig.getName() + "\t" + feature.getStrand() + "\t" + feature);
//		}
//		for (Contig contig:contigs){
//			for(JapsaFeature feature:contig.oriRep)
//				System.out.println(contig.getName() + "\t" + feature.getStrand() + "\t" + feature);
//		}
//		for (Contig contig:contigs){
//			for(JapsaFeature feature:contig.insertSeq)
//				System.out.println(contig.getName() + "\t" + feature.getStrand() + "\t" + feature);
//		}
	}

	private void readDb(String data, String type, double minCov, double minID) throws IOException, InterruptedException{
		type = type.toLowerCase();
		
		String blastn = "blastn";

		ProcessBuilder pb = new ProcessBuilder(blastn, 
			"-subject", 
			"-",
			"-query", 
			data, 
			"-outfmt", 
			"7 qseqid qlen qstart qend sseqid slen sstart send length frames pident nident gaps mismatch score bitscore sstrand");
		/////    0      1     2     3    4     5     6     7     8      9      10     11    12    13       14     15		16
		
		//6 qseqid qlen length pident nident gaps mismatch
		//    0      1    2      3      4      5    6
		Process process = pb.start();
		//Pass on the genome to blastn
		SequenceOutputStream out = new SequenceOutputStream(process.getOutputStream());
		for (Contig ctg:contigs){
			Sequence seq=ctg.contigSequence;
			seq.writeFasta(out);
		}
		out.close();

		//Read the output of blastn
		BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));
		String line;

		while ((line = br.readLine()) != null) {
			if (line.startsWith("#"))
				continue;

			String [] toks = line.trim().split("\t");
			int length = Integer.parseInt(toks[8]);
			int qlen   = Integer.parseInt(toks[1]);
			double cov = (float)length/qlen;
			if (minCov > cov){
				continue;
			}

			if (Double.parseDouble(toks[10]) < minID * 100){
				continue;
			}
			//pass			
			
			Contig ctg = getSPadesContig(toks[4]);
			if(ctg != null){
				char strand = toks[16].equals("plus")?'+':'-';
				JapsaFeature feature = new JapsaFeature(Integer.parseInt(toks[6]), Integer.parseInt(toks[7]), type, toks[0], strand, ctg.getName());

				feature.addDesc(toks[0]+ ":" + (int)(cov*100) + "% cover, " + toks[10] + "% identity");

				switch (type.toLowerCase()){
					case "resistance genes":
						ctg.resistanceGenes.add(feature);
						break;
					case "insertion sites":
						ctg.insertSeq.add(feature);
						break;
					case "origin of replication":
						ctg.oriRep.add(feature);
						break;
					default:
						System.err.println(type + " has not yet included in our analysis!");
						break;

				}
				//Rewritten the below with the above
				//switch (type.toLowerCase()){
				//	case "resistance genes":
				//		ctg.resistanceGenes.add(feature);
				//		break;
				//	case "insertion sites":
				//		ctg.insertSeq.add(feature);
				//		break;
				//	case "origin of replication":
				//		ctg.oriRep.add(feature);
				//		break;
				//	default:
				//		System.err.println(type + " has not yet included in our analysis!");
				//		break;
				//}
			}
			Collections.sort(ctg.resistanceGenes);
			Collections.sort(ctg.insertSeq);
			Collections.sort(ctg.oriRep);
		}
		br.close();
		//process.waitFor();//Do i need this???
	}
	private void readGFF(String fileName) throws IOException, InterruptedException{
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line;

		while ((line = br.readLine()) != null) {
			if (line.startsWith("#"))
				continue;
			if (line.startsWith(">"))
				break;
			String [] toks = line.trim().split("\t");
			Contig ctg = getSPadesContig(toks[0]+"_"); //get the contig from its shorten name
			if(ctg != null){
				int start = Integer.parseInt(toks[3]),
					end = Integer.parseInt(toks[4]);
				String [] des = toks[8].trim().split(";");
				String [] id = des[0].trim().split("=");
				String ID = "undefined";
				if(id[0].equals("ID"))
					ID = id[1];
				if(!toks[2].equals("gene")){
					JapsaFeature feature = new JapsaFeature(start, end, toks[2], ID, toks[6].charAt(0), ctg.getName());
					feature.addDesc(toks[8]);
					ctg.genes.add(feature);
				}
			}
		}
		br.close();
	}
	public ScaffoldGraphDFS(String sequenceFile, String graphFile) throws IOException {
		super(sequenceFile);
		//TODO: implement fastg reader for SequenceReader to have pre-assembled bridges

	}
	
	/* 
	 * Arrange contigs
	 * 
	 */
	@Override
	public synchronized void connectBridges(){
//		for(Contig ctg:contigs){
//			//System.out.println(ctg.getName());
//			if(ctg.isCircular){
//				if(ctg.bridges.size() != 2){
//					ctg.isCircular = false;
//					break;
//				}
//				for(ContigBridge brg:ctg.bridges){
//					//System.out.println("\t"+brg.hashKey+" : "+brg.getTransVector());
//					if(brg.firstContig.getIndex() != brg.secondContig.getIndex()){
//						ctg.isCircular = false;
//						break;
//					}
//				}
//			}
//		}
		
		// Start scaffolding
		if(verbose) 
			System.out.println("Starting scaffolding.......");
		
		List<LengthIndex> list = new ArrayList<LengthIndex>();
		for(int i = 0; i<contigs.size(); i++)
			list.add(new LengthIndex(contigs.get(i).length(),i));
		Collections.sort(list);
		for(LengthIndex idx:list){
			int i = idx.index;
			//TODO: include also the case of circularized scaffold (for realtime self-correction)
			if (contigs.get(i).head !=i || scaffolds[i].size() < 1 || scaffolds[i].closeBridge != null)
				continue;
			//Now extend scaffold i				
			if(	isRepeat(scaffolds[i].element()) 
				|| scaffolds[i].element().length() < minContigLength){
				if(!scaffolds[i].element().isCircular())
					continue;
			}
			//1.a extend to the first
			boolean closed = false;
			/*****************************************************************/
			////////////////////////////First try to extend to the end/////////////////////////////
			if(verbose) 
				System.out.printf("Extending %d to the rear\n",i);
			closed = walk2(i,true);

			////////////////////////////Then extend to the front///////////////////////////////////

			if (!closed){
				if(verbose) 
					System.out.printf("Extending %d to the front\n",i);
				closed = walk2(i,false);
			}
			//scaffolds[i].trim();
			/*****************************************************************/
			if(verbose && scaffolds[i].size() > 1) 
				System.out.printf("Finally, scaffold %d size %d  and is %s\n",
					i,
					scaffolds[i].getLast().rightMost() - scaffolds[i].getFirst().leftMost(),
					closed?"circular":"linear");			
		}//for

	}

	
	/* 	To check the agreement of locations between 3 contigs: marker -> prevContig -> curContig
	 *	(also their corresponding scores?). 
	 * 	Purpose: avoid false positive alignments.
	 *	
	 */
	private ContigBridge checkHang(Contig marker, ContigBridge toPrev, ContigBridge toCurrent){
		Contig 	prevContig = toPrev.secondContig,
				curContig = toCurrent.secondContig;
		ScaffoldVector 	mark2Prev = toPrev.getTransVector(),
						mark2Cur = toCurrent.getTransVector();
		ScaffoldVector prevToCur = ScaffoldVector.composition(mark2Cur,ScaffoldVector.reverse(mark2Prev));
		if(verbose) 
			System.out.printf("\texamining %s to %s:\n ",prevContig.getName(),curContig.getName());
		for(ContigBridge brg:getListOfBridgesFromContig(prevContig)){
			if(brg.secondContig.getIndex() == curContig.getIndex()){
				if(brg.consistentWith(prevToCur)){
					if(verbose) 
						System.out.printf("=> consistent between bridge vector %s and estimated vector %s\n",brg.getTransVector(),prevToCur);
					return brg; 
				}
				else{
					if(verbose) 
						System.out.printf("=> inconsistent between bridge vector %s and estimated vector %s\n",brg.getTransVector(),prevToCur);
				}
			}

		}	
		
		return null;
	}
	
	/*
	 * Walking through markers while trying to fill the gap by repeat sequences simultaneously
	 */
	private boolean walk2(int i, boolean direction ){	
		Scaffold scaffold = scaffolds[i];					
		boolean extended = true;
		boolean closed = scaffold.closeBridge!=null;

		/*****************************************************************/
		while (extended && (!closed) && scaffold.size() > 0){
			Contig ctg = direction?scaffold.getLast():scaffold.getFirst();		
			ArrayList<ContigBridge> bridges=getListOfBridgesFromContig(ctg);

			if(verbose) {
				System.out.printf(" Last of scaffold %d extention is on contig %d (%s): ",i,ctg.getIndex(),ctg.getName());
				System.out.printf("iterating among %d bridges\n",bridges.size());
			}
			int ctgEnd = direction?ctg.rightMost():ctg.leftMost();
			 
			extended = false; //only continue the while loop if extension is on the move (line 122)
			int extendDir = 0, //direction to go on the second scaffold: ScaffoldT (realtime mode)
				curStep = Integer.MAX_VALUE; //distance between singleton1 -> singleton2
			double	curScore = 0.0; //score between singleton1 -> singleton2
			ContigBridge stepBridge = null;
			
//			ArrayList<Contig> extendableContig = new ArrayList<Contig>(maxLink);
//			ArrayList<ContigBridge> extendableContigBridge = new ArrayList<ContigBridge>(maxLink);
//			ArrayList<ScaffoldVector> extendableVector = new ArrayList<ScaffoldVector>(maxLink);
//			ArrayList<Integer> distances = new ArrayList<Integer>(maxLink);
			ArrayList<ExtendCandidate> candidates = new ArrayList<ExtendCandidate>();
			
			Collections.sort(bridges);
			for (ContigBridge bridge:bridges){
				if (bridge.firstContig == bridge.secondContig) //2 identical markers ??!
					if(!bridge.firstContig.isCircular())
						continue;
				Contig nextContig = bridge.secondContig;					
				ScaffoldVector trans = bridge.getTransVector();
//				if (ctg == bridge.secondContig){
//					nextContig = bridge.firstContig;
//					trans = ScaffoldVector.reverse(trans);
//				}					
				if(verbose) 
					System.out.println("..." + nextContig.getName());
				ScaffoldVector trialTrans = ScaffoldVector.composition(trans, ctg.getVector());
				int newEnd = direction?nextContig.rightMost(trialTrans):nextContig.leftMost(trialTrans);				
				
				//see if the next contig would extend the scaffold to the right
				//only take one next singleton (with highest score possible sorted) as the marker for the next extension
				int distance = bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig);
				if (direction?(newEnd > ctgEnd):(newEnd < ctgEnd)){	
					int aDir = 0;

					if(!isRepeat(nextContig) || (ctg.isCircular() && ctg.getIndex() == nextContig.getIndex())){
						//check quality of the bridge connected 2 markers
						if(scaffolds[nextContig.head].size() > 1){
							aDir = extendDirection(ctg, bridge);
							if(aDir==0){
								if(verbose)
									System.out.println("No jump to "  + nextContig.getName());
								continue;
							}
						}
						if(distance > -maxRepeatLength && bridge.getConnections().size() >= minSupportReads && bridge.getScore() > curScore ){
							if(verbose)
								bridge.display();
							curStep = distance;
							curScore = bridge.getScore();
							stepBridge = bridge;
							extendDir = aDir;
						}else{
							if(verbose)
								System.out.printf("Cannot form unique bridge from %d to %d with %d connections and score %.2f\n", 
												ctg.getIndex(), nextContig.getIndex(), bridge.getConnections().size(), bridge.getScore());
							continue;
						}
					}
					
					if(verbose) 
						System.out.printf(" Might extend %d from %d(%d) to %d(%d) (%s direction) with score %f and distance %d\n"
										,i,ctg.index, ctgEnd, nextContig.index, newEnd, 
										(bridge.getTransVector().getDirection() > 0?"same":"opposite"), bridge.getScore(), distance);
					
//					int j = 0;
					//looking for right position to have the list sorted
//					for(j=0; j<distances.size(); j++)
//						if(distances.get(j) > distance)
//							break;
					candidates.add(new ExtendCandidate(distance, nextContig, bridge, trialTrans));
//					distances.add(j,distance);
//					extendableContig.add(j, nextContig);
//					extendableContigBridge.add(j, bridge);
//					extendableVector.add(j, trialTrans);
				}
				else if(verbose) 
						System.out.printf(" No extend %d from %d(%d) to %d(%d) with score %f and distance %d\n",i,ctg.index, ctgEnd, nextContig.index, newEnd,bridge.getScore(), distance);
				
			}//for	
			Collections.sort(candidates);
			
			int noOfUniqueContig = 0; //reset to count how many singleton will be added now
			
			if(stepBridge==null){
				if(verbose) 
					System.out.printf(" Extension of Scaffold %d toward stopped at %d due to the lack of next marker!\n", i,ctg.index);
				return false;
			}
			int curEnd = ctgEnd;
			Contig prevContig = ctg;
			ContigBridge 	prevContigBridge = null;
			ScaffoldVector prevVector = new ScaffoldVector();
			for(int index = 0; index < candidates.size(); index++){
				ExtendCandidate cdd = candidates.get(index);
				if(cdd.getDistance() > curStep)
					continue;
				Contig curContig = cdd.getContig();
				ContigBridge curContigBridge = cdd.getBridge(); //will be replaced by the bridge to prev contig later
				ScaffoldVector curVector = cdd.getVector();
				if(verbose) 
					System.out.println("Checking contig " + curContig.getName() + "...");
				if(	isRepeat(curContig) && !curContig.isCircular())
					if(checkHang(ctg, curContigBridge, stepBridge)==null)
						continue;
				prevVector = prevContig.getVector();
				boolean extendable = false;
				ScaffoldVector prevToCur = ScaffoldVector.composition(curVector,ScaffoldVector.reverse(prevVector));
				ContigBridge confirmedBridge;
				//TODO: if this happen with singleton -> chimeric happen (unique+repeat=contig) need to do smt...
				if(	isRepeat(curContig) &&
					(direction?(curContig.rightMost(curVector) < curEnd):(curContig.leftMost(curVector)) > curEnd)){
					if(verbose) 
						System.out.println(curContig.getName() + " is ignored because current end " + curEnd + 
										" cover the new end " + (direction?curContig.rightMost(curVector): curContig.leftMost(curVector)));
					continue;
				}
				else{	
					if(index >= 1 && prevContigBridge != null){
						confirmedBridge = checkHang(ctg, prevContigBridge, curContigBridge);
						if(confirmedBridge != null){
							prevContigBridge = curContigBridge;
							prevToCur =confirmedBridge.getTransVector();
							extendable = true;
						}
						else
							continue;
					}
					else{
						prevContigBridge = curContigBridge;
						confirmedBridge = curContigBridge;
						extendable = true;
					}
				}		
				if(extendable){
					// if extension is circularized
					if(curContig.getIndex() == (direction?scaffold.getFirst().getIndex():scaffold.getLast().getIndex())
						&& (!isRepeat(curContig) || curContig.isCircular())
						){
						if(verbose) 
							System.out.printf(" *****************SCAFFOLD %d CLOSED AFTER CONNECT %d ***********************\n", i,curContig.index);
						scaffold.setCloseBridge(direction?confirmedBridge:getReversedBridge(confirmedBridge));
						curContigBridge.setContigScores();
						return true;
					}
					
					if(isRepeat(curContig)){
						curContig.head = i; //must be here!
						curContig = curContig.clone();
					}else{
						//check to join 2 scaffolds and stop this round
						if (scaffolds[curContig.head].size() > 1){
							if(!joinScaffold(prevContig,confirmedBridge,direction,extendDir)){
								if(verbose)
									System.out.printf(" Skip to connect contig %d of %d to contig %d of %d\n", ctg.index,i,curContig.index, curContig.head);
								continue;
							}			
							else{
								curContigBridge.setContigScores();
								noOfUniqueContig++;
								break;
							}

						}
						curContig.head = i; //must be here!
						noOfUniqueContig++;
						curContigBridge.setContigScores();
					}

					if(verbose) {						
						System.out.printf(" Extend %d from %d(%d) to %d(%d) with score %f: ",
							i,ctg.index, ctgEnd, curContig.index, curContig.rightMost(curVector), curContigBridge.getScore());						
						System.out.printf(" curContigBridge %d -> %d\n", confirmedBridge.firstContig.getIndex(), curContigBridge.secondContig.getIndex());
					}
					curContig.myVector = ScaffoldVector.composition(prevToCur,prevContig.getVector());//from the head contig
					//confirmedBridge=confirmedBridge.clone(prevContig,curContig);

					if(direction)
						scaffolds[i].addForward(curContig, confirmedBridge);
					else
						scaffolds[i].addBackward(curContig, getReversedBridge(confirmedBridge));
					
					
					curEnd = direction?curContig.rightMost(curVector):curContig.leftMost(curVector);
					extended = true; //scaffold extension is really on the move...									
					
					prevContig = curContig;
					
					if(verbose)
						scaffolds[i].view();

				}
				
				if(cdd.getDistance() == curStep)
					break;
			}//for
			if(noOfUniqueContig < 1){
				if(verbose) 
					System.out.printf(" Extension of Scaffold %d toward stopped at %d because next marker is not reachable!\n", i,ctg.index);
				scaffolds[i].trim();
				return false;
			}
			// TO THE NEXT UNIQUE CONTIG
		}//while
		return closed;
	}
	
	class ExtendCandidate implements Comparable<ExtendCandidate>{
		int distance;
		Contig cContig;
		ContigBridge cBridge;
		ScaffoldVector cVector;
		ExtendCandidate(int dist, Contig ctg, ContigBridge brg, ScaffoldVector vect){
			distance=dist;
			cContig=ctg;
			cBridge=brg;
			cVector=vect;
		}
		int getDistance(){return distance;}
		Contig getContig(){return cContig;}
		ContigBridge getBridge(){return cBridge;}
		ScaffoldVector getVector(){return cVector;}
		@Override
		public int compareTo(ExtendCandidate o) {
			// TODO Auto-generated method stub
			return distance-o.distance;
		}
		
	}
	class LengthIndex implements Comparable<LengthIndex>{
		int length, index;
		public LengthIndex(int len, int index){
			this.length = len;
			this.index = index;
		}
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(LengthIndex o) {
			return (int) (o.length - length);

		}	
	}
}
