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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author minhduc&sonnguyen
 *
 */
public class ScaffoldGraphDFS extends ScaffoldGraph {
	boolean [] scaffoldExtended;
	/**
	 * @param sequenceFile
	 * @throws IOException
	 */
	public ScaffoldGraphDFS(String sequenceFile) throws IOException {
		super(sequenceFile);
		scaffoldExtended = new boolean[scaffolds.length];
		Arrays.fill(scaffoldExtended, false);
	}
	
	public ScaffoldGraphDFS(String sequenceFile, String graphFile) throws IOException {
		super(sequenceFile);
		//TODO: implement fastg reader for SequenceReader to have pre-assembled bridges
		scaffoldExtended = new boolean[scaffolds.length];
		Arrays.fill(scaffoldExtended, false);
	}
	
	private boolean isUniqueBridge(ContigBridge bridge){
		if(!isRepeat(bridge.firstContig) &&
			!isRepeat(bridge.secondContig))
			return true;
		else 
			return false;
	}
	/* 
	 * Arrange contigs
	 * mode == true: fill the gaps between markers with fillers.
	 * mode == false: only arrange the markers
	 */
	@Override
	public void connectBridges(boolean mode){
		System.out.println(bridgeList.size());
		//Add bridges for DFS
		//Note that the bridges will be ordered by scores in contig.bridges due to sorted bridgeList
		for (ContigBridge bridge:bridgeList){		
			bridge.firstContig.bridges.add(bridge);
			bridge.secondContig.bridges.add(bridge);
		}
		//TODO more about repeats(number of occurrences, plasmid contig self-scaffold...)
		// Start scaffolding
		System.out.println("Starting scaffolding in " + (mode?"fast":"full") + " mode:");
		
		List<LengthIndex> list = new ArrayList<LengthIndex>();
		for(int i = 0; i<contigs.size(); i++)
			list.add(new LengthIndex(contigs.get(i).length(),i));
		Collections.sort(list);
		//for (int i = 0; i < scaffolds.length;i++){
		for(LengthIndex idx:list){
			int i = idx.index;
			if (scaffoldExtended[i])
				continue;
			//Now extend scaffold i				
			if(isRepeat(scaffolds[i].element()) || scaffolds[i].element().length() < 1000){
				if(!scaffolds[i].element().isCircular)
					continue;
			}
			//1.a extend to the first
			boolean closed = false;
			/*****************************************************************/
			////////////////////////////First try to extend to the end/////////////////////////////
			System.out.printf("Extending %d to the rear\n",i);
			if(!mode)
				closed = walk2(i,true);
			else
				closed = walk(i,true);
			////////////////////////////Then extend to the front///////////////////////////////////

			if (!closed){
				System.out.printf("Extending %d to the front\n",i);
				if(!mode)
					closed = walk2(i,false);
				else
					closed = walk(i, false);
			}
			/*****************************************************************/
			System.out.printf("Finally, scaffold %d size %d  and is %s\n",
					i,
					scaffolds[i].getLast().rightMost() - scaffolds[i].getFirst().leftMost(),
					closed?"circular":"linear");			
		}//for

	}

	
	/* 	To check the agreement of locations between 3 contigs: singleton (milestone) -> prevContig -> curContig
	 *	(also their corresponding scores?). 
	 * 	Purpose: avoid false positive alignments.
	 *	
	 */
	private ContigBridge checkHang(Contig singleton, ContigBridge toPrev, ContigBridge toCurrent){
		Contig 	prevContig = toPrev.firstContig.getIndex()==singleton.getIndex()?toPrev.secondContig:toPrev.firstContig,
				curContig = toCurrent.firstContig.getIndex()==singleton.getIndex()?toCurrent.secondContig:toCurrent.firstContig;
		ScaffoldVector 	sing2Prev = toPrev.firstContig.getIndex()==singleton.getIndex()?toPrev.getTransVector():ScaffoldVector.reverse(toPrev.getTransVector()),
						sing2Cur = toCurrent.firstContig.getIndex()==singleton.getIndex()?toCurrent.getTransVector():ScaffoldVector.reverse(toCurrent.getTransVector());
		ScaffoldVector prevToCur = ScaffoldVector.composition(sing2Cur,ScaffoldVector.reverse(sing2Prev));
		System.out.printf("\texamining %s to %s:\n ",prevContig.getName(),curContig.getName());
		for(ContigBridge brg:prevContig.bridges){
			if(	brg.firstContig.getIndex() == prevContig.getIndex() &&
				brg.secondContig.getIndex() == curContig.getIndex()){
				if(brg.consistentWith(prevToCur)){
					System.out.printf("=> consistent between bridge vector %s and estimated vector %s\n",brg.getTransVector(),prevToCur);
					return brg;
				}
				else{
					System.out.printf("=> inconsistent between bridge vector %s and estimated vector %s\n",brg.getTransVector(),prevToCur);
				}
			}
					
			if(	brg.firstContig.getIndex() == curContig.getIndex() &&
				brg.secondContig.getIndex() == prevContig.getIndex()){						
				if(brg.consistentWith(ScaffoldVector.reverse(prevToCur))){
					System.out.printf("=> consistent between bridge vector %s and estimated vector %s\n",brg.getTransVector(),ScaffoldVector.reverse(prevToCur));
					return brg;
				}
				else{
					System.out.printf("=> inconsistent between bridge vector %s and estimated vector %s\n",brg.getTransVector(),ScaffoldVector.reverse(prevToCur));
				}
			}
		}	
		
		return null;
	}
	
	/*
	 * Walking without filling repeats to gap between markers
	 */
	private boolean walk(int i, boolean direction ){	
		ScaffoldDeque scaffold = scaffolds[i];					
		boolean extended = true;
		Contig farMost = direction?scaffold.getFirst():scaffold.getLast();
		/*****************************************************************/
		while (extended){
			Contig ctg = direction?scaffold.getLast():scaffold.getFirst();				
			System.out.printf(" Last of scaffold %d extention is on contig %d: ",i,ctg.index);
			System.out.printf("iterating among %d bridges\n",ctg.bridges.size());
			int ctgEnd = direction?ctg.rightMost():ctg.leftMost();
			 
			extended = false;
			for (ContigBridge bridge:ctg.bridges){

				Contig nextContig = bridge.secondContig;					
				ScaffoldVector trans = bridge.getTransVector();
				if (ctg == bridge.secondContig){
					nextContig = bridge.firstContig;
					trans = ScaffoldVector.reverse(trans);
				}					
				
				ScaffoldVector trialTrans = ScaffoldVector.composition(trans, ctg.getVector());
				int newEnd = direction?nextContig.rightMost(trialTrans):nextContig.leftMost(trialTrans);

				System.out.printf(" Might extend %d from %d(%d) to %s(%d) with score %f and distance %d\n",i,ctg.index, ctgEnd, nextContig.getName(), newEnd,
									bridge.getScore(), bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig));

				if (!isUniqueBridge(bridge))
					continue;
				if (nextContig == farMost){
					double ratio = (newEnd - (direction?farMost.rightMost():farMost.leftMost())) 
									/ (0.0 + ctgEnd - (direction?farMost.leftMost():farMost.rightMost()));

					if (ratio > 0.6){
						System.out.printf(" Yay! scaffold %d closed after rear connect %d %f!\n", i,nextContig.index,ratio);
						scaffold.setCloseBridge(bridge);
						return true;
					}
				}					
				
				if (scaffolds[head[nextContig.index]].size() != 1){
					System.out.printf(" Attempt to connect contig %d of %d to contig %d of %d failed\n", ctg.index,i,nextContig.index, head[nextContig.index]);
					continue;//for bridge						
				}
				
				//see if the next contig would extend the scaffold to the right
				if (direction?(newEnd > ctgEnd):(newEnd < ctgEnd)){
					System.out.printf(" Can extend %d from %d(%d) to %d(%d) with score %f\n",i,ctg.index, ctgEnd, nextContig.index, newEnd, bridge.getScore());						
					if(direction)
						scaffolds[i].addRear(nextContig, bridge);
					else
						scaffolds[i].addFront(nextContig, bridge);
					nextContig.myVector = trialTrans;
					head[nextContig.index] = i;
					scaffoldExtended[nextContig.index] = true;
					
					extended = true;						
					
					scaffolds[i].view();						
					break;//for
				}else{
					System.out.printf(" Not extend scaffold %d from contig %d(%d) to %d(%d) with score %f\n",i,ctg.index, ctgEnd, nextContig.index, newEnd,bridge.getScore());
				}
			}//for	
			
		}//while
		return false;
	}
	/*
	 * Walking through markers while trying to fill the gap by repeat sequences simultaneously
	 */
	private boolean walk2(int i, boolean direction ){	
		ScaffoldDeque scaffold = scaffolds[i];					
		boolean extended = true;
		boolean closed = false;

		/*****************************************************************/
		while (extended && (!closed)){
			Contig ctg = direction?scaffold.getLast():scaffold.getFirst();				
			System.out.printf(" Last of scaffold %d extention is on contig %d (%s): ",i,ctg.getIndex(),ctg.getName());
			System.out.printf("iterating among %d bridges\n",ctg.bridges.size());
			int ctgEnd = direction?ctg.rightMost():ctg.leftMost();
			 
			extended = false; //only continue the while loop if extension is on the move (line 122)
			int maxLink = ctg.bridges.size(),
				noOfUniqueContig = 0,
				step = Integer.MAX_VALUE; //distance between singleton1 -> singleton2
			ContigBridge stepBridge = null;
			
			ArrayList<Contig> extendableContig = new ArrayList<Contig>(maxLink);
			ArrayList<ContigBridge> extendableContigBridge = new ArrayList<ContigBridge>(maxLink);
			ArrayList<ScaffoldVector> extendableVector = new ArrayList<ScaffoldVector>(maxLink);
			ArrayList<Integer> distances = new ArrayList<Integer>(maxLink);
			for (ContigBridge bridge:ctg.bridges){
				if (bridge.firstContig == bridge.secondContig) //not include tandem repeat
					if(!bridge.firstContig.isCircular)
						continue;
				Contig nextContig = bridge.secondContig;					
				ScaffoldVector trans = bridge.getTransVector();
				if (ctg == bridge.secondContig){
					nextContig = bridge.firstContig;
					trans = ScaffoldVector.reverse(trans);
				}					
				System.out.println("..." + nextContig.getName());
				ScaffoldVector trialTrans = ScaffoldVector.composition(trans, ctg.getVector());
				int newEnd = direction?nextContig.rightMost(trialTrans):nextContig.leftMost(trialTrans);				
				
				//see if the next contig would extend the scaffold to the right
				//only take one next singleton (with highest score possible sorted) as the mark for the next extension
				int distance = bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig);
				if (direction?(newEnd > ctgEnd):(newEnd < ctgEnd)){	
					if(distance > step)
						continue;
					if(!isRepeat(nextContig) || (ctg.isCircular && ctg.getIndex() == nextContig.getIndex())){
						if(noOfUniqueContig < 1){
							noOfUniqueContig++;
							step = distance;
							stepBridge = bridge;
						}
						else 
							continue;
					}
					
					System.out.printf(" Might extend %d from %d(%d) to %d(%d) (%s direction) with score %f and distance %d\n"
										,i,ctg.index, ctgEnd, nextContig.index, newEnd, 
										(bridge.getTransVector().getDirection() > 0?"same":"opposite"), bridge.getScore(), distance);
					
					int j = 0;
					//looking for right position to have the list sorted
					for(j=0; j<distances.size(); j++)
						if(distances.get(j) > distance)
							break;

					distances.add(j,distance);
					extendableContig.add(j, nextContig);
					extendableContigBridge.add(j, bridge);
					extendableVector.add(j, trialTrans);
				}
				else{
					System.out.printf(" No extend %d from %d(%d) to %d(%d) with score %f and distance %d\n",i,ctg.index, ctgEnd, nextContig.index, newEnd,bridge.getScore(), distance);
				}
			}//for	
			noOfUniqueContig = 0; //reset to count how many singleton will be added now
			
			if(stepBridge==null){
				System.out.printf(" Extension of Scaffold %d toward stopped at %d due to the lack of next marker!\n", i,ctg.index);
				return false;
			}
			int curEnd = ctgEnd;
			Contig prevContig = ctg;
			ContigBridge 	prevContigBridge = null;
			ScaffoldVector prevVector = new ScaffoldVector();
			for(int index = 0; index < extendableContig.size(); index++){
				Contig curContig = extendableContig.get(index);
				ContigBridge curContigBridge = extendableContigBridge.get(index); //will be replaced by the bridge to prev contig later
				ScaffoldVector curVector = extendableVector.get(index);
				System.out.println("Checking contig " + curContig.getName() + "...");
				if(	isRepeat(curContig) && !curContig.isCircular)
					if(checkHang(ctg, curContigBridge, stepBridge)==null)
						continue;
				prevVector = prevContig.getVector();
				boolean extendable = false;
				ScaffoldVector prevToCur = ScaffoldVector.composition(curVector,ScaffoldVector.reverse(prevVector));
				
				//TODO: if this happen with singleton -> chimeric happen (unique+repeat=contig) need to do smt...
				if(	isRepeat(curContig) &&
					(direction?(curContig.rightMost(curVector) < curEnd):(curContig.leftMost(curVector)) > curEnd)){
					System.out.println(curContig.getName() + " is ignored because current end " + curEnd + 
										" cover the new end " + (direction?curContig.rightMost(curVector): curContig.leftMost(curVector)));
					continue;
				}
				else{	
					if(index >= 1 && prevContigBridge != null){
						ContigBridge confirmedBridge = checkHang(ctg, prevContigBridge, curContigBridge);
						if(confirmedBridge != null){
							prevContigBridge = curContigBridge;
							curContigBridge = confirmedBridge;
							prevToCur = confirmedBridge.firstContig.getIndex()==prevContig.getIndex()?confirmedBridge.getTransVector():ScaffoldVector.reverse(confirmedBridge.getTransVector());
							extendable = true;
						}
						else
							continue;
					}
					else{
						prevContigBridge = curContigBridge;
						extendable = true;
					}
				}		
				if(extendable){
					if(curContig.getIndex() == (direction?scaffold.getFirst().getIndex():scaffold.getLast().getIndex())
						//&& !isRepeat(curContig)
						){
						System.out.printf(" *****************SCAFFOLD %d CLOSED AFTER CONNECT %d ***********************\n", i,curContig.index);
						scaffold.setCloseBridge(curContigBridge);
						return true;
					}
					if (scaffolds[head[curContig.index]].size() > 1){
						System.out.printf(" Skip to connect contig %d of %d to contig %d of %d\n", ctg.index,i,curContig.index, head[curContig.index]);
						continue;				
					}
					System.out.printf(" Extend %d from %d(%d) to %d(%d) with score %f: ",
							i,ctg.index, ctgEnd, curContig.index, curContig.rightMost(curVector), curContigBridge.getScore());						
					System.out.printf(" curContigBridge %d -> %d\n", curContigBridge.firstContig.getIndex(), curContigBridge.secondContig.getIndex());
					if(isRepeat(curContig)){
						curContig = curContig.clone();
//						curContig = curContig.clone(curContigBridge);
						if(curContigBridge.firstContig.getIndex() == curContigBridge.secondContig.getIndex())
							scaffolds[head[curContig.index]].element().isCircular = false; // tandem!
					}

					if(direction)
						scaffolds[i].addRear(curContig, curContigBridge);
					else
						scaffolds[i].addFront(curContig, curContigBridge);
					curContig.myVector = ScaffoldVector.composition(prevToCur,prevContig.getVector());//from the head contig
					curEnd = direction?curContig.rightMost(curVector):curContig.leftMost(curVector);
					extended = true; //scaffold extension is really on the move...									
					scaffolds[i].view();
					
					prevContig = curContig;
					if(!isRepeat(curContig)){
						noOfUniqueContig++;
						head[curContig.index] = i;
						scaffoldExtended[curContig.index] = true;
						break;
					}
				}
			}
			if(noOfUniqueContig < 1){
				System.out.printf(" Extension of Scaffold %d toward stopped at %d because next marker is not reachable!\n", i,ctg.index);
				return false;
			}
			// TO THE NEXT UNIQUE CONTIG
		}//while
		return closed;
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
