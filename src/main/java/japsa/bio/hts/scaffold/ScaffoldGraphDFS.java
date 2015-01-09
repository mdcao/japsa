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
import java.util.Arrays;

/**
 * @author minhduc
 *
 */
public class ScaffoldGraphDFS extends ScaffoldGraph {

	/**
	 * @param sequenceFile
	 * @throws IOException
	 */
	public ScaffoldGraphDFS(String sequenceFile) throws IOException {
		super(sequenceFile);
		// TODO Auto-generated constructor stub
	}

	/* (non-Javadoc)
	 * @see japsa.bio.hts.scaffold.ScaffoldGraph#connectBridges()
	 */
	@Override
	public void connectBridges(){
		System.out.println(bridgeList.size());
		
		//Make
		for (ContigBridge bridge:bridgeList){
			bridge.firstContig.bridges.add(bridge);
			bridge.secondContig.bridges.add(bridge);
		}
		
		boolean [] scaffoldExtended = new boolean[scaffolds.length];
		Arrays.fill(scaffoldExtended, false);
		
		for (int i = 0; i < scaffolds.length;i++){
			if (scaffoldExtended[i])
				continue;
			//Now extend scaffold i
			
			ScaffoldDeque scaffold = scaffolds[i];					
			
			//1.a extend to the first
			boolean extended = true;
			boolean closed = false;
			/*****************************************************************/
			//First try to extend to the end
			System.out.printf("Extending %d to the rear\n",i);
			while (extended&& (!closed)){
				Contig ctg = scaffold.getLast();				
				System.out.printf(" Last %d on %d (%d): \n",i,ctg.index,ctg.bridges.size());				
				int ctgEnd = ctg.rightMost();
				 
				extended = false;
				for (ContigBridge bridge:ctg.bridges){
					Contig nextContig = bridge.secondContig;					
					ScaffoldVector trans = bridge.getTransVector();
					if (ctg == bridge.secondContig){
						nextContig = bridge.firstContig;
						trans = ScaffoldVector.reverse(trans);
					}					
					
					ScaffoldVector trialTrans = ScaffoldVector.composition(trans, ctg.getVector());
					int newEnd = nextContig.rightMost(trialTrans);
					if (nextContig == scaffold.getFirst()){
						double ratio = (newEnd - scaffold.getFirst().rightMost()) / (0.0 + ctgEnd - scaffold.getFirst().leftMost());
						if (ratio > 0.6){
							System.out.printf(" Yay! scaffold %d closed after rear connect %d %f!\n", i,nextContig.index,ratio);
							scaffold.setCloseBridge(bridge);
							closed = true;
							break;
						}
					}					
					
					if (scaffolds[head[nextContig.index]].size() != 1){
						System.out.printf(" Attempt to connect contig %d of %d to contig %d of %d failed\n", ctg.index,i,nextContig.index, head[nextContig.index]);
						continue;//for bridge						
					}
					
					//see if the next contig would extend the scaffold to the right
					if (newEnd > ctgEnd){
						System.out.printf(" Extend %d from %d(%d) to %d(%d) with score %f\n",i,ctg.index, ctgEnd, nextContig.index, newEnd, bridge.getScore());						
						scaffolds[i].addRear(nextContig, bridge);
						nextContig.myVector = trialTrans;
						head[nextContig.index] = i;
						scaffoldExtended[nextContig.index] = true;
						
						extended = true;						
						
						scaffolds[i].view();						
						break;//for
					}else{
						System.out.printf(" Not Extend %d from %d(%d) to %d(%d) with score %f\n",i,ctg.index, ctgEnd, nextContig.index, newEnd,bridge.getScore());
					}
				}//for				
			}//while
			
			extended = true;
			//extend to the front
			if (!closed)
				System.out.printf("Extending %d to the front\n",i);
			while (extended&& (!closed)){			
				Contig ctg = scaffold.getFirst();				
				System.out.printf(" First %d on %d (%d): \n",i,ctg.index,ctg.bridges.size());				
				int ctgStart = ctg.leftMost();
				 
				extended = false;
				for (ContigBridge bridge:ctg.bridges){
					Contig nextContig = bridge.secondContig;					
					ScaffoldVector trans = bridge.getTransVector();
					if (ctg == bridge.secondContig){
						nextContig = bridge.firstContig;
						trans = ScaffoldVector.reverse(trans);
					}
					
					ScaffoldVector trialTrans = ScaffoldVector.composition(trans, ctg.getVector());
					int newStart = nextContig.leftMost(trialTrans);					
					if (nextContig == scaffold.getLast()){
						double ratio = (scaffold.getLast().leftMost() - newStart) / (0.0 + scaffold.getLast().rightMost() - ctgStart);
						if (ratio > 0.6){
							System.out.printf(" Yay! scaffold %d closed after front connect %d %f!\n", 
									i,
									nextContig.index,
									ratio);
							scaffold.setCloseBridge(bridge);
							closed = false;
							break;
						}
					}	
					
					if (scaffolds[head[nextContig.index]].size() != 1){
						System.out.printf(" Attempt to connect contig %d of %d to contig %d of %d failed\n", 
								ctg.index,
								i,
								nextContig.index, 
								head[nextContig.index]);
						continue;//for bridge						
					}
					
					if (newStart < ctgStart){
						System.out.printf(" Extend %d from %d(%d) to %d(%d) with score %f\n",
								i,
								ctg.index, 
								ctgStart, 
								nextContig.index, 
								newStart,bridge.getScore());						
						scaffolds[i].addFront(nextContig, bridge);
						nextContig.myVector = trialTrans;
						head[nextContig.index] = i;
						scaffoldExtended[nextContig.index] = true;
						
						extended = true;						
						
						scaffolds[i].view();						
						break;//for
					}else{
						System.out.printf(" NOT Extend %d from %d(%d) to %d(%d) with score %f\n",
								i,
								ctg.index, 
								ctgStart, 
								nextContig.index, 
								newStart,
								bridge.getScore());
					}
				}//for				
			}//while
			/*****************************************************************/
			System.out.printf("Finally, scaffold %d size %d  and is %s\n",
					i,
					scaffold.getLast().rightMost() - scaffold.getFirst().leftMost(),
					closed?"circular":"linear");			
		}//for
	}
	
}
