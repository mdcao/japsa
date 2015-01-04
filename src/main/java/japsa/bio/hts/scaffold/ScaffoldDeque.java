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

import java.util.ArrayDeque;
import java.util.Iterator;

/**
 * Implement scaffold as an array deque, that is a linear array that can be 
 * added/removed from either end
 * 
 * @author minhduc
 */
public final class ScaffoldDeque extends ArrayDeque<Contig>{

	private static final long serialVersionUID = -4310125261868862931L;	
	ArrayDeque<ContigBridge> bridges;//the bridges should syncronysed with contigs
	boolean closed = false;
	/**
	 * invariant: the direction of the decque is the same as the main (the longest one)
	 * @param myFContig
	 */

	public ScaffoldDeque(Contig myFContig){
		super();
		add(myFContig);//the first one		
		bridges = new ArrayDeque<ContigBridge>(); 
	}	

	/**
	 * Return 1 or -1 if the contig is at the first or last of the deque. 
	 * Otherwise, return 0
	 * @param ctg
	 * @return
	 */
	public int isEnd(Contig ctg){
		if (ctg == this.peekFirst())
			return 1;
		if (ctg == this.peekLast())
			return -1;

		return 0;			
	}

	public boolean isFirst(Contig ctg){
		return ctg == this.peekFirst();
	}

	public boolean isLast(Contig ctg){
		return ctg == this.peekLast();
	}

	/**
	 * Add a contig and its bridge to the beginning of the deque
	 * @param contig
	 * @param bridge
	 */
	public void addFront(Contig contig, ContigBridge bridge){
		this.addFirst(contig);
		bridges.addFirst(bridge);
	}

	/**
	 * Add a contig and its bridge to the end of the deque
	 * @param contig
	 * @param bridge
	 */
	public void addRear(Contig contig, ContigBridge bridge){
		this.addLast(contig);
		bridges.addLast(bridge);
	}

	public void combineScaffold(ScaffoldDeque scaffold, ContigBridge bridge, int myPos, int itsPos){
		//my pos == 1: add to front, else to rear
		//itsPos == 1: taken from front, else from rear
		
		Contig ctg = (itsPos==1) ? scaffold.removeFirst():scaffold.removeLast();			
		while(true){
			if (myPos == 1){ 
				addFront(ctg,bridge);
			}else{
				addRear(ctg,bridge);
			}				
			if (scaffold.isEmpty())
				break;

			ctg = (itsPos==1) ? scaffold.removeFirst():scaffold.removeLast();
			bridge = (itsPos==1) ? scaffold.bridges.removeFirst():scaffold.bridges.removeLast();				
		}		
	}		

	public ArrayDeque<ContigBridge> getBridges(){
		return bridges;
	}

	/**
	 * @param start the start to set
	 */
	
	public void view(){
		System.out.println("========================== START =============================");
		Iterator<ContigBridge> bridIter = bridges.iterator();

		for (Contig ctg:this){				
			System.out.printf("  contig %3d  ======" + (ctg.getRelDir() > 0?">":"<") + "%6d  %6d %s ",ctg.getIndex(), ctg.leftMost(),ctg.rightMost(), ctg.getName());
			if (bridIter.hasNext()){
				ContigBridge bridge = bridIter.next();
				System.out.printf("    %d\n", bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig));					
			}else
				System.out.println();			
		}
		System.out.println("============================ END ===========================");
	}

	public void viewSequence(){
		System.out.println("========================== START =============================");
		Iterator<ContigBridge> bridIter = bridges.iterator();
		for (Contig ctg:this){				
			System.out.printf("  contig %3d  ======" + (ctg.getRelDir() > 0?">":"<") + "%6d  %6d %s ",ctg.getIndex(), ctg.leftMost(),ctg.rightMost(), ctg.getName());

			if (bridIter.hasNext()){
				ContigBridge bridge = bridIter.next();
				System.out.printf("gaps =  %d\n", bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig));					
			}else
				System.out.println();
		}

		for (ContigBridge bridge:bridges){
			System.out.printf("Bridge %3d --> %3d using %d reads\n", bridge.firstContig.getIndex(),bridge.secondContig.getIndex(), bridge.getConnections().size());				
			int fill = bridge.fill();
			System.out.printf("Fill results %3d --> %3d with %d unfilled\n" ,bridge.firstContig.getIndex(),bridge.secondContig.getIndex(),fill);
		}
		System.out.println("============================ END ===========================");
	}
}