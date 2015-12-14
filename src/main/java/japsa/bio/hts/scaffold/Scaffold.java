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
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;

import java.io.IOException;
import java.util.LinkedList;
import java.util.Iterator;
import java.util.ListIterator;

/**
 * Implement scaffold as an array deque, that is a linear array that can be 
 * added/removed from either end
 * 
 * @author minhduc
 */
public final class Scaffold extends LinkedList<Contig>{
	ContigBridge closeBridge = null;//if not null, will bridge the last and the first contig
	ScaffoldVector circle = null;
	private static final long serialVersionUID = -4310125261868862931L;	
	LinkedList<ContigBridge> bridges;
	int scaffoldIndex;
	//boolean closed = false;
	/**
	 * invariant: the direction of the decque is the same as the main (the longest one)
	 * @param myFContig
	 */
	public Scaffold(int index){
		super();
		scaffoldIndex = index;
		bridges = new LinkedList<ContigBridge>(); 
	}
	public Scaffold(Contig myFContig){
		super();
		scaffoldIndex = myFContig.index;
		add(myFContig);//the first one		
		bridges = new LinkedList<ContigBridge>(); 
	}	


	public void setCloseBridge(ContigBridge bridge){
		assert bridge.firstContig.getIndex() == this.getLast().getIndex():"Closed bridge: " + bridge.hashKey + " <-> " +this.getLast().getIndex();
		closeBridge = bridge;
		circle = ScaffoldVector.composition(this.getLast().getVector(), ScaffoldVector.reverse(this.getFirst().getVector())); //first->last
		ScaffoldVector last2first = bridge.getTransVector();
		if(this.peekFirst().getIndex() == bridge.firstContig.getIndex())
			last2first = ScaffoldVector.reverse(last2first);
		circle = ScaffoldVector.composition(last2first, circle);
		//bridge.setContigScores();
		//closed = true;
	}

	/**
	 * Return 1 or -1 if the contig is at the first or last of the list. 
	 * Otherwise, return 0
	 * @param ctg
	 * @return
	 */
	public int isEnd(Contig ctg){
		if (ctg.getIndex() == this.peekLast().getIndex())
			return -1;
		if (ctg.getIndex() == this.peekFirst().getIndex())
			return 1;

		return 0;			
	}

	public boolean isFirst(Contig ctg){
		return ctg.getIndex() == this.peekFirst().getIndex();
	}

	public boolean isLast(Contig ctg){
		return ctg.getIndex() == this.peekLast().getIndex();
	}

	/**
	 * Add a contig and its bridge to the beginning of the deque
	 * @param contig
	 * @param bridge
	 */
	public void addFront(Contig contig, ContigBridge bridge){
		this.addFirst(contig);
		assert bridge.firstContig.getIndex() == contig.getIndex(): "Front: "+ bridge.hashKey + " " + contig.getIndex();
		bridges.addFirst(bridge);
		if(ScaffoldGraph.verbose)
			System.out.printf("...adding contig %d to scaffold %d backward!\n", contig.getIndex(), scaffoldIndex);
			
	}

	/**
	 * Add a contig and its bridge to the end of the deque
	 * @param contig
	 * @param bridge
	 */
	public void addRear(Contig contig, ContigBridge bridge){
		this.addLast(contig);
		assert bridge.secondContig.getIndex() == contig.getIndex():"Rear: "+ bridge.hashKey + " " + contig.getIndex();
		bridges.addLast(bridge);
		if(ScaffoldGraph.verbose)
			System.out.printf("...adding contig %d to scaffold %d forward!\n", contig.getIndex(), scaffoldIndex);
		
	}
	
	public Contig nearestMarker(Contig ctg, boolean forward){

		if(ScaffoldGraph.isRepeat(ctg)){
			if(ScaffoldGraph.verbose)
				System.out.println("Cannot determine nearest marker of a repeat!");
			return null;
		}
		int index = this.indexOf(ctg);
		if(index < 0) return null;
		ListIterator<Contig> iterator = this.listIterator(index);

		if(ScaffoldGraph.verbose){
			System.out.printf("Tracing scaffold %d from contig %d with index %d\n", scaffoldIndex, ctg.getIndex(), index);
			//this.view();
			System.out.printf("Finding nearest %s marker of contig %d:", forward?"next":"previous", ctg.getIndex());
		}
		Contig 	marker = null; 
		while((forward?iterator.hasNext():iterator.hasPrevious())){
			marker = (forward?iterator.next():iterator.previous());
			if(ScaffoldGraph.verbose)
				System.out.print("..."+marker.getIndex());
			if(marker != null && !ScaffoldGraph.isRepeat(marker) && marker.getIndex() != ctg.getIndex())
				break;
		}
		if(closeBridge!=null && (marker == null || ScaffoldGraph.isRepeat(marker))){
			marker = forward?this.getFirst():this.getLast();
			while((forward?iterator.hasNext():iterator.hasPrevious())){
				if(ScaffoldGraph.verbose)
					System.out.print("......"+marker.getIndex());
				if(marker != null && !ScaffoldGraph.isRepeat(marker) && marker.getIndex() != ctg.getIndex())
					break;
				else
					marker = (forward?iterator.next():iterator.previous());
			}
		}
		if(ScaffoldGraph.verbose)
			System.out.println();
		return marker;

	}
		
	public void trim(){
		if(ScaffoldGraph.verbose)
			System.out.println("Trimming scaffold: " + scaffoldIndex);
		if(closeBridge != null || this.isEmpty())
			return;
		//from right
		Contig rightmost = this.peekLast();
		while(ScaffoldGraph.isRepeat(rightmost) && this.size()>=1){
			if(ScaffoldGraph.verbose)
				System.out.println("...removing contig " + rightmost.getIndex());	
			this.removeLast();
			//bridges.removeLast();
			rightmost=this.peekLast();

		}
		if(this.size() <=1){
			bridges= new LinkedList<ContigBridge>();
			return;
		}
		
		while(!bridges.isEmpty()){
			if(bridges.peekLast().isContaining(rightmost))
				break;
			else{
				if(ScaffoldGraph.verbose)
					System.out.println("...removing bridge " + bridges.peekLast().hashKey);	
				bridges.removeLast();
			}
		}
		if(bridges.size() > 1){
			if(bridges.get(bridges.size() - 2).isContaining(rightmost)){
				if(ScaffoldGraph.verbose)
					System.out.println("...removing bridge " + bridges.peekLast().hashKey);	
				this.bridges.removeLast();
			}
		}
		
		//from left
		Contig leftmost = this.peekFirst();
		while(ScaffoldGraph.isRepeat(leftmost) && this.size()>=1){
			if(ScaffoldGraph.verbose)
				System.out.println("...removing contig " + leftmost.getIndex());			
			this.removeFirst();
			//bridges.removeFirst();
			leftmost=this.peekFirst();
		}
		if(this.size() <=1){
			bridges= new LinkedList<ContigBridge>();
			return;
		}
		
		while(!bridges.isEmpty()){
			if(bridges.peekFirst().isContaining(leftmost))
				break;
			else{
				if(ScaffoldGraph.verbose)
					System.out.println("...removing bridge " + bridges.peekFirst().hashKey);	
				bridges.removeFirst();
			}
		}
		if(bridges.size() > 1){
			if(bridges.get(1).isContaining(leftmost)){
				if(ScaffoldGraph.verbose)
					System.out.println("...removing bridge " + bridges.peekFirst().hashKey);	
				this.bridges.removeFirst();
			}
		}

	}
	public void setHead(int head){
		scaffoldIndex = head;
		for (Contig ctg:this)
			ctg.head = head;
	}
	public LinkedList<ContigBridge> getBridges(){
		return bridges;
	}

	/**
	 * @param start the start to set
	 */

	public void view(){
		System.out.println("========================== START =============================");
		Iterator<ContigBridge> bridIter = bridges.iterator();
		if(closeBridge!=null){
			System.out.println("Close bridge: " + closeBridge.hashKey + " Circularized vector: " + circle);
		}
		for (Contig ctg:this){				
			System.out.printf("  contig %3d  ======" + (ctg.getRelDir() > 0?">":"<") + "%6d  %6d %s ",ctg.getIndex(), ctg.leftMost(),ctg.rightMost(), ctg.getName());
			if (bridIter.hasNext()){
				ContigBridge bridge = bridIter.next();
				System.out.printf("    %d: %s\n", bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig), bridge.hashKey);					
			}else
				System.out.println();			
		}
		System.out.println("============================ END ===========================");
	}

	public void viewSequence(SequenceOutputStream out) throws IOException{		


		System.out.println("========================== START =============================");
		if(closeBridge!=null){
			System.out.println("Close bridge: " + closeBridge.hashKey + " Circularized vector: " + circle);
		}
		Iterator<ContigBridge> bridIter = bridges.iterator();
		for (Contig ctg:this){				
			System.out.printf("  contig %3d  ======" + (ctg.getRelDir() > 0?">":"<") + "%6d  %6d %s ",ctg.getIndex(), ctg.leftMost(),ctg.rightMost(), ctg.getName());

			if (bridIter.hasNext()){
				ContigBridge bridge = bridIter.next();
				System.out.printf("gaps =  %d\n", bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig));					
			}else
				System.out.println();
		}

		System.out.println("Size = " + size() + " sequence");
		
		SequenceBuilder seq = new SequenceBuilder(Alphabet.DNA16(), 1024*1024,  "Scaffold" + scaffoldIndex);
		JapsaAnnotation anno = new JapsaAnnotation();


		ContigBridge.Connection bestCloseConnection = null;				
		Contig leftContig, rightContig;		
/*
 * Nanopore reads:
 * 	====================================                ==========================================
 *                           |      |                     |       |                |       |
 *                           |      |                     |       |                |       |
 *        <fillFrom>         |      |     leftContig      |       |   <fillFrom>   |       |    rightContig
 * Contigs:   ...		~~~~~*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*~~~	...		~~~*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
 * 				          startLeft                            endLeft
 * 
 *                                that's what happens below!
 * 
 */
		
		rightContig = getFirst();		

		//startLeft: the leftPoint of leftContig, endLeft: rightPoint of left Contig
		int startLeft = (rightContig.getRelDir() > 0)?1:rightContig.length(); //starting point after the last fillFrom
		int endLeft   = (rightContig.getRelDir() < 0)?1:rightContig.length();
		
		if (closeBridge != null){
			bestCloseConnection = closeBridge.fewestGapConnection();
			leftContig = closeBridge.firstContig;
			startLeft = bestCloseConnection.filling(null, null); //adjust the starting point

			anno.addDescription("Circular");

		}else
			anno.addDescription("Linear");


		bridIter = bridges.iterator();
		Iterator<Contig> ctgIter = this.iterator();
		leftContig = ctgIter.next();//The first

		for (ContigBridge bridge:bridges){
			//System.out.println("------------------------------------ START ------------------------------------");
			rightContig = ctgIter.next();

			ContigBridge.Connection connection = bridge.fewestGapConnection();

			endLeft = (leftContig.getRelDir()>0)?(connection.firstAlignment.refEnd):
				(connection.firstAlignment.refStart);

			if (startLeft<endLeft){
				JapsaFeature feature = 
						new JapsaFeature(seq.length() + 1, seq.length() + endLeft - startLeft + 1,
								"CONTIG",leftContig.getName(),'+',"");
				feature.addDesc(leftContig.getName() + "+("+startLeft +"," + endLeft+")");
				anno.add(feature);				
				seq.append(leftContig.contigSequence.subSequence(startLeft - 1, endLeft));

				leftContig.portionUsed += (1.0 + endLeft - startLeft + 1) / leftContig.length();
			}else{

				JapsaFeature feature = 
						new JapsaFeature(seq.length() + 1, seq.length() + startLeft - endLeft + 1,
								"CONTIG",leftContig.getName(),'+',"");
				feature.addDesc(leftContig.getName() + "-("+endLeft +"," + startLeft+")");
				anno.add(feature);

				seq.append(Alphabet.DNA.complement(leftContig.contigSequence.subSequence(endLeft - 1, startLeft)));
				leftContig.portionUsed += (1.0 - endLeft + startLeft + 1) / leftContig.length();
			}			
			//Fill in the connection
			startLeft = connection.filling(seq, anno);					
			leftContig = rightContig;			
					
		}//for

		//leftContig = lastContig in the queue
		if (bestCloseConnection != null){
			endLeft = (leftContig.getRelDir()>0)?(bestCloseConnection.firstAlignment.refEnd):
				(bestCloseConnection.firstAlignment.refStart);	
		}
		else
			endLeft = (rightContig.getRelDir() < 0)?1:rightContig.length();
		
		if (startLeft<endLeft){
			JapsaFeature feature = 
					new JapsaFeature(seq.length() + 1, seq.length() + endLeft - startLeft,
							"CONTIG",leftContig.getName(),'+',"");
			feature.addDesc(leftContig.getName() + "+("+startLeft +"," + endLeft+")");
			anno.add(feature);				
			seq.append(leftContig.contigSequence.subSequence(startLeft - 1, endLeft));
			leftContig.portionUsed += (1.0 + endLeft - startLeft + 1) / leftContig.length();
		}else{

			JapsaFeature feature = 
					new JapsaFeature(seq.length() + 1, seq.length() + startLeft - endLeft,
							"CONTIG",leftContig.getName(),'+',"");
			feature.addDesc(leftContig.getName() + "-("+endLeft +"," + startLeft+")");
			anno.add(feature);

			seq.append(Alphabet.DNA.complement(leftContig.contigSequence.subSequence(endLeft - 1, startLeft)));
			leftContig.portionUsed += (1.0 - endLeft + startLeft + 1) / leftContig.length();
		}
		if (bestCloseConnection != null){	
			System.out.printf("Append bridge %d -- %d\n",closeBridge.firstContig.index,  closeBridge.secondContig.index);
			bestCloseConnection.filling(seq, anno);	
		}

		System.out.println("============================ END ===========================");
		//JapsaAnnotation.write(seq.toSequence(), anno, out); //uncomment this line and comment next line for debug
		seq.writeFasta(out);
	}
}