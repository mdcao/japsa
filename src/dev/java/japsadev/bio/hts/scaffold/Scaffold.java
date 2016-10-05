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

package japsadev.bio.hts.scaffold;

import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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
	int len=-1;
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
		assert bridge.firstContig.getIndex() == contig.getIndex(): "Front prob: "+ bridge.hashKey + " not connect " + contig.getIndex() + " and " + this.getFirst().getIndex();
		this.addFirst(contig);
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
		assert bridge.secondContig.getIndex() == contig.getIndex():"Rear prob: "+ bridge.hashKey + " not connect " + this.getLast().getIndex() + " and " + contig.getIndex();
		this.addLast(contig);
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
	//reset prevScore or nextScore to 0 according to removed bridges.
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
				//bridges.peekLast();
				bridges.removeLast().resetContigScores();
			}
		}
		if(bridges.size() > 1){
			if(bridges.get(bridges.size() - 2).isContaining(rightmost)){
				if(ScaffoldGraph.verbose)
					System.out.println("...removing bridge " + bridges.peekLast().hashKey);	
				bridges.removeLast().resetContigScores();
			}
		}
		
		//from left
		Contig leftmost = this.peekFirst();
		while(ScaffoldGraph.isRepeat(leftmost) && this.size()>=1){
			if(ScaffoldGraph.verbose)
				System.out.println("...removing contig " + leftmost.getIndex());			
			this.removeFirst();
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
				//bridges.peekFirst();
				bridges.removeFirst().resetContigScores();
			}
		}
		if(bridges.size() > 1){
			if(bridges.get(1).isContaining(leftmost)){
				if(ScaffoldGraph.verbose)
					System.out.println("...removing bridge " + bridges.peekFirst().hashKey);
				bridges.removeFirst().resetContigScores();
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
	
	/**
	 * Return the length of this scaffold
	 * Check out quast (https://github.com/ablab/quast)
	 */
	public int length(){
		if(isEmpty())
			return 0;
		if(len > 0)
			return len;
		int len = getLast().rightMost() - getFirst().leftMost();
		if(circle!=null)
			len = Math.abs(circle.getMagnitute());
		return len;
	}

	public synchronized void viewSequence(SequenceOutputStream fout, SequenceOutputStream jout) throws IOException{		


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
		
		// Synchronize positions of 2 contigs (myVector) of a bridge based on the real list of contigs
		// TODO: do the same with viewAnnotation()
		assert this.size()==bridges.size()+1:"Number of contigs ("+this.size()+")" + " doesn't agree with number of bridges ("+bridges.size()+"!";
		for(int i=0;i<bridges.size();i++){
			ContigBridge brg=bridges.get(i).clone(this.get(i), this.get(i+1));
			bridges.set(i, brg);
		}
		if(closeBridge!=null)
			closeBridge=closeBridge.clone(this.get(this.size()-1), this.get(0));
		
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
 * 				          startLeft(1)                          endLeft(1)       startLeft(2)                 ...
 * 
 *                                that's what happens below!
 * 
 */
		
		rightContig = getFirst();		

		//startLeft: the leftPoint of leftContig, endLeft: rightPoint of left Contig
		int startLeft = (rightContig.getRelDir() > 0)?0:rightContig.length(); //starting point after the last fillFrom
		int endLeft   = (rightContig.getRelDir() < 0)?0:rightContig.length();
		//TODO: re-check the coordinate of two ends (0 or 1, inclusive/exclusive)
//		int startLeft = (rightContig.getRelDir() > 0)?1:rightContig.length(); //starting point after the last fillFrom
//		int endLeft   = (rightContig.getRelDir() < 0)?1:rightContig.length();
		
		
		/* uncomment for illumina-based */
		//int closeDis = 0;
		
		if (closeBridge != null){
			bestCloseConnection = closeBridge.fewestGapConnection();
			leftContig = closeBridge.firstContig;
			/* uncomment for longread-based */
			startLeft = bestCloseConnection.filling(null, null); 
			/* uncomment for illumina-based */
//			closeDis =closeBridge.getTransVector().distance(closeBridge.firstContig, closeBridge.secondContig);		
//			if(closeDis > 0)
//				startLeft = bestCloseConnection.filling(null, null); //adjust the starting point

			anno.addDescription("Circular");

		}else
			anno.addDescription("Linear");


		Iterator<Contig> ctgIter = this.iterator();
		leftContig = ctgIter.next();//The first

		for (ContigBridge bridge:bridges){
			rightContig = ctgIter.next();
			ContigBridge.Connection connection = bridge.fewestGapConnection();
			/* uncomment for longread-based */
			endLeft = 	(leftContig.getRelDir()>0)?(connection.firstAlignment.refEnd):
				 		(connection.firstAlignment.refStart);
			/* uncomment for illumina-based */
//			int distance = bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig);
//			if(distance < 0){
//				endLeft = (leftContig.getRelDir()>0)?(leftContig.length()-Math.abs(distance)):Math.abs(distance);
//			}else{
//				endLeft = (leftContig.getRelDir()>0)?(connection.firstAlignment.refEnd):
//					(connection.firstAlignment.refStart);
//			}
			
			if (startLeft<endLeft){
				JapsaFeature feature = 
						new JapsaFeature(seq.length() + 1, seq.length() + endLeft - startLeft,
								"CONTIG",leftContig.getName(),'+',"");
				feature.addDesc(leftContig.getName() + "+["+startLeft +"," + endLeft+")");
				anno.add(feature);		
				if(ScaffoldGraph.verbose)
					System.out.println("Append " + leftContig.getName() + ": " + startLeft + " to " + (endLeft-1));
				
				seq.append(leftContig.contigSequence.subSequence(startLeft, endLeft));
				//seq.append(leftContig.contigSequence.subSequence(startLeft-1, endLeft));

			}else{

				JapsaFeature feature = 
						new JapsaFeature(seq.length() + 1, seq.length() + startLeft - endLeft,
								"CONTIG",leftContig.getName(),'+',"");
				feature.addDesc(leftContig.getName() + "-["+endLeft +"," + startLeft+")");
				anno.add(feature);
				
				if(ScaffoldGraph.verbose)
					System.out.println("Append RC of " + leftContig.getName() + ": " + endLeft + " to " + (startLeft-1));

				seq.append(Alphabet.DNA.complement(leftContig.contigSequence.subSequence(endLeft, startLeft)));
				//seq.append(Alphabet.DNA.complement(leftContig.contigSequence.subSequence(endLeft-1, startLeft)));
			}			
			/* uncomment for longread-based */
			startLeft = connection.filling(seq, anno);
			/* uncomment for illumina-based */
//			if(distance < 0){
//				startLeft = (rightContig.getRelDir() > 0)?1:rightContig.length();
//			}else{
//				//Fill in the connection
//				startLeft = connection.filling(seq, anno);	
//			}
			leftContig = rightContig;			
					
		}//for

		//leftContig = lastContig in the queue
		if (bestCloseConnection != null){
			/* uncomment for longread-based */
			endLeft = (leftContig.getRelDir()>0)?(bestCloseConnection.firstAlignment.refEnd):
				(bestCloseConnection.firstAlignment.refStart);
			/* uncomment for illumina-based */
//			if(closeDis > 0)
//				endLeft = (leftContig.getRelDir()>0)?(bestCloseConnection.firstAlignment.refEnd):
//				(bestCloseConnection.firstAlignment.refStart);	
//			else{ 
//				endLeft = (rightContig.getRelDir() < 0)?Math.abs(closeDis):rightContig.length()-Math.abs(closeDis);
//			}
		}
		else
			endLeft = (rightContig.getRelDir() < 0)?1:rightContig.length();
		
		if (startLeft<endLeft){
			JapsaFeature feature = 
					new JapsaFeature(seq.length() + 1, seq.length() + endLeft - startLeft,
							"CONTIG",leftContig.getName(),'+',"");
			feature.addDesc(leftContig.getName() + "+["+startLeft +"," + endLeft+")");
			anno.add(feature);				
			if(ScaffoldGraph.verbose)
				System.out.println("Append " + leftContig.getName() + ": " + startLeft + " to " + (endLeft-1));

			seq.append(leftContig.contigSequence.subSequence(startLeft, endLeft));
			//seq.append(leftContig.contigSequence.subSequence(startLeft - 1, endLeft));
		}else{

			JapsaFeature feature = 
					new JapsaFeature(seq.length() + 1, seq.length() + startLeft - endLeft,
							"CONTIG",leftContig.getName(),'+',"");
			feature.addDesc(leftContig.getName() + "-["+endLeft +"," + startLeft+")");
			anno.add(feature);
			
			if(ScaffoldGraph.verbose)
				System.out.println("Append RC of " + leftContig.getName() + ": " + endLeft + " to " + (startLeft-1));
			seq.append(Alphabet.DNA.complement(leftContig.contigSequence.subSequence(endLeft, startLeft)));
			//seq.append(Alphabet.DNA.complement(leftContig.contigSequence.subSequence(endLeft - 1, startLeft)));
		}
		if (bestCloseConnection != null){	
			System.out.printf("Append bridge %d -- %d\n",closeBridge.firstContig.index,  closeBridge.secondContig.index);
			/* uncomment for longread-based */
			bestCloseConnection.filling(seq, anno);	
			/* uncomment for illumina-based */
//			if(closeDis >0 )
//				bestCloseConnection.filling(seq, anno);	
		}

		System.out.println("============================ END ===========================");
		len = seq.length();
		JapsaAnnotation.write(seq.toSequence(), anno, jout); 
		seq.writeFasta(fout);
	}
	/* Output annotation of this scaffold
	 * TODO: output annotations from the filling sequences
	 */
	public synchronized void viewAnnotation(SequenceOutputStream out) throws IOException{		

		SequenceBuilder seq = new SequenceBuilder(Alphabet.DNA16(), 1024*1024,  "Scaffold" + scaffoldIndex);
		JapsaAnnotation anno = new JapsaAnnotation();

		ContigBridge.Connection bestCloseConnection = null;				
		Contig leftContig, rightContig;		
/*
 * Nanopore reads:
 * 	====================================                ==========================================
 *                           |      |                     |       |                |       |
 *                           |      |                     |       |                |       |
 *         <filling>         |      |     leftContig      |       |    <filling>   |       |    rightContig
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
			startLeft = bestCloseConnection.filling(null, null); 

			anno.addDescription("Circular");

		}else
			anno.addDescription("Linear");

		Iterator<Contig> ctgIter = this.iterator();
		leftContig = ctgIter.next();//The first

		JapsaFeature lastGene=null; 
		if(leftContig.resistanceGenes.size()>0)
			lastGene = startLeft<endLeft?leftContig.resistanceGenes.get(0):leftContig.resistanceGenes.get(leftContig.resistanceGenes.size()-1);

		for (ContigBridge bridge:bridges){
			rightContig = ctgIter.next();
			ContigBridge.Connection connection = bridge.fewestGapConnection();

			endLeft = 	(leftContig.getRelDir()>0)?(connection.firstAlignment.refEnd):
				 		(connection.firstAlignment.refStart);
			
			/**********************************************************************************************/
			ArrayList<JapsaFeature> resistLeft = leftContig.getFeatures(leftContig.resistanceGenes, startLeft, endLeft),
					insertLeft = leftContig.getFeatures(leftContig.insertSeq, startLeft, endLeft),
					oriLeft = leftContig.getFeatures(leftContig.oriRep, startLeft, endLeft);

			for(JapsaFeature resist:resistLeft){
				resist.setStart(resist.getStart()+seq.length());
				resist.setEnd(resist.getEnd()+seq.length());
				anno.add(resist);
			}
			Collections.sort(resistLeft);
			if(resistLeft.size() > 0){
				JapsaFeature leftGene = startLeft<endLeft?resistLeft.get(0):resistLeft.get(resistLeft.size()-1),
							rightGene = startLeft>endLeft?resistLeft.get(0):resistLeft.get(resistLeft.size()-1);
				//extract first and last gene here
				if(lastGene!=null)
					System.out.println(lastGene.getID() + "...next to..." + leftGene.getID());
				lastGene = rightGene;
			}
			for(JapsaFeature insert:insertLeft){
				insert.setStart(insert.getStart()+seq.length());
				insert.setEnd(insert.getEnd()+seq.length());			
				anno.add(insert);
			}
			for(JapsaFeature ori:oriLeft){
				ori.setStart(ori.getStart()+seq.length());
				ori.setEnd(ori.getEnd()+seq.length());
				anno.add(ori);
			}
			/**********************************************************************************************/

			if (startLeft<endLeft)			
				seq.append(leftContig.contigSequence.subSequence(startLeft, endLeft));
			else
				seq.append(Alphabet.DNA.complement(leftContig.contigSequence.subSequence(endLeft, startLeft)));

			startLeft = connection.filling(seq, new JapsaAnnotation());
			leftContig = rightContig;			
					
		}//for

		if (bestCloseConnection != null){
			endLeft = (leftContig.getRelDir()>0)?(bestCloseConnection.firstAlignment.refEnd):
				(bestCloseConnection.firstAlignment.refStart);

		}
		else
			endLeft = (rightContig.getRelDir() < 0)?1:rightContig.length();
		
		/**********************************************************************************************/
		ArrayList<JapsaFeature> resistLeft = leftContig.getFeatures(leftContig.resistanceGenes, startLeft, endLeft),
				insertLeft = leftContig.getFeatures(leftContig.insertSeq, startLeft, endLeft),
				oriLeft = leftContig.getFeatures(leftContig.oriRep, startLeft, endLeft);

		for(JapsaFeature resist:resistLeft){
			resist.setStart(resist.getStart()+seq.length());
			resist.setEnd(resist.getEnd()+seq.length());
			anno.add(resist);
		}
		if(resistLeft.size() > 0){
			JapsaFeature leftGene = startLeft<endLeft?resistLeft.get(0):resistLeft.get(resistLeft.size()-1);
			if(lastGene!=null)
				System.out.println(lastGene.getID() + "...next to..." + leftGene.getID());
		}
		for(JapsaFeature insert:insertLeft){
			insert.setStart(insert.getStart()+seq.length());
			insert.setEnd(insert.getEnd()+seq.length());			
			anno.add(insert);
		}
		for(JapsaFeature ori:oriLeft){
			ori.setStart(ori.getStart()+seq.length());
			ori.setEnd(ori.getEnd()+seq.length());
			anno.add(ori);
		}
		/**********************************************************************************************/

		if (startLeft<endLeft)			
			seq.append(leftContig.contigSequence.subSequence(startLeft, endLeft));

		else
			seq.append(Alphabet.DNA.complement(leftContig.contigSequence.subSequence(endLeft, startLeft)));
			
		if (bestCloseConnection != null)	
			bestCloseConnection.filling(seq, new JapsaAnnotation());	
		
		anno.setSequence(seq.toSequence());
		JapsaAnnotation.write(null, anno, out); 
		len = seq.length();
	}
}
