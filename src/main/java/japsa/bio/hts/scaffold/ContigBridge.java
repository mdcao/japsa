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
 * 19/12/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.hts.scaffold;



import htsjdk.samtools.CigarElement;
import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.util.Pair;


/**
 * Create a bridge that connects two contigs. The bridged can be ranked based
 * on the confidence so that more confident bridge is used before.
 * Note that two contigs can have more than one bridges from circular 
 * sequence or false positives. 
 * @author minhduc
 *
 */

public class ContigBridge implements Comparable<ContigBridge>{

	final Contig firstContig, secondContig;
	final String hashKey;
	final int orderIndex;

	private double score = 0;//more is better
	private ScaffoldVector transVector = null;
	private ArrayList<Connection> connections;//a list of connections that make up this

	public ContigBridge(Contig c1, Contig c2, int ind){
		firstContig = c1;
		secondContig = c2;
		orderIndex = ind;		
		hashKey = makeHash(c1.index,c2.index, orderIndex);

		connections = new  ArrayList<Connection>();
		//TODO: if firstContig = secondContig: circular or not
	}

	public static String makeHash(int aIndex, int bIndex, int order){		
		return aIndex+"#"+bIndex + "#" + order;
	}

	public boolean consistentWith(ScaffoldVector aVector){
		return (aVector.direction == transVector.direction)
				&& (aVector.magnitude * 1.0 / transVector.magnitude > 0.9)
				&& (aVector.magnitude * 1.0 / transVector.magnitude < 1.1)
				;
	}

	public double addConnection(ReadFilling readSequence, 
			AlignmentRecord firstAlignment, 
			AlignmentRecord secondAlignment, 
			ScaffoldVector trans, 
			double sc){

		//NB: firstAlignment for firstContig, secondAlignment for secondContig

		if (transVector == null){
			transVector = trans;
			score = sc;			
			connections.add(new Connection(readSequence, firstAlignment,secondAlignment,trans));			
		}else{
			boolean alreadyIn = false;
			for (int i = 0; i < connections.size();i++){
				Connection oldConnect = connections.get(i); 
				if (connections.get(i).readID / 3 == firstAlignment.readID /3){
					if (sc > oldConnect.score){
						//replace only if better						
						Connection newConnect = new Connection(readSequence, firstAlignment,secondAlignment,trans);
						connections.set(i, newConnect);
						transVector.magnitude = (transVector.magnitude * connections.size() - oldConnect.trans.magnitude + trans.magnitude)/connections.size();
					}
					alreadyIn = true;
					break;
				}//if
			}//for
			if (alreadyIn)
				return score;

			//if not already in			
			Connection newConnect = new Connection(readSequence, firstAlignment,secondAlignment,trans);
			connections.add(newConnect);			
			transVector.magnitude = (transVector.magnitude * connections.size() + trans.magnitude) / (connections.size() + 1);

			score += sc;
		}
		return score;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	public void setScore(double s) {
		score = s;
	}

	/**
	 * @return the transVector
	 */
	public ScaffoldVector getTransVector() {
		return transVector;
	}
	/**
	 * @return the connections
	 */
	public ArrayList<Connection> getConnections() {
		return connections;
	}	
	/**
	 * Assume connections has been sorted
	 * @param contig
	 * @return
	 */
	private HashMap<Contig, ScaffoldVector> fillFrom(Contig fromContig, HashMap<Contig, Range> contigRanges){
		HashMap<Contig, ScaffoldVector> mainVectors = new HashMap<Contig, ScaffoldVector>();
		mainVectors.put(fromContig, new ScaffoldVector());

		Contig toContig = null;
		if (fromContig == firstContig)
			toContig = secondContig;
		else if (fromContig == secondContig)
			toContig = firstContig;
		else
			throw new RuntimeException("Need to pass in a contig part of the bridge");


		for (Connection connection:connections){
			ReadFilling readFilling = connection.read;
			readFilling.sortAlignment();
			System.out.printf("  Fill using read %s to connect %d (%b) to %d(%b)\n", 
					readFilling.readSequence.getName(),
					connection.firstAlignment.contig.index, 
					connection.firstAlignment.strand,
					connection.secondAlignment.contig.index,
					connection.secondAlignment.strand);

			int fromStart = 0, toStart = 0;
			if (fromContig == firstContig){
				fromStart = connection.firstAlignment.readAlignmentStart();
				toStart = connection.secondAlignment.readAlignmentStart();				
			}else{
				//from is the second
				toStart = connection.firstAlignment.readAlignmentStart();
				fromStart  = connection.secondAlignment.readAlignmentStart();
			}

			if (toStart < fromStart){
				System.out.printf("  -----> Reverse read\n");
				readFilling = readFilling.reverse();
				readFilling.sortAlignment();
			}
			int startIndex = -1, endIndex = -1; 
			for (int i = 0; i < readFilling.alignments.size();i++){
				if (readFilling.alignments.get(i).contig == fromContig)
					startIndex = i;

				if (readFilling.alignments.get(i).contig == toContig)
					endIndex = i;
			}

			if (startIndex < 0 || endIndex < 0)
				throw new RuntimeException("Error in fillFrom()");

			if (startIndex >= endIndex){
				System.out.println("The two overlap, good news, but need to do some thing " + startIndex + " vs " + endIndex);
				continue;
			}
			//assert startIndex < endIndex

			//Prim algorithm			
			for (int i = 0; i < readFilling.alignments.size();i++){
				AlignmentRecord inAlignment = readFilling.alignments.get(i);
				Range range = contigRanges.get(inAlignment.contig);
				if (range == null){
					range = new Range(inAlignment.refStart,inAlignment.refEnd);
					contigRanges.put(inAlignment.contig, range);
				}else{
					if (inAlignment.refStart < range.start )
						range.start = inAlignment.refStart;

					if (inAlignment.refEnd > range.end )
						range.end = inAlignment.refEnd;
				}
			}//for i


			for (int i = 0; i < readFilling.alignments.size();i++){
				AlignmentRecord inAlignment = readFilling.alignments.get(i);
				ScaffoldVector inVector = mainVectors.get(inAlignment.contig);
				//skip if this alignment has not been added
				if (inVector == null)			
					continue;


				//for simplicity, try the righ-most position
				int inRefPos = inAlignment.strand? inAlignment.refEnd:inAlignment.refStart;
				int inReadPos = inAlignment.strand? inAlignment.readEnd:inAlignment.readStart;//essentially, it = readAlignmentEnd().

				for (int j = i + 1; j < readFilling.alignments.size();j++){
					AlignmentRecord outAlignment = readFilling.alignments.get(j);
					System.out.printf("     Tempting [%d %d] (%d %d) on contig %d(%b) score %d vs [%d %d] (%d %d) on contig %d(%b) score %d\n",
							inAlignment.refStart,
							inAlignment.refEnd,
							inAlignment.readStart,
							inAlignment.readEnd,
							inAlignment.contig.index,
							inAlignment.strand,
							inAlignment.score,
							outAlignment.refStart,
							outAlignment.refEnd,
							outAlignment.readStart,
							outAlignment.readEnd,
							outAlignment.contig.index,
							outAlignment.strand,
							outAlignment.score
							);

					if (inReadPos <  outAlignment.readAlignmentStart())
						break;//no point going further
					if (inReadPos >  outAlignment.readAlignmentEnd()){
						//TODO: may still get some information
						continue;
					}

					int outRefPos = positionOnRef(inReadPos, outAlignment);
					System.out.printf("         Found %d on %d matches with %d on %d through %d on %s\n",
							outRefPos,
							outAlignment.contig.index,
							inRefPos,
							inAlignment.contig.index,
							inReadPos,
							readFilling.readSequence.getName()
							);


					//outRefPos == inRefPos
					if (outRefPos <=0){
						throw new RuntimeException("Error code = 100");
					}
					ScaffoldVector newVector = new ScaffoldVector();
					if (inAlignment.strand == outAlignment.strand){
						newVector.direction = 1;
						newVector.magnitude = inRefPos - outRefPos;
					}else{
						newVector.direction = -1;
						newVector.magnitude = inRefPos + outRefPos;
					}	
					newVector = ScaffoldVector.composition(newVector, inVector);
					ScaffoldVector outVector = mainVectors.get(outAlignment.contig);
					if (outVector == null){
						mainVectors.put(outAlignment.contig, newVector);
						System.out.printf("         Putting %d (%b) with vector (%s)\n",
								outAlignment.contig.index,
								outAlignment.strand,
								newVector.toString()
								);
					}else{//check for consistency
						if ((outVector.direction != newVector.direction))
							System.out.printf("         Fatal inconsistency direction\n");
						else if (Math.abs(outVector.getMagnitute() - newVector.getMagnitute()) > 10){
							System.out.printf("         Fatal inconsistency magnitute %d %d \n", outVector.getMagnitute(), newVector.getMagnitute());
						}else{
							System.out.printf("         In-Out consistent  %d %d \n",outVector.getMagnitute(), newVector.getMagnitute());
						}
					}
				}//for
			}
		}
		return mainVectors;		
	}

	public Connection fewestGapConnection(){
		Collections.sort(connections);

		//Find the best connections (has fewest gaps)
		Connection gapsBestConnection = null;
		int gapsBest = Integer.MAX_VALUE;

		for (Connection connection:connections){
			int gapsBt = connection.gapsBetween();
			if (gapsBt < gapsBest){
				gapsBest = gapsBt;
				gapsBestConnection = connection;
			}
		}	
		//System.out.println("                  Min gaps between = " + gapsBest + " "  + gapsBestConnection.read.getReadSequence().getName());
		return gapsBestConnection;

	}


	public void fillConnection(Contig left, Contig right){
		Collections.sort(connections);
		int gaps = Integer.MAX_VALUE;
		Pair<Integer, ArrayList<FillRecord>> bestPair = null;

		Connection gapsBestConnection = null;
		int gapsBest = Integer.MAX_VALUE;


		for (Connection connection:connections){
			int gapsBt = connection.gapsBetween();
			if (gapsBt < gapsBest){
				gapsBest = gapsBt;
				gapsBestConnection = connection;
			}

			Pair<Integer, ArrayList<FillRecord>>  pair = connection.fillConnection(left);
			if (pair.getKey() == 0){
				bestPair = pair;
				gaps = 0;
				break;

			}
			if (pair.getKey() <= gaps){
				gaps = pair.getKey(); 
				bestPair = pair;
			}
		}
		System.out.println("                  Min gaps = " + gaps);

		for (FillRecord record:bestPair.getValue()){
			System.out.println("                   Record " + record.contigSequence.getName() + ":" + record.start + "-" + record.end);
		}

		System.out.println("                  Min gaps between = " + gapsBest + " "  + gapsBestConnection.read.getReadSequence().getName());
	}

	/**
	 * Fill in the gap between contig a and contig b, note that they are the
	 * first and second contigs but not neccesarily in that order
	 * 
	 * @param a
	 * @param b
	 * @param startLeft: 1-index position
	 * @return a pair of numbers, indicating the rightmost position of the left contig, and
	 * the left-most position of the right contig
	 */
	public void fillGap(Contig left, Contig right){
		System.out.printf("\n\n"
				+ "Contig %3d (%d) -> Contig %3d (%d) Vector (%s) score = %f distance = %d\n",
				this.firstContig.index,
				this.firstContig.length(),
				this.secondContig.index,
				this.secondContig.length(),
				transVector.toString(),
				this.score,
				transVector.distance(firstContig, secondContig)				
				);		

		Collections.sort(connections);
		//int leftContigEnd = 0, rightContigStart = 0;

		/**
		 * Algorithm:
		 * For each connection, look at every pair-wise inter-connected alignments, 
		 * and try to find the relative position between the two contigs of the two alignment
		 * 
		 *  Alignments are sorted in the direction of reads, so fill in the direction of read for each read. In the end it doesnt matter, does it?
		 */


		HashMap<Contig, Range> contigRanges  = new HashMap<Contig, Range>();

		//Map the vector to the contigs that are connect to first contig		
		HashMap<Contig, ScaffoldVector> leftVectors = this.fillFrom(left, contigRanges);		
		for (Contig contig:leftVectors.keySet()){
			Range range = contigRanges.get(contig);
			ScaffoldVector vector = leftVectors.get(contig);
			System.out.printf("               Contig %d (%d,%d) of %d and vector (%s)\n",contig.index, range.start, range.end, contig.length(),vector.toString());
		}		
		if (leftVectors.containsKey(right)){
			System.out.printf("             CONNECTED %d %d\n",left.index, right.index);

			//Range range = contigRanges.get(left);
			//leftContigEnd = (left.getRelDir()>0) ? range.end:range.start;
			//range = contigRanges.get(right);
			//rightContigStart = (right.getRelDir()>0) ? range.start:range.end;
			//TODO: This is NOT good			


			//return new Pair<Integer, Integer>(leftContigEnd, rightContigStart);			
		}else{			
			HashMap<Contig, ScaffoldVector> rightVectors = this.fillFrom(right, contigRanges);		
			for (Contig contig:rightVectors.keySet()){
				Range range = contigRanges.get(contig);
				ScaffoldVector vector = rightVectors.get(contig);
				System.out.printf("               Contig %d (%d,%d) of %d and vector (%s)\n",contig.index, range.start, range.end, contig.length(),vector.toString());
			}	
			if (!rightVectors.containsKey(left)){
				System.out.printf("             DISCONNECTED %d %d\n",left.index, right.index);
			}			

			int smallestGaps = Integer.MAX_VALUE;
			Connection bestConnection = null;

			for (Connection connection:connections){
				ReadFilling read = connection.read;
				//assert. read.alignments sorted
				System.out.printf("     Read %s\n",read.readSequence.getName());


				int leftStart = 0, rightStart = 0;
				AlignmentRecord rightMostOfLeft, leftMostOfRight;
				if (left == firstContig){
					leftStart = connection.firstAlignment.readAlignmentStart();
					rightStart = connection.secondAlignment.readAlignmentStart();				
					rightMostOfLeft  = connection.firstAlignment;
					leftMostOfRight  = connection.secondAlignment;					
				}else{
					//from is the second
					rightStart = connection.firstAlignment.readAlignmentStart();
					leftStart  = connection.secondAlignment.readAlignmentStart();
					leftMostOfRight  = connection.firstAlignment;
					rightMostOfLeft= connection.secondAlignment;
				}
				//leftStart	

				for (AlignmentRecord record:read.alignments){
					Range range = contigRanges.get(record.contig);

					if (leftVectors.containsKey(record.contig)){
						if (leftStart < rightStart){//left is on the low, so get the highest
							if (record.readAlignmentEnd() > rightMostOfLeft.readAlignmentEnd())
								rightMostOfLeft = record;
						}else{
							if (record.readAlignmentStart() < rightMostOfLeft.readAlignmentStart()){
								rightMostOfLeft = record;
							}
						}
					}

					if (rightVectors.containsKey(record.contig)){
						if (leftStart > rightStart){//left is on the low, so get the highest
							if (record.readAlignmentEnd() > leftMostOfRight.readAlignmentEnd())
								leftMostOfRight = record;
						}else{
							if (record.readAlignmentStart() < leftMostOfRight.readAlignmentStart()){
								leftMostOfRight = record;
							}
						}
					}

					System.out.printf("       [%6d %6d] Contig %4d [%7d %7d] [%7d %7d] %s %s %s\n",							
							record.readStart,
							record.readEnd,
							record.contig.index,
							record.refStart,
							record.refEnd,
							range.start,
							range.end,							
							leftVectors.containsKey(record.contig)?"L":" ",
									rightVectors.containsKey(record.contig)?"R":" ",
											(leftVectors.containsKey(record.contig) && rightVectors.containsKey(record.contig))?"B":" "		
							);
				}//for
				int gaps = Math.min(Math.abs(rightMostOfLeft.readAlignmentEnd() -leftMostOfRight.readAlignmentStart()),
						Math.abs(leftMostOfRight.readAlignmentEnd() -rightMostOfLeft.readAlignmentStart()) ); 
				System.out.printf("          Gaps %d  : [%6d %6d] Contig %4d [%7d %7d] to [%6d %6d] Contig %4d [%7d %7d]\n",	
						gaps,								
						rightMostOfLeft.readStart,
						rightMostOfLeft.readEnd,
						rightMostOfLeft.contig.index,
						rightMostOfLeft.refStart,
						rightMostOfLeft.refEnd,
						leftMostOfRight.readStart,
						leftMostOfRight.readEnd,
						leftMostOfRight.contig.index,
						leftMostOfRight.refStart,
						leftMostOfRight.refEnd						
						);

				if (gaps < smallestGaps){
					smallestGaps = gaps;
					bestConnection = connection;
				}				
			}
			System.out.printf("                ---> best gaps = %d %d\n", smallestGaps, bestConnection.readID);



			//To get the position to return
			//Range range = contigRanges.get(left);
			//leftContigEnd = (left.getRelDir()>0) ? range.end:range.start;
			//range = contigRanges.get(right);
			//rightContigStart = (right.getRelDir()>0) ? range.start:range.end;
			//return new Pair<Integer, Integer>(leftContigEnd, rightContigStart);
		}
		//return null;
	}

	static class Range{
		int start = 0, end = 0;
		Range(int s, int e){
			start = s;
			end = e;
		}
	}


	/**
	 * Return the position on the reference that corresponds to a given position
	 * on read.
	 *  
	 * @param posInRead
	 * @param record
	 * @return
	 */
	static int positionOnRef(int posOnRead, AlignmentRecord record){
		if (posOnRead < record.readAlignmentStart() || posOnRead > record.readAlignmentEnd())
			return 0;

		if (!record.strand)
			posOnRead = record.readLength - posOnRead + 1;



		int pos = record.strand?record.readStart:(record.readLength + 1 - record.readStart);
		int posOnRef = record.refStart;

		//assert pos <= posOnRead
		for (final CigarElement e : record.alignmentCigars) {
			final int  length = e.getLength();
			switch (e.getOperator()) {
			case H :
			case S :					
			case P :
				break; // ignore pads and clips
			case I :				
				//insert
				if (pos + length < posOnRead){
					pos += length;				
				}else{
					return posOnRef;
				}
				break;
			case M ://match or mismatch				
			case EQ://match
			case X ://mismatch
				if (pos + length < posOnRead){
					pos += length;
					posOnRef += length;
				}else{
					return posOnRef + posOnRead - pos;
				}
				break;
			case D :
				posOnRef += length;
				break;
			case N :	
				posOnRef += length;
				break;								
			default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}//casse
		}//for		

		return 0;
	}



	public int fill(){
		System.out.printf("##################START########################\n"
				+ "Contig %3d (%d) -> Contig %3d (%d) Vector (%s) score = %f distance = %d\n",
				this.firstContig.index,
				this.firstContig.length(),
				this.secondContig.index,
				this.secondContig.length(),
				transVector.toString(),
				this.score,
				transVector.distance(firstContig, secondContig)				
				);

		Collections.sort(connections);
		int ret = Integer.MAX_VALUE;

		for (Connection connect:connections){
			int f = connect.fill();

			if (f <= 0)
				return 0;

			if (f < ret)
				ret = f;
		}
		return ret;
	}

	public void display(){
		System.out.printf("##################START########################\n"
				+ "Contig %3d (%d) -> Contig %3d (%d) Vector (%s) score = %f distance = %d\n",
				this.firstContig.index,
				this.firstContig.length(),
				this.secondContig.index,
				this.secondContig.length(),
				transVector.toString(),
				this.score,
				transVector.distance(firstContig, secondContig)				
				);

		Collections.sort(connections);
		for (Connection connect:connections)
			connect.display();
		System.out.println("##################END########################");
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(ContigBridge o) {
		return (int) (o.score - score);

	}	

	static class FillRecord{
		Sequence contigSequence;
		int start;
		int end;

		public String toString(){
			return contigSequence.getName() +":"+start+"-"+end;
		}
	}

	public class Connection implements Comparable<Connection>{
		ReadFilling read;
		int readID;

		int score;
		ScaffoldVector trans;
		AlignmentRecord firstAlignment, secondAlignment;

		int distanceOnRead = 0;

		Connection(ReadFilling mRead, AlignmentRecord a, AlignmentRecord b, ScaffoldVector trans){
			this.read = mRead;
			this.readID = a.readID;
			this.firstAlignment = a;
			this.secondAlignment = b;

			int aAlign = Math.abs(a.refStart - a.refEnd);
			int bAlign = Math.abs(b.refStart - b.refEnd);

			score = aAlign * bAlign / (aAlign  +bAlign);
			this.trans = trans;			

			distanceOnRead = Math.max(Math.min(b.readStart, b.readEnd) - Math.max(a.readEnd,a.readStart),
					Math.min(a.readStart, a.readEnd) - Math.max(b.readEnd,b.readStart));											
		}

		public AlignmentRecord getAlignment(Contig contig){
			if (contig == firstContig)
				return firstAlignment;
			if (contig == secondContig)
				return secondAlignment;

			return null;
		}

		void display (){
			System.out.printf("[%6d %6d] -> [%6d %6d] : [%6d %6d] -> [%6d %6d] (%s) score=%d Read %s ==> %d [%d]\n", 
					firstAlignment.refStart, firstAlignment.refEnd, secondAlignment.refStart, secondAlignment.refEnd,
					firstAlignment.readStart, firstAlignment.readEnd, secondAlignment.readStart, secondAlignment.readEnd,
					trans.toString(),					
					score, read.readSequence.getName(),
					trans.distance(firstContig, secondContig),
					distanceOnRead);
		}

		@Deprecated
		int fill (){			
			System.out.printf("FILL %s for %4d and %4d : [%4d -> %4d]  [%4d -> %4d]\n",					
					read.readSequence.getName(),
					firstContig.index,
					secondContig.index,
					this.firstAlignment.readStart,
					this.firstAlignment.readEnd,
					this.secondAlignment.readStart,
					this.secondAlignment.readEnd					
					);
			int start = 0, end = read.readSequence.length();

			if (firstAlignment.readAlignmentStart() < secondAlignment.readAlignmentStart()){
				start = firstAlignment.readAlignmentEnd();
				end = secondAlignment.readAlignmentStart();
			}else{
				start = secondAlignment.readAlignmentEnd();
				end = firstAlignment.readAlignmentStart();
			}

			return read.fill(start, end);
		}
		/**
		 * Count the number of gaps (that is the number of bases that are not
		 * aligned to a contig) between 2 main contigs 
		 * @return
		 */
		public int gapsBetween(){
			int start = 0, end = read.readSequence.length();

			if (firstAlignment.readAlignmentStart() < secondAlignment.readAlignmentStart()){
				start = firstAlignment.readAlignmentEnd();
				end = secondAlignment.readAlignmentStart();
			}else{
				start = secondAlignment.readAlignmentEnd();
				end = firstAlignment.readAlignmentStart();
			}	

			if (start >= end)
				return 0;

			BitSet bitSet = new BitSet(end);
			bitSet.set(start, end);
			for (AlignmentRecord record:read.alignments){
				bitSet.clear(record.readAlignmentStart(), record.readAlignmentEnd());
			}

			return bitSet.cardinality();
		}

		public void fillFrom(Contig fromContig, SequenceBuilder seqBuilder, JapsaAnnotation anno){
			Contig toContig = 
					(fromContig == firstContig)?secondContig:firstContig;

			AlignmentRecord fromAlignment = 
					(fromContig == firstContig)?firstAlignment:secondAlignment;

			AlignmentRecord toAlignment = 
					(fromContig == firstContig)?secondAlignment:firstAlignment;

			ReadFilling readFilling = read;
			if ( fromContig.getRelDir()>0 != fromAlignment.strand){
				//swap
				readFilling = read.reverse();
				readFilling.sortAlignment();
				for (AlignmentRecord record:readFilling.alignments){
					if (record.contig == fromContig)
						fromAlignment = record;

					if (record.contig == toContig)
						toAlignment = record;	
				}
			}
			//now readFilling is good to go
			int posReadEnd   = fromAlignment.readAlignmentEnd();
			int posReadFinal = toAlignment.readAlignmentStart();// I need as far as posReadFinal

			for (AlignmentRecord record:readFilling.alignments){
				Contig contig = record.contig;
				if (contig == fromContig)
					continue;
				if (contig == toContig)
					continue;//could break
				contig.portionUsed += (1.0 + record.refEnd - record.refStart) / contig.length();

				if (posReadEnd >= posReadFinal -1)
					continue;//I can break here, but want to get porstionUsed of other contigs


				if (record.readAlignmentEnd() < posReadEnd)
					continue;				

				//assert:  posReadEnd < readEnd				
				if (record.readAlignmentStart() > posReadEnd){
					//Really need to fill in using read information
					int newPosReadEnd = Math.min(posReadFinal - 1, record.readAlignmentStart() -1);
					if (newPosReadEnd > posReadEnd){
						JapsaFeature feature = 
								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() + newPosReadEnd - posReadEnd,
										"CONTIG",readFilling.readSequence.getName(),'+',"");

						//TODO: P=0 get the orignial read name and position
						feature.addDesc(readFilling.readSequence.getName() + "+("+(posReadEnd + 1) +"," + newPosReadEnd+")");
						anno.add(feature);
						seqBuilder.append(readFilling.readSequence.subSequence(posReadEnd, newPosReadEnd));
						posReadEnd = newPosReadEnd;
					}
					if (posReadEnd + 1 >= posReadFinal)
						continue;//Done

					//Now get information on the contig from start
					if (record.strand){
						int refLeft = record.refStart;
						int refRight = record.refEnd;

						if (posReadFinal <= record.readAlignmentEnd()){
							refRight = positionOnRef(posReadFinal, record) -1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						
						JapsaFeature feature = 
								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() + refRight - refLeft +1,
										"CONTIG",contig.getName(),'+',"");
						feature.addDesc(contig.getName() + "+("+(refLeft ) +"," + refRight+")");
						anno.add(feature);
						
						seqBuilder.append(contig.contigSequence.subSequence(refLeft - 1, refRight));
					}else{//neg strain
						int refRight = record.refStart;
						int refLeft = record.refEnd;

						if (posReadFinal <= record.readAlignmentEnd()){
							refLeft = positionOnRef(posReadFinal, record) + 1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						
						JapsaFeature feature = 
								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() - refRight + refLeft +1,
										"CONTIG",contig.getName(),'+',"");
						feature.addDesc(contig.getName() + "-("+(refRight ) +"," + refLeft+")");
						anno.add(feature);
						
						seqBuilder.append(Alphabet.DNA.complement(contig.contigSequence.subSequence(refRight - 1, refLeft)));
					}
				}//if record.readAlignmentStart() > posReadEnd
				else{//Now get information on the contig from start
					if (record.strand){
						int refLeft = positionOnRef(posReadEnd, record) + 1;						
						int refRight = record.refEnd;

						if (posReadFinal <= record.readAlignmentEnd()){
							refRight = positionOnRef(posReadFinal, record) -1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						
						JapsaFeature feature = 
								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() + refRight - refLeft +1,
										"CONTIG",contig.getName(),'+',"");
						feature.addDesc(contig.getName() + "+("+(refLeft ) +"," + refRight+")");
						anno.add(feature);
						
						seqBuilder.append(contig.contigSequence.subSequence(refLeft - 1, refRight));
					}else{//neg strain						
						int refLeft = positionOnRef(posReadEnd, record) + 1;		
						int refRight = record.refStart;

						if (posReadFinal <= record.readAlignmentEnd()){
							refLeft = positionOnRef(posReadFinal, record) + 1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						
						JapsaFeature feature = 
								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() - refRight + refLeft +1,
										"CONTIG",contig.getName(),'+',"");
						feature.addDesc(contig.getName() + "-("+(refRight ) +"," + refLeft+")");
						anno.add(feature);
						
						seqBuilder.append(Alphabet.DNA.complement(contig.contigSequence.subSequence(refRight - 1, refLeft)));
					}
				}
			}			
		}
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(Connection o) {
			// TODO Auto-generated method stub
			return o.score - score;
		}

		public Pair<Integer, ArrayList<FillRecord>> fillConnection(Contig fromContig){			
			ArrayList<FillRecord> fillRecords = new  ArrayList<FillRecord> ();
			ReadFilling readFilling = this.read;
			readFilling.sortAlignment();

			System.out.printf("  Fill using read %s to connect %d (%b) to %d(%b)\n", 
					readFilling.readSequence.getName(),
					firstAlignment.contig.index, 
					firstAlignment.strand,
					secondAlignment.contig.index,
					secondAlignment.strand);

			Contig toContig = null;
			int fromStart = 0, toStart = 0;
			{
				AlignmentRecord fromAlignment, toAlignment;

				if (fromContig == firstContig){
					fromAlignment = firstAlignment;
					toAlignment = secondAlignment;						
					toContig = secondContig;
				}else{
					//from is the second
					fromAlignment = secondAlignment;
					toAlignment   = firstAlignment;				
					toContig = firstContig;
				}
				fromStart = fromAlignment.readAlignmentStart();
				toStart = toAlignment.readAlignmentStart();
			}

			if (toStart < fromStart){
				System.out.printf("  -----> Reverse read\n");
				readFilling = readFilling.reverse();
				readFilling.sortAlignment();
			}

			int startIndex = -1, endIndex = -1; 
			for (int i = 0; i < readFilling.alignments.size();i++){
				if (readFilling.alignments.get(i).contig == fromContig)
					startIndex = i;

				if (readFilling.alignments.get(i).contig == toContig)
					endIndex = i;
			}

			if (startIndex < 0 || endIndex < 0)
				throw new RuntimeException("Error in fillFrom()");

			if (startIndex >= endIndex){
				System.out.println("The two overlap, good news, but need to do some thing " + startIndex + " vs " + endIndex);
				//continue;
			}
			//assert startIndex < endIndex

			FillRecord fillRecord = new FillRecord();
			fillRecord.contigSequence = readFilling.alignments.get(startIndex).contig.contigSequence;

			fillRecord.start = readFilling.alignments.get(startIndex).refStart;
			fillRecord.end = readFilling.alignments.get(startIndex).refEnd;

			fillRecords.add(fillRecord);//Record from the from contig

			int posOnRead = readFilling.alignments.get(startIndex).readAlignmentEnd();

			for (int i = 0; i < readFilling.alignments.size();i++){
				AlignmentRecord alignment = readFilling.alignments.get(i);
				if (alignment.readAlignmentStart() > posOnRead)
					break;

				if (alignment.readAlignmentEnd() <= posOnRead)//want to add something that positively extend the sequence
					continue;

				//assert alignmentStart =<posOnRead<alignmentEnd
				int outRefPos = positionOnRef(posOnRead, alignment);

				//outRefPos == inRefPos
				if (outRefPos <=0)
					throw new RuntimeException("Error code = 100");


				fillRecord = new FillRecord();
				fillRecord.contigSequence = alignment.contig.contigSequence;
				if (alignment.strand){
					//contig same direction
					fillRecord.start = outRefPos + 1;
					fillRecord.end = alignment.refEnd;					
				}else{
					fillRecord.start = outRefPos - 1;
					fillRecord.end = alignment.refStart;
				}					
				posOnRead = alignment.readAlignmentEnd(); 
				fillRecords.add(fillRecord);
				if (alignment.contig == toContig){
					return new Pair<Integer, ArrayList<FillRecord>> (0, fillRecords);
				}
			}
			//if you get here meaning not connected yet

			//reverse if not reverted before
			if (readFilling != this.read){
				readFilling = read;
			}else{
				readFilling = read.reverse();
				readFilling.sortAlignment();
			}

			ArrayList<FillRecord> rightFillRecords = new  ArrayList<FillRecord> ();

			startIndex = -1;
			endIndex = -1; 
			for (int i = 0; i < readFilling.alignments.size();i++){
				if (readFilling.alignments.get(i).contig == toContig)
					startIndex = i;

				if (readFilling.alignments.get(i).contig == fromContig)
					endIndex = i;
			}

			if (startIndex < 0 || endIndex < 0)
				throw new RuntimeException("Error in fillFrom()");

			if (startIndex >= endIndex){
				System.out.println("The two overlap, good news, but need to do some thing " + startIndex + " vs " + endIndex);
				//continue;
			}


			fillRecord = new FillRecord();
			fillRecord.contigSequence = readFilling.alignments.get(startIndex).contig.contigSequence;

			fillRecord.start = readFilling.alignments.get(startIndex).refEnd;
			fillRecord.end = readFilling.alignments.get(startIndex).refStart;
			rightFillRecords.add(fillRecord);//Record from the from contig

			int posOnReadFromRight = readFilling.alignments.get(startIndex).readAlignmentEnd();

			for (int i = 0; i < readFilling.alignments.size();i++){
				AlignmentRecord alignment = readFilling.alignments.get(i);
				if (alignment.readAlignmentStart() > posOnReadFromRight)
					break;

				if (alignment.readAlignmentEnd() <= posOnReadFromRight)//want to add something that positively extend the sequence
					continue;

				//assert alignmentStart =<posOnRead<alignmentEnd
				int outRefPos = positionOnRef(posOnReadFromRight, alignment);

				//outRefPos == inRefPos
				if (outRefPos <=0)
					throw new RuntimeException("Error code = 100");


				fillRecord = new FillRecord();
				fillRecord.contigSequence = alignment.contig.contigSequence;
				if (alignment.strand){
					//contig same direction
					fillRecord.end = outRefPos + 1;
					fillRecord.start = alignment.refEnd;					
				}else{
					fillRecord.end = outRefPos -1;
					fillRecord.start = alignment.refStart;
				}					
				posOnReadFromRight = alignment.readAlignmentEnd();
				rightFillRecords.add(fillRecord);				
			}
			posOnRead = read.readSequence.length() - posOnRead + 1;

			for (int i = rightFillRecords.size() - 1; i >=0;i--){
				fillRecords.add(rightFillRecords.get(i));			
			}

			return new Pair<Integer, ArrayList<FillRecord>> (Math.abs(posOnRead - posOnReadFromRight), fillRecords);
		}
	}
}
