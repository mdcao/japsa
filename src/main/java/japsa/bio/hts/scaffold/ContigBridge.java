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


import java.util.ArrayList;
import java.util.Collections;




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
	final int order;

	private double score = 0;//more is better
	private ScaffoldVector transVector = null;
	private ArrayList<Connection> connections;//a list of connections that make up this


	public ContigBridge(Contig c1, Contig c2, int ind){
		firstContig = c1;
		secondContig = c2;
		order = ind;		
		hashKey = makeHash(c1.index,c2.index, order);

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
			AlignmentRecord a, 
			AlignmentRecord b, 
			ScaffoldVector trans, 
			double sc){

		if (transVector == null){
			transVector = trans;
			score = sc;			
			connections.add(new Connection(readSequence, a,b,trans));			
		}else{
			boolean alreadyIn = false;
			for (int i = 0; i < connections.size();i++){
				Connection oldConnect = connections.get(i); 
				if (connections.get(i).readID / 3 == a.readID /3){
					if (sc > oldConnect.score){
						//replace only if better						
						Connection newConnect = new Connection(readSequence, a,b,trans);
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
			Connection newConnect = new Connection(readSequence, a,b,trans);
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

	public class Connection implements Comparable<Connection>{
		ReadFilling read;
		int readID;
		//int aReadStart, aReadEnd, bReadStart, bReadEnd;
		//int aRefStart, aRefEnd, bRefStart, bRefEnd;
		int score;
		ScaffoldVector trans;
		AlignmentRecord a,b;

		int distanceOnRead = 0;

		Connection(ReadFilling mRead, AlignmentRecord a, AlignmentRecord b, ScaffoldVector trans){
			this.read = mRead;
			this.readID = a.readID;
			this.a = a;
			this.b = b;

			//aReadStart = a.readStart;
			//aReadEnd = a.readEnd;
			//aRefStart = a.refStart;
			//aRefEnd = a.refEnd;

			//bReadStart = b.readStart;
			//bReadEnd = b.readEnd;
			//bRefStart = b.refStart;
			//bRefEnd = b.refEnd;

			int aAlign = Math.abs(a.refStart - a.refEnd);
			int bAlign = Math.abs(b.refStart - b.refEnd);

			score = aAlign * bAlign / (aAlign  +bAlign);
			this.trans = trans;			

			distanceOnRead = Math.max(Math.min(b.readStart, b.readEnd) - Math.max(a.readEnd,a.readStart),
					Math.min(a.readStart, a.readEnd) - Math.max(b.readEnd,b.readStart));											
		}

		void display (){
			System.out.printf("[%6d %6d] -> [%6d %6d] : [%6d %6d] -> [%6d %6d] (%s) score=%d Read %s ==> %d [%d]\n", 
					a.refStart, a.refEnd, b.refStart, b.refEnd,
					a.readStart, a.readEnd, b.readStart, b.readEnd,
					trans.toString(),					
					score, read.readSequence.getName(),
					trans.distance(firstContig, secondContig),
					distanceOnRead);
		}

		int fill (){			
			System.out.printf("FILL %s for %4d and %4d : [%4d -> %4d]  [%4d -> %4d]\n",					
					read.readSequence.getName(),
					firstContig.index,
					secondContig.index,
					this.a.readStart,
					this.a.readEnd,
					this.b.readStart,
					this.b.readEnd					
					);
			int start = 0, end = read.readSequence.length();

			if (a.readAlignmentStart() < b.readAlignmentStart()){
				start = a.readAlignmentEnd();
				end = b.readAlignmentStart();
			}else{
				start = b.readAlignmentEnd();
				end = a.readAlignmentStart();
			}

			return read.fill(start, end);						
		}

		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(Connection o) {
			// TODO Auto-generated method stub
			return o.score - score;
		}
	}
}
