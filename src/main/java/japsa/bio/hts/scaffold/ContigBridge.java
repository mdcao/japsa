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

import japsa.seq.Sequence;


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
	}
	
	public static String makeHash(int aIndex, int bIndex, int order){		
		return aIndex+"#"+bIndex + "#" + order;
	}
	
	public boolean consistentWith(ScaffoldVector aVector){
		return (aVector.dir == transVector.dir)
				&& (aVector.pos * 1.0 / aVector.pos > 0.9)
				&& (aVector.pos * 1.0 / aVector.pos < 1.1)
				;
	}
	
	public double addConnection(){
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

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(ContigBridge o) {
		return (int) (o.score - score);
	}
	
	
	public class Connection{
		Sequence read;	
		//need a read, two contigs, and positions of mapping
		ScaffoldVector vector(){			
			return null;
		}
		
	}
}
