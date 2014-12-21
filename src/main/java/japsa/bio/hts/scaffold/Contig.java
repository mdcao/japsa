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

import japsa.seq.Sequence;

public class Contig implements Comparable<Contig>{
	int index;
	ScaffoldVector myVector;//relative position to the head contig of my scaffold
	Sequence contigSequence;//the sequence of the contig	
	double coverage = 1.0;
	
	int start;//relative to the head contig of the scaffold
	int end;//relative to the head contig of the scaffold
	
	public Contig(int index, Sequence seq){
		this.index = index;
		contigSequence = seq;
		myVector = new ScaffoldVector(0,1);
	}		
	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */		
	@Override
	public int compareTo(Contig o) {
		return start - o.start;
	}
	
	public String getName(){
		return contigSequence.getName();
	}
	
	public int getIndex(){
		return index;
	}
	
	public void composite(ScaffoldVector aVector){
		myVector = ScaffoldVector.composition(myVector, aVector);		
	}	
	/**
	 * Relative position to the head of the scaffold
	 * @return
	 */
	public int getRelPos(){
		return myVector.pos;
	}
	
	/**
	 * @return the start
	 */
	public int getStartPosition() {
		return start;
	}
	/**
	 * @return the end
	 */
	public int getEndPosition() {
		return end;
	}
	
	
	public int getRelDir(){
		return myVector.dir;
	}
	
	
	public ScaffoldVector getVector(){
		return myVector;
	}
	
	public int length(){
		return contigSequence.length();
	}		
	
	public double getCoverage(){
		return coverage;
	}
	public void setCoverage(double cov){
		coverage =  cov;
	}
	
	public void computePosition(){
		if (myVector.dir > 0){
			start = getRelPos();
			end = getRelPos() + length();
		}else{
			end   = getRelPos();
			start = getRelPos() - length();
		}
		
	}
	
}