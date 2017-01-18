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


/**
 * Implementation of a vector of relative position of a contig in its scaffold 
 * @author minhduc
 *
 */
public class ScaffoldVector{
	int magnitude = 0; //distance of two contigs' starting points *---> <---* (with sign follows the first contig's direction)
	int direction = 1; //relative direction of those two (+/- for same/opposite directions)
	
	public ScaffoldVector(){
		//magnitude = 0
		//direction = 1;
	}
	public ScaffoldVector(int p, int d){
		magnitude = p;
		direction = d;
	}
	// reverse vector a -> b to b -> a	(a and b are contigs)
	public static ScaffoldVector reverse(ScaffoldVector v){
		ScaffoldVector rev = new ScaffoldVector();
		if (v.direction > 0){
			rev.direction = 1;
			rev.magnitude = - v.magnitude;
		}else{
			rev.direction = -1;
			rev.magnitude = v.magnitude;
		}
		return rev;		
	}
	
	/**
	 * Return the distance between two closest tips of two contigs
	 * that relative position of fContig with regard to tContig (tContig->fContig) is represented by *this*
	 * A negative value indicate the overlap of two contigs.
	 * @param fContig
	 * @param tContig
	 * @return
	 */
	
	public int distance(Contig tContig, Contig fContig){
		
		int tS = 0, tE = tContig.length(),
			fS, fE;
		
		if (direction > 0){
			fS = magnitude;
			fE = magnitude + fContig.length();
		}else{
			fE = magnitude;
			fS = magnitude - fContig.length();
		}		
		//System.out.printf("tS=%d tE=%d fS=%d fE=%d fS-tE=%d tS-fE = %d ret=%d\n",tS, tE, fE, fE, fS-tE,tS-fE,Math.max(fS - tE, tS - fE));
		
		//FIXME: not handle the case that contig A contain contigB and via verse		
		return Math.max(fS - tE, tS - fE);
	}

	/**
	 * Compose two vectors: a -> b is v2, b -> c is v1. returned a -> c is v1 * v2
	 * Warning: the parameters' order doesn't follow normal intuition. USE WITH CARE!!!
	 * @param v1
	 * @param v2
	 * @return
	 */
	public static ScaffoldVector composition(ScaffoldVector v1, ScaffoldVector v2){
		ScaffoldVector ret = new ScaffoldVector();

		ret.magnitude = v2.magnitude + v2.direction * v1.magnitude; 
		ret.direction = v1.direction * v2.direction;	

		return ret;		
	}
	
	public String toString(){
		return "<" + magnitude + ", " + direction + ">";
	}
	/**
	 * @return the magnitude
	 */
	public int getMagnitute() {
		return magnitude;
	}
	/**
	 * @return the direction
	 */
	public int getDirection() {
		return direction;
	}
	
	/**
	 * Set the new magnitude
	 * @param magnitude
	 */
	public void setMagnitute(int magnitude){
		this.magnitude=magnitude;
		
	}
}