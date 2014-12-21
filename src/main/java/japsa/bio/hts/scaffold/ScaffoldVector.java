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

/**
 * Implementation of a vector of relative position of a contig in its scaffold 
 * @author minhduc
 *
 */
public class ScaffoldVector{
	int pos = 0;
	int dir = 1;
	
	public ScaffoldVector(){

	}
	public ScaffoldVector(int p, int d){
		pos = p;
		dir = d;
	}
		
	public static ScaffoldVector reverse(ScaffoldVector v){
		ScaffoldVector rev = new ScaffoldVector();
		if (v.dir > 0){
			rev.dir = 1;
			rev.pos = - v.pos;
		}else{
			rev.dir = -1;
			rev.pos = v.pos;
		}		

		return rev;		
	}
	/**
	 * Compose two vectors: a to b is v1, b to c is v2. a to c is v1 * v2
	 * @param v1
	 * @param v2
	 * @return
	 */
	public static ScaffoldVector composition(ScaffoldVector v1, ScaffoldVector v2){
		ScaffoldVector ret = new ScaffoldVector();

		ret.pos = v2.pos + v2.dir * v1.pos; 
		ret.dir = v1.dir * v2.dir;	

		return ret;		
	}
	/**
	 * @return the pos
	 */
	public int getPos() {
		return pos;
	}
	/**
	 * @return the dir
	 */
	public int getDir() {
		return dir;
	}
	

}