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

/*                           Revision History                                
 * 28/08/2016 - Minh Duc Cao: Created                                        
 ****************************************************************************/
package japsa.bio.sim;

import java.util.Random;

import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;

/**
 * @author minhduc
 *
 */
public class SequencingSimulation {
	/**
	 * Simulate a read from the start of a fragment
	 * @param fragment
	 * @param len
	 * @param snp
	 * @param indel
	 * @param ext
	 * @param rnd
	 * @return
	 */
	public static SequenceBuilder simulateRead(Sequence fragment, int len,  double snp, double del, double ins, double ext, Random rnd){		
		//accumulative prob
		double aSNP = snp;
		double aDel = aSNP + del;
		double aIns = aDel + ins;

		//why - 10?/
		len = Math.min(len, fragment.length()); 

		//Sequence read = new Sequence(fragment.alphabet(), len);
		SequenceBuilder sb = new SequenceBuilder(fragment.alphabet(), len);

		//int mIndex = 0; 
		int fIndex = 0;
		//mIndex < len &&
		for (; sb.length() < len && fIndex < fragment.length();){
			byte base = fragment.getBase(fIndex);
			double r = rnd.nextDouble(); 
			if (r < aSNP){
				//simulate a SNP aka mismatch
				sb.append((byte) ((base + rnd.nextInt(3)) % 4));
				//read.setBase(mIndex, (byte) ((base + rnd.nextInt(3)) % 4));
				//mIndex ++;
				fIndex ++;				
			}else if (r < aDel){
				do{
					fIndex ++;
				}while (rnd.nextDouble() < ext);
			}else if (r < aIns){
				//insertion
				do{
					sb.append((byte) (rnd.nextInt(4)));
					//mIndex ++;
				}while (rnd.nextDouble() < ext);
			}else{
				sb.append(base);
				//read.setBase(mIndex, base);
				//mIndex ++;
				fIndex ++;
			}//else
		}//for
		//for (;mIndex < len;mIndex ++){
		//	//pad in random to fill in
		//	read.setBase(mIndex, (byte) rnd.nextInt(4));
		//}
		return sb;
	}

}
