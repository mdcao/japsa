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

import java.io.IOException;
import java.util.Random;

import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;

/**
 * Implement PacBio sequencing
 * @author minhduc
 *
 */
public class PacBioSequencing{
	static int PACBIO_ADAPTER_LENGTH = 46;


	public static void simulatePacBio(Sequence fragment, int pblen, SequenceOutputStream o, Random rnd) throws IOException{
		double snp = 0.01;
		double ins = 0.1;	
		double del = 0.04;
		double ext = 0.4;

		//int len = (int) (fragment.length() * .9);
		String name = fragment.getName();		
		SequenceBuilder read	=  SequencingSimulation.simulateRead(fragment, fragment.length(), snp, del, ins, ext, rnd);		

		o.print("@");
		o.print(name);						
		o.print("\n");
		for (int i = 0; i < read.length();i++)
			o.print(read.charAt(i));
		o.print("\n+\n");

		for (int i = 0; i < read.length();i++)
			o.print("E");
		o.print("\n");

	}
}