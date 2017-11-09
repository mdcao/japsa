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

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;

/**
 * Class represent Illumina sequencing
 * 
 * @author minhduc
 *TODO: make full class
 */
public class IlluminaSequencing{

	/**
	 * Simulate MiSeq
	 * @param fragment
	 * @param o1
	 * @param o2
	 * @param rnd
	 * @throws IOException
	 */
	public static void simulatePaired(Sequence fragment, SequenceOutputStream o1 , SequenceOutputStream o2, Random rnd) throws IOException{		
		//double r = 
		int len = Math.min(default_illen, fragment.length());
		simulatePaired(fragment, len, o1 , o2, rnd);
	}
	
	private final static int default_illen = 250;
	private final static int maximum_illen = 300;
	
	public static int IlluminaReadDefaultLength() {
		return default_illen;
	}
	
	public static int IlluminaReadMaximumLength() {
		return maximum_illen;
	}
	
	// can specify read length
	public static void simulatePaired(Sequence fragment, int illen, SequenceOutputStream o1 , SequenceOutputStream o2, Random rnd) throws IOException{		

		double snp = 0.01;
		double del = 0.0001;				
		double ins = 0.0001;
		double ext = 0.2;

		//double r = 
		int len = Math.min(maximum_illen, Math.min(illen, fragment.length()));
		
		String name = fragment.getName();				

		SequenceBuilder read1 = SequencingSimulation.simulateRead(fragment, len, snp, del, ins, ext, rnd);
		SequenceBuilder read2 = SequencingSimulation.simulateRead(Alphabet.DNA.complement(fragment), len, snp, del, ins, ext, rnd);

		if (rnd.nextBoolean()){
			o1.print("@");
			o1.print(name);						
			o1.print("\n");
			for (int i = 0; i < read1.length();i++)
				o1.print(read1.charAt(i));
			o1.print("\n+\n");

			for (int i = 0; i < read1.length();i++)
				o1.print("I");
			o1.print("\n");

			o2.print("@");
			o2.print(name);						
			o2.print("\n");
			for (int i = 0; i < read2.length();i++)
				o2.print(read2.charAt(i));
			o2.print("\n+\n");

			for (int i = 0; i < read2.length();i++)
				o2.print("I");
			o2.print("\n");			
		}else{
			o2.print("@");
			o2.print(name);						
			o2.print("\n");
			for (int i = 0; i < read1.length();i++)
				o2.print(read1.charAt(i));
			o2.print("\n+\n");

			for (int i = 0; i < read1.length();i++)
				o2.print("I");
			o2.print("\n");

			o1.print("@");
			o1.print(name);						
			o1.print("\n");
			for (int i = 0; i < read2.length();i++)
				o1.print(read2.charAt(i));
			o1.print("\n+\n");

			for (int i = 0; i < read2.length();i++)
				o1.print("I");
			o1.print("\n");
		}
	}
	
}