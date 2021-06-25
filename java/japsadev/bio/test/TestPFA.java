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
 * 09/12/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.bio.test;

import java.util.Random;

import japsa.bio.alignment.ProbFSM;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;

/**
 * @author minhduc
 *
 */
public class TestPFA {
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int length = 1000;
		Alphabet dna = Alphabet.DNA4();
		Random rand = new Random(1);
		Sequence mSeq, gSeq;
		ProbFSM fa;
		ProbFSM.Emission emission;
		
		/********************************************************************/
		System.out.println("==================================================");
		mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rand);
		fa = new ProbFSM.ProbThreeSM(mSeq);
		
		gSeq = fa.generate(rand);
		
		emission = fa.align(gSeq);
		
		System.out.println(emission.myCost);		
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();		
		
		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		fa.showProb();
		/********************************************************************/
		System.out.println("==================================================");
		mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rand);
		fa = new ProbFSM.ProbThreeSM(mSeq);
		
		gSeq = fa.generate(rand);
		
		emission = fa.align(gSeq);
		
		System.out.println(emission.myCost);		
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();		
		
		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		fa.showProb();
		/********************************************************************/
		System.out.println("==================================================");
		mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rand);
		fa = new ProbFSM.ProbThreeSM(mSeq);
		
		gSeq = fa.generate(rand);
		
		emission = fa.align(gSeq);
		
		System.out.println(emission.myCost);		
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();		
		
		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		fa.showProb();
		/********************************************************************/
		System.out.println("==================================================");
		mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rand);
		fa = new ProbFSM.ProbThreeSM(mSeq);
		
		gSeq = fa.generate(rand);
		
		emission = fa.align(gSeq);
		
		System.out.println(emission.myCost);		
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();		
		
		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		fa.showProb();
		
		/********************************************************************/
		System.out.println("==================================================");
		mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rand);
		fa = new ProbFSM.ProbThreeSM(mSeq);
		
		gSeq = fa.generate(rand);
		
		emission = fa.align(gSeq);
		
		System.out.println(emission.myCost);		
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();		
		
		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		fa.showProb();
		/********************************************************************/
		System.out.println("==================================================");
		mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rand);
		fa = new ProbFSM.ProbThreeSM(mSeq);
		
		gSeq = fa.generate(rand);
		
		emission = fa.align(gSeq);
		
		System.out.println(emission.myCost);		
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();		
		
		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		fa.showProb();
		/********************************************************************/
		System.out.println("==================================================");
		mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rand);
		fa = new ProbFSM.ProbThreeSM(mSeq);
		
		gSeq = fa.generate(rand);
		
		emission = fa.align(gSeq);
		
		System.out.println(emission.myCost);		
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();		
		
		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		fa.showProb();
		/********************************************************************/
		System.out.println("==================================================");
		mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rand);
		fa = new ProbFSM.ProbThreeSM(mSeq);
		
		gSeq = fa.generate(rand);
		
		emission = fa.align(gSeq);
		
		System.out.println(emission.myCost);		
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();		
		
		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		fa.showProb();
		
		/********************************************************************/
		System.out.println("==================================================");
		mSeq = Sequence.random(dna, length, new double[]{.25,.25,.25,.25}, rand);
		fa = new ProbFSM.ProbThreeSM(mSeq);
		
		gSeq = fa.generate(rand);
		
		emission = fa.align(gSeq);
		
		System.out.println(emission.myCost);		
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();		
		
		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		

		emission = fa.align(gSeq);		
		System.out.println(emission.myCost);
		fa.updateCount(emission);
		fa.reEstimate();
		fa.resetCount();
		
		fa.showProb();
		
	}

}
