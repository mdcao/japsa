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
 * 31/12/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.bio.hts.scaffold;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import java.util.ArrayList;
import java.util.Collections;


public class ReadFilling{
	/**
	 * The read sequence
	 */
	Sequence readSequence;		
	private boolean sorted = false;

	ArrayList<AlignmentRecord> alignments;	

	public ReadFilling(Sequence read, ArrayList<AlignmentRecord> alignmentList){
		readSequence = read;
		alignments = alignmentList;
	}

	public Sequence getReadSequence(){
		return readSequence;
	}

	public ArrayList<AlignmentRecord> getAlignmentRecords(){
		return alignments;
	}
	
	public void print(){
		System.out.printf(" Read filling: %s with %ssorted alignments: \n",readSequence.getName(),sorted?"":"un");
		for(AlignmentRecord rec: alignments)
			System.out.println("\t" + rec.readID + ": " + rec);
	}
	
	public ReadFilling reverse(){
		//return an (conceptually the same) read filling with the a reverse read
		Sequence revRead = Alphabet.DNA.complement(readSequence);
		revRead.setName("REV"+readSequence.getName());
		ArrayList<AlignmentRecord> revAlignments = new ArrayList<AlignmentRecord>(); 
		
		for (AlignmentRecord alignment:alignments)
			revAlignments.add(alignment.reverseRead());
		
		ReadFilling revFilling = new ReadFilling(revRead, revAlignments);		
		return revFilling;
	}
	
	public void sortAlignment(){
		if (!sorted){				
			Collections.sort(alignments);
			sorted = true;
		}
	}
}