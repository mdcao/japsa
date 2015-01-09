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

package japsa.bio.hts.scaffold;

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

	
	public ReadFilling reverse(){
		//return an (concentually the same) read filling with the a reverse read
		Sequence revRead = Alphabet.DNA.complement(readSequence);
		revRead.setName("REV"+readSequence.getName());
		ArrayList<AlignmentRecord> revAlignments = new ArrayList<AlignmentRecord>(); 
		
		for (AlignmentRecord alignment:alignments){
			AlignmentRecord revAlign = new AlignmentRecord();
			revAlign.readID = alignment.readID;
			revAlign.contig = alignment.contig;
			revAlign.refStart = alignment.refStart;
			revAlign.refEnd = alignment.refEnd;
			
			revAlign.readLength = alignment.readLength;
			revAlign.readStart = alignment.readLength - alignment.readStart + 1;//1-index
			revAlign.readEnd = alignment.readLength - alignment.readEnd + 1;//1-index
			revAlign.strand = !alignment.strand;
			revAlign.useful = alignment.useful;			
			revAlign.alignmentCigars = alignment.alignmentCigars;
			revAlign.contig = alignment.contig;
			revAlign.score = alignment.score;
			//System.out.println("      switch " + alignment.readStart + " and " + alignment.readEnd
			//		+ " to " + revAlign.readStart + " and " + revAlign.readEnd);
			revAlignments.add(revAlign);
		}
		
		ReadFilling revFilling = new ReadFilling(revRead, revAlignments);		
		return revFilling;
	}
	
	public void sortAlignment(){
		if (!sorted){				
			Collections.sort(alignments);
			sorted = true;
		}
	}
	
	/**
	 * Return left over (not filled) bases
	 * @param start
	 * @param end
	 * @return
	 */
	public int fill(int start, int end){
		System.out.println("Read " + readSequence.getName() + " length " + readSequence.length());
		if (start >= end){
			System.out.println("========FILLED==========");
			return 0;
		}

		sortAlignment();		

		int gaps = 0;
		int     myStart = start; 
				//myEnd = end;

		for (AlignmentRecord record:alignments){
			if (record.readAlignmentEnd() < start)
				continue;

			if (record.readAlignmentStart() > end)
				break;			

			System.out.printf("============== %5d -> %5d Contig %d: %5d %5d\n",
					record.readAlignmentStart(),
					record.readAlignmentEnd(),
					record.contig.index,
					record.strand?record.refStart:record.refEnd,
							record.strand?record.refEnd:record.refStart
					);

			if (record.readAlignmentStart() > myStart)
				gaps += record.readAlignmentStart() - myStart;

			myStart = Math.max(myStart, record.readAlignmentEnd());
		}		
		if (myStart < end)
			gaps += end - myStart;

		System.out.printf("Left %d bases\n",gaps);
		return gaps;
	}
}