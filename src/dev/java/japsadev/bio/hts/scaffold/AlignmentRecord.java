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

import java.util.ArrayList;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class AlignmentRecord implements Comparable<AlignmentRecord> {
	static final double matchCost = 0;
	int score;

	public String readID;
	Contig contig;

	public int refStart, refEnd;  //1-based position on ref of the start and end of the alignment
	
	//Position on read of the start and end of the alignment (using the direction of read) 
	//readStart map to refStart, readEnd map to refEnd. 
	//readStart < readEnd if strand = true, else readStart > readEnd
	public int readStart = 0, readEnd = 0;	
	
	//read length
	public int readLength = 0;		

	public boolean strand = true;//positive
	public boolean useful = false;
	//SAMRecord mySam;
	
	ArrayList<CigarElement> alignmentCigars = new ArrayList<CigarElement>();
	

	//public int readLeft, readRight, readAlign, refLeft, refRight, refAlign;
	//left and right are in the direction of the reference sequence
	
	public AlignmentRecord(String readID, int refStart, int refEnd, int readLength, 
			int readStart, int readEnd, boolean strand, boolean useful, Contig contig, int score){
		this.readID = readID;
		this.contig = contig;
		this.refStart = refStart;
		this.refEnd = refEnd;
		
		this.readLength = readLength;
		this.readStart = readStart;//1-index
		this.readEnd = readEnd;//1-index
		this.strand = strand;
		this.useful = useful;			
		this.contig = contig;
		this.score = score;
	}
	public AlignmentRecord(SAMRecord sam, Contig ctg) {
//		readID = Integer.parseInt(sam.getReadName().split("_")[0]);
		readID = sam.getReadName();

		contig = ctg;

		//mySam = sam;
		refStart = sam.getAlignmentStart();
		refEnd = sam.getAlignmentEnd();
		
		Cigar cigar = sam.getCigar();			
		boolean enterAlignment = false;						
		//////////////////////////////////////////////////////////////////////////////////

		for (final CigarElement e : cigar.getCigarElements()) {
			alignmentCigars.add(e);
			final int  length = e.getLength();
			switch (e.getOperator()) {
			case H :
			case S :					
			case P : //pad is a kind of clipped
				if (enterAlignment)
					readEnd = readLength;
				readLength += length;
				break; // soft clip read bases
			case I :	                	
			case M :					
			case EQ :
			case X :
				if (!enterAlignment){
					readStart = readLength + 1;
					enterAlignment = true;
				}
				readLength += length;
				break;
			case D :
			case N :
				if (!enterAlignment){
					readStart = readLength + 1;
					enterAlignment = true;
				}
				break;				
			default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}//case
		}//for
		if (readEnd == 0)
			readEnd = readLength;
		//these temporary variable to determine usefulness
		int readLeft = readStart -1;
		int readRight = readLength - readEnd;

		int refLeft = refStart - 1;
		int refRight = contig.length() - refEnd;
		score = refEnd + 1 - refStart;
		if (sam.getReadNegativeStrandFlag()){			
			strand = false;
			//need to convert the alignment position on read the correct direction 
			readStart = 1 + readLength - readStart;
			readEnd = 1 + readLength - readEnd;
		}

		if (
				(readLeft < ScaffoldGraph.marginThres || refLeft < ScaffoldGraph.marginThres) &&
				(readRight  < ScaffoldGraph.marginThres || refRight < ScaffoldGraph.marginThres) &&
				score > ScaffoldGraph.minContigLength
			)
			useful = true;

	}
	
	
	public int readAlignmentStart(){
		return Math.min(readStart,readEnd);
	
	}
	
	public int readAlignmentEnd(){
		return Math.max(readStart,readEnd);
	}

	public String toString() {
		return contig.index    
				+ " " + refStart 
				+ " " + refEnd
				+ " " + contig.length()
				+ " " + readStart 
				+ " " + readEnd
				+ " " + readLength				
				+ " " + strand;
	}
	
	public String pos() {			
		return  
				refStart 
				+ " " + refEnd
				+ " " + contig.length()
				+ " " + readStart
				+ " " + readEnd
				+ " " + readLength
				+ " " + score
				+ " " + strand
				;
	}
	// return same alignment but with reversed read
	//TODO: change to object self-editing function?
	public AlignmentRecord reverseRead(){
		AlignmentRecord revAlign = new AlignmentRecord(readID, refStart, refEnd, readLength, 
		readLength - readStart + 1, readLength - readEnd + 1, !strand, useful, contig, score);
		
		revAlign.alignmentCigars = alignmentCigars;

		return revAlign;
	}
	public AlignmentRecord clones(){
		AlignmentRecord align = new AlignmentRecord(readID, refStart, refEnd, readLength,
				readStart, readEnd, strand, useful, contig, score);
		
		align.alignmentCigars = alignmentCigars;

		return align;
	}
	public void copy(AlignmentRecord rec){
		readID = rec.readID;
		contig = rec.contig;
		refStart = rec.refStart;
		refEnd = rec.refEnd;
		
		readLength = rec.readLength;
		readStart = rec.readStart;//1-index
		readEnd = rec.readEnd;//1-index
		strand = rec.strand;
		useful = rec.useful;			
		alignmentCigars = rec.alignmentCigars;
		contig = rec.contig;
		score = rec.score;
	}
	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(AlignmentRecord o) {			
		return readAlignmentStart() - o.readAlignmentStart();
	}
}