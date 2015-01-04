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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class AlignmentRecord implements Comparable<AlignmentRecord> {
	static final double matchCost = 0;

	public int readID;
	Contig contig;

	public int refStart, refEnd;  //position on ref of the start and end of the alignment

	//Position on read of the start and end of the alignment 
	public int readStart = 0, readEnd = 0;
	//readStart map to refStart, readEnd map to refEnd. 
	//readStart < readEnd if strand = true, else readStart > readEnd

	//read length
	public int readLength = 0;		

	public boolean strand = true;//positive
	public boolean useful = false;


	public int readLeft, readRight, readAlign, refLeft, refRight, refAlign;
	//left and right are in the direction of the reference sequence
	
	public AlignmentRecord(SAMRecord sam, Contig ctg) {
		readID = Integer.parseInt(sam.getReadName().split("_")[0]);
		contig = ctg;

		refStart = sam.getAlignmentStart();
		refEnd = sam.getAlignmentEnd();		
		
		Cigar cigar = sam.getCigar();			
		boolean enterAlignment = false;						
		//////////////////////////////////////////////////////////////////////////////////

		for (final CigarElement e : cigar.getCigarElements()) {
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
			}//casse
		}//for
		if (readEnd == 0)
			readEnd = readLength;

		readLeft = readStart;
		readRight = readLength - readEnd;
		readAlign = readEnd + 1 - readStart;

		refLeft = refStart;
		refRight = contig.length() - refEnd;
		refAlign = refEnd + 1 - refStart;

		if (sam.getReadNegativeStrandFlag()){			
			strand = false;
			//need to convert the aligment position on read the the right direction
			readStart = readLength - readStart;
			readEnd = readLength - readEnd;
		}

		int gaps = 500;
		int extend = 250;
		//only useful if
		if ((readLeft > refLeft + gaps || readRight > gaps + refRight)
				&& (readLeft < extend || refLeft < extend)
				&& (readRight  < extend || refRight < extend)
				)
			useful = true;
	}
	
	
	public int readAlignmentStart(){
		return strand?readStart:readEnd;
	}
	
	public int readAlignmentEnd(){
		return strand?readEnd:readStart;
	}

	public String toString() {
		return contig.index    
				+ " " + refStart 
				+ " " + refEnd
				+ " " + contig.length()
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
				+ " " + refLeft
				+ " " + refAlign
				+ " " + refRight
				+ " " + readLeft
				+ " " + readAlign
				+ " " + readRight
				+ " " + strand;
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(AlignmentRecord o) {			
		return readAlignmentStart() - o.readAlignmentStart();
	}
}