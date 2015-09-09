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
 * 19 Mar 2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.util;

import java.util.ArrayList;
import java.util.Arrays;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;

/**
 * A collection of utilities to analyse HTS, based on HTS library
 * @author minhduc
 *
 */
public class HTSUtilities {

	
	/**
	 * Get the read subsequence that spans the gene. The method look at an alignment,
	 * estimated the position on reads that might have mapped to the start and the end
	 * of the gene, and extracts the subsequence mapped to the whole gene plus the flank 
	 * 
	 * 
	 * @param record
	 * @param readSequence: The actual read sequence (not from bam/sam file as this may be reverse complemented)
	 * @param refLength
	 * @return
	 */
	public static Sequence spanningSequence(SAMRecord record, Sequence readSequence, int refLength, int flank){		
		//int flank = 0;
		try{
			int refStart = record.getAlignmentStart();
			int refEnd   = record.getAlignmentEnd();

			int left = (int) (refStart * 1.05) + flank;
			int right = (int) ((refLength - refEnd) * 1.05 + flank);

			int readLength = 0,	readStart = 0, readEnd = 0;	

			boolean enterAlignment = false;		
			for (final CigarElement e : record.getCigar().getCigarElements()) {				
				final int  length = e.getLength();
				switch (e.getOperator()) {
				case H :
				case P : //pad is a kind of clipped
					//throw new RuntimeException("Hard clipping is not supported for this read");
				case S :
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
				readEnd = readLength;//1-index

			if (readLength != readSequence.length()){
				Logging.error("Error0 " + record.getReadName() + " " + readSequence.length() + " vs estimated " + readLength + " Flag = " + record.getFlags());
				return null;
			}

			//start point of the extracted region			
			int start = readStart - left;
			if (start <= 0)
				start = 1;//I am still live in 1-index world

			if (readEnd > readSequence.length()){
				Logging.error("Error1 " + record.getReadName() + " " + record.getReadLength() + " vs " + readEnd);
				return null;
			}
			int end = readEnd + right;

			if (end > readSequence.length())
				end = readSequence.length();

			if (start >= end){
				Logging.error("Error2 " + record.getReadName() + " " + record.getReadLength() + " " + start + " " + end);
				return null;
			}

			if (record.getReadNegativeStrandFlag()){
				//Need to complement the read sequence before calling subsequence = calling sub l-e, l-s then complementing
				//return Alphabet.DNA.complement(readSequence.subSequence(readSequence.length() - end,  readSequence.length() - start + 1));
				Sequence seq = Alphabet.DNA.complement(readSequence).subSequence(start - 1, end);
				seq.setName(readSequence.getName()+"_r_"+readStart+"_"+readEnd + "_"+start+"_"+end);

				return seq;
			}else{
				Sequence seq = readSequence.subSequence(start - 1, end);
				seq.setName(readSequence.getName()+"_"+readStart+"_"+readEnd +"_"+start+"_"+end);
				return seq;
			}

		}catch(Exception e){
			Logging.warn(e.getMessage());
			e.printStackTrace();
			//continue;//while
			return null;
		}

	}
	
	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public static IdentityProfile identity(Sequence refSeq, Sequence readSeq,  SAMRecord sam){
		IdentityProfile profile = new IdentityProfile();
				
		int readPos = 0;//start from 0					
		int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index				
		
		profile.readClipped = 0;
		profile.refClipped = sam.getAlignmentStart() + refSeq.length() - sam.getAlignmentEnd();
		profile.baseDel = 0;
		profile.baseIns = 0;
		profile.numDel = 0;
		profile.numIns = 0;
		profile.match = 0;
		profile.mismatch = 0;
		profile.refBase = 0;
		profile.readBase = 0;//the number of bases from ref and read		

		for (final CigarElement e : sam.getCigar().getCigarElements()) {
			final int  length = e.getLength();
			switch (e.getOperator()) {
			case H :
				//nothing todo
				profile.readClipped += length;
				break; // ignore hard clips
			case P : 
				profile.readClipped += length;
				//pad is a kind of hard clipped ?? 					
				break; // ignore pads	                
			case S :
				//advance on the reference
				profile.readClipped += length;
				readPos += length;
				break; // soft clip read bases	                	
			case N : 
				refPos += length; 
				profile.refClipped += length;
				break;  // reference skip

			case D ://deletion      	
				refPos += length;
				profile.refBase += length;

				profile.baseDel += length;
				profile.numDel ++;
				break; 	

			case I :	                	
				readPos += length;
				profile.readBase += length;

				profile.baseIns += length;
				profile.numIns ++;
				break;
			case M :
				for (int i = 0; i < length; i++){
					if (refSeq.getBase(refPos + i) == readSeq.getBase(readPos + i))
						profile.match ++;
					else
						profile.mismatch ++;
				}
				profile.readBase += length;
				profile.refBase += length;

				readPos += length;
				refPos  += length;
				break;

			case EQ :
				readPos += length;
				refPos  += length;

				profile.readBase += length;
				profile.refBase += length;
				profile.match += length;
				break;

			case X :
				readPos += length;
				refPos  += length;

				profile.readBase += length;
				profile.refBase += length;

				profile.mismatch += length;
				break;
			default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}//case
		}//for			

		return profile;
				
	}
	
	public static class IdentityProfile{
		public int match, mismatch, baseIns, baseDel, numIns, numDel, refClipped, readClipped, refBase, readBase;
		
	}
	
	/**
	 * Compute the N50 statistics of an assembly
	 * @param seqs: List of sequences
	 * @return
	 */
	public static double n50(ArrayList<Sequence> seqs){
		int [] lengths = new int[seqs.size()];
		double sum = 0;
		for (int i = 0;i < lengths.length;i++){
			int l = seqs.get(i).length();
			lengths[i] = l;
			sum += l;
		}		
		Arrays.sort(lengths);
		
		int index = lengths.length;
		double contains = 0;
		while (contains < sum/2){
			index --;
			contains += lengths[index];
		}
		
		return lengths[index];		
	}	
}
