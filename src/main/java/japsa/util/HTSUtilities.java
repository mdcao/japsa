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
import java.util.List;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A collection of utilities to analyse HTS, based on HTS library
 * @author minhduc
 *
 */
public class HTSUtilities {
	private static final Logger LOG = LoggerFactory.getLogger(HTSUtilities.class);



	/**
	 * Extract read between start and end (on ref)
	 * @param record
	 * @param readSequence
	 * @param fromPos
	 * @param toPos
	 * @return
	 */
	public static Sequence readSequence(SAMRecord record, Sequence readSequence, int fromPos, int toPos){
		int refStart = record.getAlignmentStart();
		int refEnd   = record.getAlignmentEnd();

		if(refStart > fromPos || refEnd < toPos)
			return null;

		int refPos = refStart;
		int readPos = 0;
		int readFrom = 0, readTo = 0;

		for (final CigarElement e : record.getCigar().getCigarElements()){
			if (readTo > 0)
				break;

			final int  length = e.getLength();
			switch (e.getOperator()) {
			case H :
				break; // ignore hard clips
			case P :				 					
				break; // ignore pads	                
			case S :
				readPos += length;
				break; // soft clip read bases	                	
			case N : 
				refPos += length;				
				break;  // reference skip
			case D ://deletion 
				refPos += length;

				if (refPos >=fromPos && readFrom == 0){					
					readFrom = readPos;
				}

				if (refPos >=toPos && readTo == 0){					
					readTo = readPos + 1;
				}

				break;
			case I :	                	
				readPos += length;
				break;
			case M :
			case EQ :				
			case X :				
				readPos += length;
				refPos  += length;

				if (refPos >=fromPos && readFrom == 0){					
					readFrom = readPos - refPos + fromPos;
				}

				if (refPos >=toPos && readTo == 0){			
					readTo = readPos - refPos + toPos;					
				}

				break;
			default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}//case
		}//for	

		if (readFrom ==0 || readTo ==0){
			LOG.error("Error at HTSUtilities.readSequence " + readFrom + " " + readTo, 1);
			System.exit(1);
		}
		if (record.getReadNegativeStrandFlag()){
			//Need to complement the read sequence before calling subsequence = calling sub l-e, l-s then complementing
			//return Alphabet.DNA.complement(readSequence.subSequence(readSequence.length() - end,  readSequence.length() - start + 1));
			Sequence seq = Alphabet.DNA.complement(readSequence).subSequence(readFrom - 1, readTo);
			seq.setName(readSequence.getName()+"_r_"+readFrom+"_"+readTo);

			return seq;
		}else{
			Sequence seq = readSequence.subSequence(readFrom - 1, readTo);
			seq.setName(readSequence.getName()+"_"+readFrom+"_"+readTo);
			return seq;
		}
	}
	
	
	public static Sequence getReadPosition(SAMRecord rec, int startRef, int endRef){
		byte[]  seqRead = rec.getReadBases();//
		if (seqRead.length <= 1)
			return null;

		int startRead = -1, endRead = -1;

		int refPos = rec.getAlignmentStart();
		int readPos = 0;		
		//currentRefPos <= startRead				

		for (final CigarElement e : rec.getCigar().getCigarElements()) {
			int length = e.getLength();
			switch (e.getOperator()) {
			case H:
				break; // ignore hard clips
			case P:
				break; // ignore pads
			case S:
				readPos += e.getLength();								
				break; // soft clip read bases
			case N: // N ~ D
			case D:
				refPos += length;

				if (startRead < 0  && refPos >= startRef){					
					startRead = readPos;
				}

				if (endRead < 0  && refPos >= endRef){					
					endRead = readPos;
				}

				break;// case
			case I:				
				readPos += length;						
				break;

			case M:
			case EQ:
			case X:				
				if ((startRead < 0) && refPos + length >= startRef) {
					startRead = readPos + startRef - refPos;					
				}

				if ((endRead < 0) && (refPos + length >= endRef)){
					endRead = readPos + endRef - refPos;
				}

				refPos += length;
				readPos += length;				
				break;
			default:
				throw new IllegalStateException(
						"Case statement didn't deal with cigar op: "
								+ e.getOperator());
			}// case
			if (refPos >= endRef)
				break;//for

		}// for
		if (startRead < 0 || endRead < 0){
			LOG.warn(" " + refPos + "  " + readPos + " " + startRead + " " + endRead);
			return null;
		}		

		Alphabet alphabet = Alphabet.DNA16();
		Sequence retSeq = new Sequence(alphabet, endRead - startRead + 1, rec.getReadName() + "/" + startRead + "_" + endRead);
		for (int i = 0; i < retSeq.length();i++){
			retSeq.setBase(i, alphabet.byte2index(seqRead[startRead + i]));			
		}
		return retSeq;

	}

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
			//readStart = position of alignment start of read
			//readEnd = position of alignment end of read
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
				LOG.error("Error0 " + record.getReadName() + " " + readSequence.length() + " vs estimated " + readLength + " Flag = " + record.getFlags());
				return null;
			}

			//start point of the extracted region			
			int start = readStart - left;
			if (start <= 0)
				start = 1;//I am still live in 1-index world

			if (readEnd > readSequence.length()){
				LOG.error("Error1 " + record.getReadName() + " " + record.getReadLength() + " vs " + readEnd);
				return null;
			}
			int end = readEnd + right;

			if (end > readSequence.length())
				end = readSequence.length();

			if (start >= end){
				LOG.error("Error2 " + record.getReadName() + " " + record.getReadLength() + " " + start + " " + end);
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
			LOG.warn(e.getMessage());
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
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++){
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
	
	
	public static void  identity(SAMRecord sam, double[] v){
		List<AlignmentBlock> li = sam.getAlignmentBlocks();
		int start = li.get(0).getReadStart();
		int end = li.get(0).getLength()+start;
		int start0 = start;
		int read_gaps=0;
		
		int start_ref = li.get(0).getReferenceStart();
		int start_ref0 = start_ref;
		int end_ref = li.get(0).getLength()+start_ref;
		int ref_gaps=0;
		
		
		for(int i=1; i<li.size(); i++){
			start=li.get(i).getReadStart();
			int gap = li.get(i).getReadStart() - end;
			end = start + li.get(i).getLength();
			read_gaps+=gap;
			
			start_ref=li.get(i).getReferenceStart();
			int gap_ref = li.get(i).getReferenceStart() - end_ref;
			end_ref = start_ref + li.get(i).getLength();
			ref_gaps+=gap_ref;
			
		}
		AlignmentBlock last = li.get(li.size()-1);
		int end0 = last.getReadStart()+last.getLength();
		int end_ref0 = last.getReferenceStart()+last.getLength();
//"read_gap_first\tread_gaps\tread_gap_last\tref_gaps\talignment_blocks\talignment_blocks_rel
		v[0] = (double)start0/(double)sam.getReadLength();
		v[1] = (double) read_gaps/(double)(end0 - start0);
		v[2] = (double) (sam.getReadLength()-end0)/ (double) sam.getReadLength();
		v[3] = (double) ref_gaps/(double)(end_ref0 - start_ref0);
		v[4] = li.size();
		v[5] = (double) li.size()/(double) sam.getReadLength();
		
	}
	
	
	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * Note that this cannot distinguish match and mismatch  since reflen not available
	 * @param seqLen
	 * @param sam
	 * @return
	 */
	public static IdentityProfile identity(int refLen, SAMRecord sam){
		IdentityProfile profile = new IdentityProfile();
		//Sequence readSeq = new Sequence(Alphabet.DNA16(),sam.getReadBases());

		int readPos = 0;//start from 0					
		int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index				

		profile.readClipped = 0;
		profile.refClipped = sam.getAlignmentStart() + refLen - sam.getAlignmentEnd();
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
				for (int i = 0; i < length && refPos + i < refLen; i++){
					//if (refSeq.getBase(refPos + i) == readSeq.getBase(readPos + i))
						profile.match ++;
				//	else
					//	profile.mismatch ++;
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


	/**
	 * Get the list of positions in reads corresponding the the positions in reference
     * @param sam
	 * @param refPositions
	 * @return
	 */
	public static int[] positionsInRead(SAMRecord sam, final int [] refPositions){
		int readPos = 0;//start from 0					
		int refPos = sam.getAlignmentStart();//convert to 0-based index		
		int [] readPositions = new int[refPositions.length];
		int index = 0;

		while (index < refPositions.length && refPositions[index] <= refPos)
			index ++;

		if (index >= refPositions.length) 
			return readPositions;		

		for (final CigarElement e : sam.getCigar().getCigarElements()) {
			//assert: refPositions[index] > refPos

			final int  length = e.getLength();
			switch (e.getOperator()) {
			case H :
				//nothing todo				
				break; // ignore hard clips
			case P :				
				//pad is a kind of hard clipped ?? 					
				break; // ignore pads	                
			case S :
				//soft clip: advance on the reference				
				readPos += length;
				break; // soft clip read bases	                	
			case N : 
				refPos += length;

				//advance index
				while (index < refPositions.length && refPositions[index] <= refPos)
					index ++;

				if (index >= refPositions.length) 
					return readPositions;				
				break;  // reference skip

			case D ://deletion      	
				refPos += length;
				while (index < refPositions.length && refPositions[index] <= refPos){
					readPositions[index] = readPos;					
					index ++;					
				}
				if (index >= refPositions.length) 
					return readPositions;				
				break;
			case I :	                	
				readPos += length;
				break;
			case M :
			case EQ:
			case X:
				while (index < refPositions.length && refPositions[index] <= refPos +length){
					readPositions[index] = readPos + refPositions[index] - refPos;
					index ++;					
				}
				if (index >= refPositions.length) 
					return readPositions;		

				readPos += length;
				refPos  += length;
				break;
			default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}//case
		}//for			

		return readPositions;

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
	
	public static double n50(ArrayList<Sequence> seqs, long genomeSize){
		int [] lengths = new int[seqs.size()];
		
		for (int i = 0;i < lengths.length;i++){
			int l = seqs.get(i).length();
			lengths[i] = l;		
		}		
		Arrays.sort(lengths);

		int index = lengths.length;
		double contains = 0;
		while (contains < genomeSize/2){
			index --;
			contains += lengths[index];
		}

		return lengths[index];		
	}	

}
