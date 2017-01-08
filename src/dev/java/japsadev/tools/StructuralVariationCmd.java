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
 * 8 Jan 2017 - Minh Duc Cao: Created                                        
 ****************************************************************************/
package japsadev.tools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 * Analysis of structural variation. The idea is to find alignments of 
 * different fragments on the same reads that contradict each other.
 *
 */
@Deployable(
		scriptName = "jsa.dev.svstream",
		scriptDesc = "Structural variation detection"
		)
public class StructuralVariationCmd extends CommandLine {
	public StructuralVariationCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("input",null,"Sam/bam file (freshed from bwa or sorted by name)",true);
		addString("reference",null,"Reference sequence",true);
		addInt("quality",0,"Minimum alignment quality score");

		addStdHelp();		
	}
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new StructuralVariationCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String input = cmdLine.getStringVal("input");
		String reference = cmdLine.getStringVal("reference");
		int qual = 0;


		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = null;
		if ("-".equals(input))
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			reader = SamReaderFactory.makeDefault().open(new File(input));

		ArrayList<Sequence> refs = SequenceReader.readAll(reference, Alphabet.DNA16());

		SAMRecordIterator iter = reader.iterator();
		String readID = "";		
		ArrayList<AlignmentRec> samList = new ArrayList<AlignmentRec>(); 
		long currentReadCount = 0,
				currentBaseCount = 0;


		while (iter.hasNext()) {
			SAMRecord rec = iter.next();

			if (rec.getReadUnmappedFlag() || rec.getMappingQuality() < qual){
				//not count
				if (!readID.equals(rec.getReadName())){
					samList.clear();
					readID = rec.getReadName();
					currentReadCount ++;
					currentBaseCount += rec.getReadLength();				
				}
				continue;		
			}

			AlignmentRec myRec = new AlignmentRec(rec, refs.get(rec.getReferenceIndex()));

			if (readID.equals(myRec.readID)) {				
			}else {
				if (samList.size() > 1){
					System.out.println(readID);
					for (AlignmentRec aRec:samList){
						System.out.printf("%8d %8d %6s %10d %10d %c %d\n",aRec.readStart, aRec.readEnd,aRec.sequence.getName(),aRec.refStart, aRec.refEnd, aRec.strand?'+':'-',aRec.score);	
					}
				}
				samList.clear();
				readID = myRec.readID;				
				currentReadCount ++;
				currentBaseCount += rec.getReadLength();				
			}
			samList.add(myRec);

		}// while

		reader.close();
	}

	public static class AlignmentRec{
		String readID;
		Sequence sequence;
		int refStart, refEnd, readStart, readEnd;
		public int readLength = 0;
		public boolean strand = true;//positive				
		ArrayList<CigarElement> alignmentCigars = new ArrayList<CigarElement>();
		int score = 0;

		public AlignmentRec(SAMRecord sam, Sequence seq) {
			//		readID = Integer.parseInt(sam.getReadName().split("_")[0]);
			readID = sam.getReadName();

			sequence = seq;

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
			
			score = refEnd + 1 - refStart;
			if (sam.getReadNegativeStrandFlag()){			
				strand = false;
				//need to convert the alignment position on read the correct direction 
				readStart = 1 + readLength - readStart;
				readEnd = 1 + readLength - readEnd;
			}


		}
	}

}
