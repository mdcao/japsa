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
 * 30/06/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.hts.scaffold;

import japsa.seq.Alphabet;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

import java.io.IOException;




/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.dev.gapcloser", scriptDesc = "Gap closer")
public class GapCloser {
	static Alphabet alphabet = Alphabet.DNA();
	//static boolean hardClip = false;

	public static void main(String[] args) throws IOException,
	InterruptedException {
		/*********************** Setting up script ****************************/
		Deployable annotation = GapCloser.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options]",
				annotation.scriptDesc());

		cmdLine.addString("bamFile", null, "Name of the bam file", true);
		cmdLine.addString("sequenceFile", null, "Name of the assembly file (sorted by length)",true);
		cmdLine.addString("output", "-",
				"Name of the output file, -  for stdout");		
		cmdLine.addInt("threshold", 500, "Threshold");
		cmdLine.addDouble("cov", 0, "Expected average coverage of Illumina, <=0 to estimate");
		cmdLine.addInt("qual", 5, "Minimum quality");

		//cmdLine.addBoolean("hardclip", false,"True for last, false for bwa");
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		String output = cmdLine.getStringVal("output");
		String bamFile = cmdLine.getStringVal("bamFile");
		String sequenceFile = cmdLine.getStringVal("sequenceFile");
		//int stage = cmdLine.getIntVal("stage");
		int threshold = cmdLine.getIntVal("threshold");
		double cov = cmdLine.getDoubleVal("cov");
		int qual = cmdLine.getIntVal("qual");
		/**********************************************************************/

		ScaffoldGraph graph = new ScaffoldGraph(sequenceFile);
		if (cov <=0)
			cov = graph.estimatedCov;

		graph.makeConnections(bamFile, cov, threshold, qual);
		graph.connectBridges();
		graph.viewStatus();

		SequenceOutputStream outOS = SequenceOutputStream.makeOutputStream(output);
		graph.printScaffoldSequence(outOS);
		outOS.close();		
	}
	
	public static class AlignmentRecord {
		static final double matchCost = 0;

		public String name;
		public int refIndex;

		public int refStart,  //position on ref of the start of the alignment
		refEnd;    //position on ref at the start of the
		public int refLength;

		//Position on read of the start and end of the alignment 
		public int readStart = 0, readEnd = 0;
		
		
		//read length
		public int readLength = 0;		
		
		public boolean strand = true;//positive
		SAMRecord sam;
		public boolean useful = false;

		
		public int readLeft, readRight, readAlign, refLeft, refRight, refAlign;

		public AlignmentRecord(SAMRecord sam, int rLength) {
			name = sam.getReadName();
			refIndex = sam.getReferenceIndex();

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


			this.sam = sam;
			this.refLength = rLength;

			readLeft = readStart;
			readRight = readLength - readEnd;
			readAlign = readEnd + 1 - readStart;

			refLeft = refStart;
			refRight = refLength - refEnd;
			refAlign = refEnd + 1 - refStart;

			if (sam.getReadNegativeStrandFlag()){
				strand = false;
				readStart = readLength - readStart;
				readEnd = readLength - readEnd;
			}

			//only useful if
			if ((readLeft > refLeft + 250 || readRight > 250 + refRight)
					&& (readLeft < 250 || refLeft < 250)
					&& (readRight  < 250 || refRight < 250)
					)
				useful = true;
		}

		public String toString() {
			return refIndex  
					//			+ " " + pStart 
					//			+ " " + pAlignStart 
					//			+ " " + pAlignEnd 
					//			+ " " + pEnd  
					+ " " + refStart 
					+ " " + refEnd 
					//					+ " " + leftClipped 
					//					+ " " + rightClipped 
					//					+ " " + (leftClipped + rightClipped + pAlignEnd - pAlignStart)
					+ " " + refLength
					//+ " " + (pAlignEnd - pAlignStart)
					+ " " + strand;
		}
		public String pos() {			
			return  
					refStart 
					+ " " + refEnd
					+ " " + refLength
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
	}	
}
