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

/****************************************************************************
 *                           Revision History                                
 * 01/04/2014 - Minh Duc Cao: Started                                 
 *  
 ****************************************************************************/

package japsadev.tools;

import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTextWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


/**
 * 
 */
@Deployable(
	scriptName = "jsa.dev.filterPE", 
	scriptDesc = "Filter concordance PE reads (keeps only reads having mate mapped within a distance)"
)
public class FilterPEConcordance extends CommandLine{ 

	public FilterPEConcordance(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("input", "-", "Name of input sam/bam file (- for from standard input)", true);
		addString("output", "-", "Name of output file, (- for from standard out)");
		addInt("max", 700, "The maxmum size of a fragment");
		
		addStdHelp();	
	}

	static int checkPoint = 2000000;

	public static void main(String[] args) throws Exception {
		CommandLine cmdLine = new FilterPEConcordance();
		args = cmdLine.stdParseLine(args);
		
		
		String output = cmdLine.getStringVal("output");
		String samFile = cmdLine.getStringVal("input");
		filter(samFile, output, cmdLine.getIntVal("max"));
	}

	static void filter(String inFile, String outFile, int max)
			throws IOException {
		HashMap<String, SAMRecord> hashRec = new HashMap<String, SAMRecord>();
		LinkedList<SAMRecord> listRec = new LinkedList<SAMRecord>();

		///////////////////////////////////////////////////////////
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(inFile));						
			
		
		SAMFileHeader samHeader = samReader.getFileHeader();		
		SAMTextWriter bamWriter = new SAMTextWriter(new File(outFile));
		
		//samWriter.setSortOrder(SortOrder.unsorted, false);		
		bamWriter.writeHeader( samHeader.getTextHeader());	
		///////////////////////////////////////////////////////////


		int readIn = 0, readMapped = 0, readOut = 0;

		SAMRecordIterator samIter = samReader.iterator();
		int currentRefIndex = -1;
		int dup = 0;

		while (samIter.hasNext()) {
			SAMRecord sam = samIter.next();
			readIn++;

			if (readIn % 10000000 == 0) {
				Date date = new Date();

				Logging.info("No. of reads processed : " + readIn + " at "
						+ date.toString());
				Logging.info("Statistics so far : " + " Total Read    = "
						+ readIn + " Mapped Read   = " + readMapped
						+ " Out Reads     = " + readOut);
				Logging.info("Hash = " + hashRec.size() + " List = "
						+ listRec.size() + " dup = " + dup + " at "
						+ currentRefIndex);
				Runtime.getRuntime().gc();
			}

			// if the read is not mapped: continue
			if (sam.getReadUnmappedFlag())
				continue;

			readMapped++;
			int order = 0;
			if (sam.getFirstOfPairFlag())
				order = 1;

			int referenceIndex = sam.getReferenceIndex();
			int pos = sam.getAlignmentStart();

			// clean if ref differ to the current one
			if (referenceIndex != currentRefIndex) {
				while (!listRec.isEmpty()) {// clear queue and hash
					SAMRecord head = listRec.remove();
					hashRec.remove(head.getReadName() + "_"
							+ (head.getFirstOfPairFlag() ? 1 : 0));
					if (head.getProperPairFlag()) {
						bamWriter.addAlignment(head);
						readOut++;
					}
				}
				currentRefIndex = referenceIndex;
				dup = 0;
			} else {
				while (!listRec.isEmpty()) {
					SAMRecord head = listRec.getFirst();
					if (pos - head.getAlignmentStart() > max) {
						listRec.remove();
						hashRec.remove(head.getReadName() + "_"
								+ (head.getFirstOfPairFlag() ? 1 : 0));
						if (head.getProperPairFlag()) {
							bamWriter.addAlignment(head);
							readOut++;
						}
					} else
						break;// while
				}// while
			}

			// String refID = sam.getReferenceName();
			String readName = sam.getReadName();
			String readKey = readName + "_" + order;

			SAMRecord mate = hashRec.get(readName + "_" + (1 - order));
			if (mate != null) {
				// found: both reads are properly mapped
				// pair them together
				mate.setMateAlignmentStart(sam.getAlignmentStart());
				sam.setMateAlignmentStart(mate.getAlignmentStart());

				mate.setMateReferenceIndex(referenceIndex);
				sam.setMateReferenceIndex(referenceIndex);

				mate.setMateNegativeStrandFlag(sam.getReadNegativeStrandFlag());
				sam.setMateNegativeStrandFlag(mate.getReadNegativeStrandFlag());

				mate.setMateUnmappedFlag(false);
				sam.setMateUnmappedFlag(false);

				mate.setProperPairFlag(true);
				sam.setProperPairFlag(true);

				mate.setInferredInsertSize(sam.getAlignmentEnd()
						- mate.getAlignmentStart());
				sam.setInferredInsertSize(-mate.getInferredInsertSize());

			} else {
				sam.setProperPairFlag(false);
			}

			listRec.add(sam);
			if (hashRec.put(readKey, sam) != null)
				dup++;
		}// while

		// finally clear
		while (!listRec.isEmpty()) {// clear queue and hash
			SAMRecord head = listRec.remove();
			hashRec.remove(head.getReadName() + "_"
					+ (head.getFirstOfPairFlag() ? 1 : 0));
			if (head.getProperPairFlag()) {
				bamWriter.addAlignment(head);
				readOut++;
			}
		}

		samReader.close();
		bamWriter.close();

		Date date = new Date();

		Logging.info("Finally No. of reads processed : " + readIn + " at "
				+ date.toString());
		Logging.info("Statistics so far : " + " Total Read    = " + readIn
				+ " Mapped Read   = " + readMapped + " Out Reads     = "
				+ readOut);
	}
}

