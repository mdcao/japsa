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
 * 13/07/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.tools;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * A tool to repair concordant pairs information
 * 
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.dev.pairrepair", 
scriptDesc = "Repair concordant pairs from multiple alignment (such as bwa mem -a). This bamfile should be sorted by read name (or from running as a single-threaded program).")
public class PairedEndRepair extends CommandLine{
	public PairedEndRepair(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("bamFile", null, "Name of the s/bam file", true);		
		addString("tr", null, "Name of TR file");
		
		addStdHelp();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new PairedEndRepair();		
		args = cmdLine.stdParseLine(args);	
		
		/*********************** Setting up script ****************************/

		String bamFile = cmdLine.getStringVal("bamFile");
		String trFile = cmdLine.getStringVal("tr");

		if (trFile != null) {
			BufferedReader bf = SequenceReader.openFile(trFile);
			trList = TandemRepeat.readFromFile(bf, null);
		}

		ArrayList<SAMRecord> first = new ArrayList<SAMRecord>(16), second = new ArrayList<SAMRecord>(
				16);

		ArrayList<Integer> disList = new ArrayList<Integer>();
		
		///////////////////////////////////////////////////////////
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFile));						
		

		// SAMFileWriterFactory factory = new SAMFileWriterFactory();
		// SAMFileWriter bamWriter =
		// factory.makeSAMOrBAMWriter(reader.getFileHeader(), true, new
		// File(output));

		String currentID = null;
		SAMRecordIterator iter = samReader.iterator();

		while (iter.hasNext()) {

			SAMRecord rec = iter.next();
			String readName = rec.getReadName();
			// start a new

			if (currentID != null && !readName.equals(currentID)) {
				disList.clear();
				for (SAMRecord f : first) {
					for (SAMRecord s : second) {
						int d = distance(f, s);
						if (d < 2500)
							disList.add(d);
					}// for
				}// for

				if (trList != null) {
					call(currentID, disList, false);
				} else {
					System.out.print(currentID);
					for (int x : disList) {
						System.out.print("\t" + x);
					}// for
					System.out.println();
				}

				first.clear();
				second.clear();
			}// fi

			currentID = readName;

			if (rec.getFirstOfPairFlag())
				first.add(rec);
			else
				second.add(rec);
		}// while

		disList.clear();
		for (SAMRecord f : first) {
			for (SAMRecord s : second) {
				int d = distance(f, s);
				if (d < 2500)
					disList.add(d);
			}// for
		}// for

		if (trList != null) {
			call(currentID, disList, true);// the last call

		} else {
			System.out.print(currentID);
			for (int x : disList) {
				System.out.print("\t" + x);
			}// for
			System.out.println();
		}

		first.clear();
		second.clear();

		samReader.close();
	}

	public static int distance(SAMRecord rec1, SAMRecord rec2) {

		if (rec1.getReadUnmappedFlag() || rec2.getReadUnmappedFlag())
			return Integer.MAX_VALUE;

		if (rec1.getReferenceIndex().intValue() != rec2.getReferenceIndex()
				.intValue())
			return Integer.MAX_VALUE;

		int dis = Math.abs(rec1.getAlignmentStart() - rec2.getAlignmentEnd());

		int d = Math.abs(rec2.getAlignmentStart() - rec1.getAlignmentEnd());
		if (d > dis)
			dis = d;

		return dis + 1;
	}

	static TandemRepeat currentTR = null;

	static ArrayList<TandemRepeat> trList = null;
	static int index = -1;
	static ArrayList<Stat> stats = new ArrayList<Stat>();

	public static void call(String name, ArrayList<Integer> intList,
			boolean last) {
		String[] toks = name.split("#");
		if (index < 0) {
			index = 0;
			currentTR = trList.get(0);
			System.out.print(currentTR.getChr() + "\t" + currentTR.getID()
					+ "\t" + currentTR.getPeriod());
			stats.clear();
		}

		while ((!currentTR.getID().equals(toks[0]))
				|| !currentTR.getChr().equals(toks[2])) {
			index++;

			// Collections.sort(stats);

			for (Stat stat : stats) {
				System.out.printf("\t%.2f,%d", stat.sum / stat.count,
						stat.count);
				// System.out.printf("\t%.2f",stat.sum/stat.count);
			}
			System.out.println();

			if (index >= trList.size()) {
				// return;
				throw new RuntimeException("index " + index
						+ " run out of bound");
			}

			currentTR = trList.get(index);
			System.out.print(currentTR.getChr() + "\t" + currentTR.getID()
					+ "\t" + currentTR.getPeriod());
			stats.clear();
		}

		int e = Integer.parseInt(toks[5]);

		for (int d : intList) {
			double var = ((d - e) * 1.0 / currentTR.getPeriod());
			if (var < 50) {
				// System.out.printf("\t%8.2f", var);
				boolean notuse = true;
				for (Stat stat : stats) {
					// System.out.print(var + "  " + Math.abs(var -
					// stat.sum/stat.count) + " " + stat.sum + " " +
					// stat.count);
					if (Math.abs(var - stat.sum / stat.count) < 0.5) {
						// System.out.println(true);
						stat.sum += var;
						stat.count++;
						notuse = false;
						break;// for
					}// if
						// else
						// System.out.println(false);
				}// for
				if (notuse) {
					stats.add(new Stat(1, var));
				}// if
			}// if
		}// for

		// flush
		if (last) {
			while (true) {
				index++;

				// Collections.sort(stats);
				for (Stat stat : stats) {
					System.out.printf("\t%.2f,%d", stat.sum / stat.count,
							stat.count);
				}
				System.out.println();

				if (index >= trList.size()) {
					return;
				}

				currentTR = trList.get(index);
				System.out.print(currentTR.getChr() + "\t" + currentTR.getID()
						+ "\t" + currentTR.getPeriod());
				stats.clear();
			}
		}
	}

	public static class Stat implements Comparable<Stat> {
		int count = 0;
		double sum = 0;

		Stat(int c, double s) {
			count = c;
			sum = s;
		}

		public String toString() {
			return (sum / count) + "," + count;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(Stat oo) {
			// Stat oo = (Stat) o;
			return oo.count - count;
		}

	}
}
