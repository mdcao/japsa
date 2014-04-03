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
 * 08/03/2012 - Minh Duc Cao: Started
 * 03/05/2012 - Minh Duc Cao: 
 *                 - Standardise commandLine
 * 16/11/2013 - Minh Duc Cao                                   
 *  
 ****************************************************************************/

package japsa.bio.tr;


import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;



/**
 * FIXME: Need to test
 * @author minhduc Standardise SAM format, the data from the following are
 *         supported: - JGIHeazlewood2011: (1001 Arabidopsis Genomes)
 * 
 */
@Deployable(scriptName = "jsa.trv.sam2fragment",
            scriptDesc = "Convert a sam file to list of fragment sizes")
public class Sam2FragmentSize {
	static int checkPoint = 2000000;
	public static void main(String[] args) throws Exception {	

		/*********************** Setting up script ****************************/		 
		String scriptName = "jsa.trv.sam2fragment";
		String desc = "Convert a sam file to list of fragment sizes\n";		
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [options]");
		/**********************************************************************/

		cmdLine.addString("input", "-", "Name of input sam/bam file (- for from standard input)",true);
		cmdLine.addString("output", "-", "Name of output file, (- for from standard out)");
		cmdLine.addBoolean("recalc", false, "Recalculate fragment size instead of using the default");

		cmdLine.addInt("min",100,"The minimum size of a fragment");
		cmdLine.addInt("max",1000,"The maxmum size of a fragment");		

		cmdLine.addBoolean("2", true, "Whether filter out reads with flag 0x002 turned off");

		cmdLine.addStdHelp();


		/**********************************************************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")){
			System.out.println(desc + cmdLine.usage());			
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usage());
			System.exit(-1);
		}	
		/**********************************************************************/		

		String output = cmdLine.getStringVal("output");
		String samFile = cmdLine.getStringVal("input");

		//If there are errors
		if (args == null) {
			System.err.println(cmdLine.usage());
			System.exit(-1);
		}

		final long last = System.currentTimeMillis();


		//System.err.println("in = "  + samFile + " out = " + output + cmdLine.getBooleanVal("2"));		
		//getInsertSizeModerate(samFile, output, cmdLine.getBooleanVal("2"),cmdLine.getIntVal("min"), cmdLine.getIntVal("max"), cmdLine.getBooleanVal("recalc"));		
		getInsertSizeModerate2(samFile, output, cmdLine.getBooleanVal("2"),cmdLine.getIntVal("min"), cmdLine.getIntVal("max"), cmdLine.getBooleanVal("recalc"));


		final long now = System.currentTimeMillis();
		System.err.println("Stop after  " + (now - last) / 1000.0 + " seconds)");
	}

	/****************************************
	@Deprecated
	static void getInsertSizeModerate(String inFile, String outFile, boolean f2, int min, int max, boolean recalc)
	throws IOException {

		int lineNo = 0;		

		if (min < 0) min = 0;//cant have any length less than 0

		HashMap<String, PEFragment> hashPair = new HashMap<String, PEFragment>();

		BufferedReader in = null;
		if ("-".equals(inFile))
			in = new BufferedReader(new InputStreamReader(System.in));
		else
			in = SequenceReader.openFile(inFile);

		BufferedOutputStream ps;
		if (outFile.equals("-"))
			ps = new BufferedOutputStream(System.out);
		else
			ps = new BufferedOutputStream(new FileOutputStream(outFile));

		String line = "";
		String[] toks;

		int totalRead = 0, pairedRead = 0, pairs = 0, outPair = 0;

		int numF = 0, numI = 0, maxF = 0, maxI = 0, minF = 10000, minI = 10000;
		double sumF = 0, sumFSQ = 0, sumI = 0, sumISQ = 0;

		while ((line = in.readLine()) != null) {
			lineNo++;
			// headers
			if (line.startsWith("@")){
				ps.write((line+'\n').getBytes());
				continue;
			}

			totalRead++;

			if (totalRead % checkPoint == 0) {
				Date date = new Date();
				System.err.println("No. of reads processed : " + totalRead
						+ " at " + date.toString());
				System.err.println("Statistics so far : "
						+ "\n Total Read    = " + totalRead
						+ "\n Paired Read   = " + pairedRead
						+ "\n Printed pairs = " + outPair
						+ "\n No Inserts    = " + pairs  
						+ "\n Hash Contains = " + hashPair.size());
				double mean = sumF/numF;

				System.err.printf("#Stat Fragment: min = %4d max = %4d n = %d mean = %6.2f std=%6.3f\n",minF, maxF, numF, mean, Math.sqrt((sumFSQ - numF * mean * mean)/numF));

				mean = sumI/numI;
				System.err.printf("#Stat Insert: min = %4d max = %4d n = %d mean = %6.2f std=%6.3f\n",minI, maxI, numI, mean, Math.sqrt((sumISQ - numI * mean * mean)/numI));

				System.err.println("##"+totalRead+"#"+pairedRead+"#"+pairs+"#"+outPair);

				Runtime.getRuntime().gc();
			}
			toks = //IOTools.TAB_SEPARATOR.split(line.trim());
			line.trim().split("\\t");
			if (toks.length < 11) {
				System.err.println("ERROR: SAM format corrupted at line "
						+ lineNo + ":\n  " + line + "\nOnly " + toks.length
						+ " fields");
				//continue;
				System.exit(-1);
			}

			int flags = Integer.parseInt(toks[1]);// flags

			//Make sure the read has a mate and both of them are mapped
			if (((flags & 1) == 0) || ((flags & 8) != 0) || ((flags & 4) != 0 ) ) {
				continue;
			}

			//Also filter out not proper map			
			if (f2 && ((flags & 2) == 0) ) {
				continue;
			}

			pairedRead++;
			String refID = toks[2];			
			int startPos = Integer.parseInt(toks[3]);
			int iSize = Integer.parseInt(toks[8]);
			int readLen = toks[9].length();

			if (iSize > max || iSize < -max)
				continue;
			if (iSize > -min && iSize < min)
				continue;

			// Check if its mate is in the queue
			PEFragment pair = hashPair.remove(toks[0] + "#" + refID);

			if (pair != null) {// already in the queue
				// get more info for the read
				pair.start2 = startPos;
				pair.len2 = readLen;

				// do a sanity check
				if (pair.iSize + iSize != 0) {
					System.err.println("Warning  at line " + lineNo + ":\n  "
							+ line + "\n Ref = " + refID
							+ " insert size not match : " + pair.iSize + " vs " + iSize);
				}

				pair.fragmentBegin = pair.start2;
				if (pair.fragmentBegin > pair.start1) 
					pair.fragmentBegin = pair.start1;

				pair.fragmentEnd = pair.start2 + pair.len2;
				if (pair.fragmentEnd < pair.start1 + pair.len1) 
					pair.fragmentEnd = pair.start1 + pair.len1;

				pair.fragmentEnd --;//inclusive				
				// 

				pair.length = pair.fragmentEnd +1 - pair.fragmentBegin; 
				if (recalc)
					pair.iSize = pair.length;

				pair.println(ps);
				numF ++;

				int tmp = (pair.fragmentEnd - pair.fragmentBegin) + 1; 
				sumF += tmp;
				sumFSQ += tmp * tmp;

				if (tmp > maxF) 
					maxF = tmp;
				if (tmp < minF) 
					minF = tmp;

				if (pair.iSize != 0){
					numI ++;
					tmp = Math.abs(pair.iSize);
					sumI += tmp;
					sumISQ += tmp * tmp;

					if (tmp > maxI) maxI = tmp;
					if (tmp < minI) minI = tmp;
				}		

				outPair ++;
			} else {// not yet in the queue
				pair = new PEFragment(refID, toks[0]);
				pair.start1 = startPos;
				pair.len1 = readLen;

				pair.iSize = iSize;
				pair.direction = (flags & 0x30) >> 4;// both

				// Put in the hash and first of the list
				hashPair.put(toks[0]+"#"+refID, pair);
				pairs++;
			}

		}// while

		if (!hashPair.isEmpty()){
			System.err.println("WARN: There are " + hashPair.size() + " pairs left");
			System.err.println ("#####"+hashPair.values().iterator().next().readID);
		}

		ps.write (("#param : min = " + min + " max = " +max+"\n").getBytes());
		double mean = sumF/numF;

		ps.write(("#Fragment: min = " + minF + " max = " + maxF + " = " + numF + " mean = " + mean + " std=" + (Math.sqrt((sumFSQ - numF * mean * mean)/numF)) + "\n\n").getBytes());

		mean = sumI/numI;
		ps.write(("#Insert: min = " + minI + " max = " + maxI + " = " + numI + " mean = " + mean + " std=" + (Math.sqrt((sumISQ - numI * mean * mean)/numI)) + "\n\n").getBytes());		
		ps.write(("##"+totalRead+"#"+pairedRead+"#"+pairs+"#"+outPair+"\n").getBytes());


		in.close();		
		ps.close();

		Date date = new Date();

		System.err.println("Finally No. of reads processed : " + totalRead
				+ " at " + date.toString());
		System.err.println("Statistics so far : " + 
				"\n Total Read    = " 	+ totalRead + 
				"\n Paired Read   = " + pairedRead +				 
				"\n Printed pairs = " + outPair
				+ "\n No Inserts    = " + pairs);

	}

	/****************************************/

	static void getInsertSizeModerate2(String inFile, String outFile, boolean f2, int min, int max, boolean recalc)
			throws IOException {					
		int lineNo = 0;					
		if (min < 0) min = 0;//cant have any length less than 0

		HashMap<String, PEFragment> hashPair = new HashMap<String, PEFragment>();

		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader samReader = new  SAMFileReader(new File(inFile));




		SequenceOutputStream ps = SequenceOutputStream.makeOutputStream(outFile);		
		SAMFileHeader samHeader = samReader.getFileHeader();
		ps.write(samHeader.getTextHeader());

		int totalRead = 0, pairedRead = 0, pairs = 0, outPair = 0;				
		int numF = 0, numI = 0, maxF = 0, maxI = 0, minF = 10000, minI = 10000;
		double sumF = 0, sumFSQ = 0, sumI = 0, sumISQ = 0;

		SAMRecordIterator samIter = samReader.iterator();				

		while (samIter.hasNext()) {
			SAMRecord sam = samIter.next();
			lineNo++;
			totalRead++;

			if (totalRead % checkPoint == 0) {
				Date date = new Date();
				System.err.println("No. of reads processed : " + totalRead
						+ " at " + date.toString());
				System.err.println("Statistics so far : "
						+ "\n Total Read    = " + totalRead
						+ "\n Paired Read   = " + pairedRead
						+ "\n Printed pairs = " + outPair
						+ "\n No Inserts    = " + pairs  
						+ "\n Hash Contains = " + hashPair.size());
				double mean = sumF/numF;

				System.err.printf("#Stat Fragment: min = %4d max = %4d n = %d mean = %6.2f std=%6.3f\n",minF, maxF, numF, mean, Math.sqrt((sumFSQ - numF * mean * mean)/numF));

				mean = sumI/numI;
				System.err.printf("#Stat Insert: min = %4d max = %4d n = %d mean = %6.2f std=%6.3f\n",minI, maxI, numI, mean, Math.sqrt((sumISQ - numI * mean * mean)/numI));

				System.err.println("##"+totalRead+"#"+pairedRead+"#"+pairs+"#"+outPair);

				Runtime.getRuntime().gc();
			}					

			int flags = sam.getFlags();				

			//Make sure the read has a mate and both of them are mapped
			if (((flags & 1) == 0) || ((flags & 8) != 0) || ((flags & 4) != 0 ) ) {
				continue;
			}

			//Also filter out not proper map			
			if (f2 && ((flags & 2) == 0) ) {
				continue;
			}					

			//filter out: not primary alignment,  
			//            read fails platform/vendor quality checks
			//            read is PCR or optical duplicate
			//            supplementary alignment
			if ((flags & 3840) != 0)
				continue;

			pairedRead++;					

			//String refID = sam.getReferenceName();
			String readName = sam.getReadName();

			int refIndex = sam.getReferenceIndex();

			int startPos = sam.getAlignmentStart();					
			//							sam.getReferencePositionAtReadPosition(1);

			int iSize = sam.getInferredInsertSize();
			int readLen = sam.getReadLength();

			//if (iSize > max || iSize < -max)
			//	continue;
			//if (iSize > -min && iSize < min)
			//	continue;

			// Check if its mate is in the queue
			PEFragment pair = hashPair.remove(sam.getReadName() + "#" + refIndex);

			if (pair != null) {// already in the queue						
				//check if duplicated
				if (pair.iSize == iSize && pair.start1 == startPos){						
					hashPair.put(readName+"#"+refIndex, pair);
					System.err.println("Put back " + pair.readID);	
					continue;
				}					

				// get more info for the read
				pair.start2 = startPos;
				pair.len2 = readLen;

				// do a sanity check
				if (pair.iSize + iSize != 0) {
					System.err.println("Warning  at read " + lineNo + ":\n  "
							+ "\n Ref = " + sam.getReferenceName()
							+ " ID = " + sam.getReadName()
							+ " insert size not match : " + pair.iSize + " vs " + iSize);
				}

				pair.fragmentBegin = pair.start2;
				if (pair.fragmentBegin > pair.start1) 
					pair.fragmentBegin = pair.start1;

				pair.fragmentEnd = pair.start2 + pair.len2;
				if (pair.fragmentEnd < pair.start1 + pair.len1) 
					pair.fragmentEnd = pair.start1 + pair.len1;

				pair.fragmentEnd --;//inclusive				
				// 

				pair.length = pair.fragmentEnd +1 - pair.fragmentBegin; 

				if (recalc)
					pair.iSize = pair.length;

				if (pair.iSize < 0)
					pair.iSize = - pair.iSize;

				if (pair.iSize > max || pair.iSize < min)
					continue;

				pair.println(ps);
				numF ++;

				int tmp = (pair.fragmentEnd - pair.fragmentBegin) + 1; 
				sumF += tmp;
				sumFSQ += tmp * tmp;

				if (tmp > maxF) 
					maxF = tmp;
				if (tmp < minF) 
					minF = tmp;

				if (pair.iSize != 0){
					numI ++;
					tmp = Math.abs(pair.iSize);
					sumI += tmp;
					sumISQ += tmp * tmp;

					if (tmp > maxI) maxI = tmp;
					if (tmp < minI) minI = tmp;
				}		

				outPair ++;
			} else {// not yet in the queue
				pair = new PEFragment(refIndex, readName );
				pair.start1 = startPos;
				pair.len1 = readLen;

				pair.iSize = iSize;
				pair.direction = (flags & 0x30) >> 4;// both

					// Put in the hash and first of the list
					hashPair.put(readName+"#"+refIndex, pair);
					pairs++;
			}

		}// while

		if (!hashPair.isEmpty()){
			System.err.println("WARN: There are " + hashPair.size() + " pairs left");
			System.err.println ("#####"+hashPair.values().iterator().next().readID);
		}

		ps.write ("#param : min = " + min + " max = " +max+"\n");
		double mean = sumF/numF;

		ps.write("#Fragment: min = " + minF + " max = " + maxF + " = " + numF + " mean = " + mean + " std=" + (Math.sqrt((sumFSQ - numF * mean * mean)/numF)) + "\n\n");

		mean = sumI/numI;
		ps.write("#Insert: min = " + minI + " max = " + maxI + " = " + numI + " mean = " + mean + " std=" + (Math.sqrt((sumISQ - numI * mean * mean)/numI)) + "\n\n");		
		ps.write("##"+totalRead+"#"+pairedRead+"#"+pairs+"#"+outPair+"\n");


		samReader.close();		
		ps.close();

		Date date = new Date();

		System.err.println("Finally No. of reads processed : " + totalRead
				+ " at " + date.toString());
		System.err.println("Statistics so far : " + 
				"\n Total Read    = " 	+ totalRead + 
				"\n Paired Read   = " + pairedRead +				 
				"\n Printed pairs = " + outPair
				+ "\n No Inserts    = " + pairs);

	}

}


