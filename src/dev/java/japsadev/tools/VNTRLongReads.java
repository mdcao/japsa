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
 * 15/05/2014 - Minh Duc Cao: Started
 *  
 ****************************************************************************/
package japsadev.tools;

import japsa.bio.alignment.ProfileDP;
import japsa.bio.alignment.ProfileDP.EmissionState;
import japsa.bio.tr.TandemRepeat;
import japsa.bio.tr.TandemRepeatVariant;

import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.SequenceOutputStream;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.seq.XAFReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * VNTR typing using long reads
 * 
 */

@Deployable(scriptName = "jsa.dev.longreads_oldversion", scriptDesc = "VNTR typing using long reads, old version. Pls check jsa.tr.longreads")
public class VNTRLongReads {
	static Alphabet dna = Alphabet.DNA16();

	public static void main(String[] args) throws Exception,
	InterruptedException {
		/*********************** Setting up script ****************************/
		Deployable annotation = VNTRLongReads.class
				.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options]",
				annotation.scriptDesc());

		cmdLine.addStdInputFile();
		cmdLine.addString("bamFile", null, "Name of the bam file", true);
		cmdLine.addString("output", "-",
				"Name of the output file, -  for stdout");
		cmdLine.addString("xafFile", null, "Name of the regions file in xaf",
				true);
		cmdLine.addString("prefix", "",
				"Prefix of temporary files, if not specified, will be automatically generated");
		cmdLine.addInt("flanking", 30, "Size of the flanking regions");
		cmdLine.addInt("qual", 0, "Minimum quality");
		cmdLine.addInt(
				"nploidy",
				2,
				"The ploidy of the genome 1 =  happloid, 2 = diploid. Currenly only support up to 2-ploidy");

		cmdLine.addString("msa", "kalign",
				"Name of the msa method, support kalign and clustalo");
		// cmdLine.addBoolean("contained", false,
		// "true: Reads contained in the region; false: reads overlap with the region");

		args = cmdLine.stdParseLine_old(args);
		/**********************************************************************/
		// Get options

		String prefix = cmdLine.getStringVal("prefix");

		if (prefix == null || prefix.length() == 0) {
			prefix = "p" + System.currentTimeMillis();
		}

		int flanking = cmdLine.getIntVal("flanking");
		if (flanking < 20)
			flanking = 20;

		int qual = cmdLine.getIntVal("qual");

		int np = cmdLine.getIntVal("nploidy");
		if (np > 2) {
			System.err
			.println("The program currenly only support haploid and diployd. Enter nploidy of 1 or 2");
			System.exit(1);
		}
		/**********************************************************************/

		String cmd;
		if ("kalign".equals(cmdLine.getStringVal("msa"))) {
			cmd = "kalign -gpo 60 -gpe 10 -tgpe 0 -bonus 0 -q -i " + prefix
					+ "i.fasta -o " + prefix + "o.fasta";
		} else {
			cmd = "clustalo --force -i " + prefix + "i.fasta -o " + prefix
					+ "o.fasta";
		}

		SequenceOutputStream outOS 
		= SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("output"));

		String[] headers = TandemRepeatVariant.SIMPLE_HEADERS;
		if (np > 1) {
			headers = TandemRepeatVariant.SIMPLE_HEADERS2;
		}

		TandemRepeatVariant.printHeader(outOS, headers);

		String strFile = cmdLine.getStringVal("xafFile");

		// TODO: make it multiple sequence
		//Sequence seq = SequenceReader.getReader(cmdLine.getStringVal("input"))
		//		.nextSequence(dna);

		Logging.info("Read genome begins");
		HashMap <String, Sequence> genome = new HashMap <String, Sequence>();
		SequenceReader seqReader = SequenceReader.getReader(cmdLine.getStringVal("input"));
		Sequence seq;
		while ((seq = seqReader.nextSequence(dna)) != null){
			genome.put(seq.getName(), seq);
		}
		seqReader.close();
		Logging.info("Read genome done");

		/**********************************************************************/

		XAFReader xafReader = new XAFReader(strFile);

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(cmdLine.getStringVal("bamFile")));	


		//Random random = new Random();

		//int _tIndex = 0;
		while (xafReader.next() != null){	
			//_tIndex ++;
			TandemRepeat str = TandemRepeat.read(xafReader);

			int start = Integer.parseInt(xafReader.getField("start")) - flanking;			
			int end = Integer.parseInt(xafReader.getField("end")) + flanking;
			String chrom = xafReader.getField("chrom");					

			if (seq == null || !seq.getName().equals(chrom)){
				seq = genome.get(chrom);
			}
			if (seq == null){
				xafReader.close();
				//sos.close();
				Logging.exit("Chrom in line " + xafReader.lineNo() + " not found!!!", 1);
			}

			if (end > seq.length())
				end = seq.length();

			if (start < 1)
				start = 1;

			int hmmFlank = flanking;
			int hmmPad = 
					(int)((str.getUnitNo() - Math.floor(str.getUnitNo())) * str.getPeriod()) ;

			//System.out.println("###" + str.getPeriod() + " " + str.getUnitNo() + "   " + hmmPad);			
			Sequence hmmSeq = new Sequence(dna, hmmFlank * 2 + hmmPad + str.getPeriod());
			int i = 0;

			for (;i < hmmFlank + hmmPad + str.getPeriod(); i++)
				hmmSeq.setBase(i, seq.getBase(str.getStart() - hmmFlank + i -1));

			for (;i < hmmSeq.length();i++){
				byte base = seq.getBase(str.getEnd() + i - (hmmFlank + hmmPad + str.getPeriod()) );//no need to -1
				hmmSeq.setBase(i,base);				
			}

			ProfileDP dp = new ProfileDP(hmmSeq, hmmFlank + hmmPad, hmmFlank + hmmPad + str.getPeriod() - 1);//-1 for 0-index, inclusive

			//System.out.println("Lengths: " + hmmFlank + ", " + hmmPad + " " + str.getPeriod() + " " + hmmSeq.length() );			
			//System.out.println("CHECKING BEGIN");

			outOS.print("##"+str.getID()+"\n## ");
			for (int x = 0; x < hmmSeq.length();x++){
				outOS.print(hmmSeq.charAt(x));
				if (x == hmmFlank + hmmPad -1 || x ==  hmmFlank + hmmPad + str.getPeriod() - 1)
					outOS.print("==");
			}
			outOS.println();			
			//for (int x = 0; x < 10; x++ ){			
			//	Sequence genSeq = dp.generate(random.nextInt(5) + 5);
			//	double alignScore = dp.align(genSeq);
			//	System.out.println("Best myCost: " + alignScore + " vs " + genSeq.length() * 2 + " (" + alignScore/genSeq.length() +")");
			//}
			//System.out.println("CHECKING END");

			SAMRecordIterator iter 
			= reader.query(str.getParent(), start, end, false);


			int maxAlign = 300;

			MultipleAlignment ma = new MultipleAlignment(maxAlign, seq);
			while (iter.hasNext()) {
				SAMRecord rec = iter.next();
				// Check qualilty
				if (rec.getMappingQuality() < qual) {
					continue;
				}

				// Only reads that fully span the repeat and flankings
				if (rec.getAlignmentStart() > start)
					continue;
				if (rec.getAlignmentEnd() < end)
					continue;

				ma.addRead(rec);
			}// while
			iter.close();
			// os.close();

			double var = 0;
			TandemRepeatVariant trVar = new TandemRepeatVariant();
			trVar.setTandemRepeat(str);

			if (ma.printFasta(start, end, prefix + "i.fasta") > 0) {
				Logging.info("Running " + cmd);
				Process process = Runtime.getRuntime().exec(cmd);
				process.waitFor();
				Logging.info("Done " + cmd);

				SequenceReader hmmSeqReader 
				= FastaReader.getReader(prefix	+ "i.fasta");
				Sequence readSeq;

				//IntArray intArray = new IntArray();
				//DoubleArray doubleArray = new DoubleArray();
				//ByteArray byteArray = new ByteArray(); 

				while ((readSeq = hmmSeqReader.nextSequence(dna)) != null ){					
					//System.out.println(readSeq.getName() + " : " + readSeq.length());

					EmissionState bestState = dp.align(readSeq);
					double alignScore = bestState.getScore();
					//System.out.println("Score " + alignScore + " vs " + readSeq.length()*2 + " (" + alignScore/readSeq.length() +")");
					int bestIter = bestState.getIter();

					outOS.print("==================================================================\n");
					outOS.print("##" + readSeq.getName()+"\t"+bestIter+"\t"+readSeq.length() +"\t" +alignScore+"\t" + alignScore/readSeq.length() + '\n');

					/*********************************************************
					intArray.clear();
					doubleArray.clear();
					byteArray.clear();


					//double oldCost = bestState.score;
					Emission lastState = bestState;

					bestState = bestState.bwdState;

					while (bestState != null){
						intArray.add(bestState.profilePos);						
						doubleArray.add(lastState.score - bestState.score);
						if (bestState.seqPos == lastState.seqPos)
							byteArray.add((byte)Alphabet.DNA.GAP);
						else
							byteArray.add(readSeq.getBase(lastState.seqPos));						

						lastState = bestState;
						bestState = bestState.bwdState;
					}					

					double costL = 0, costR = 0;

					for (int x = intArray.size() - 1; x >=0; x--){
						outOS.print(Alphabet.DNA().int2char(byteArray.get(x)));						
						 if (x <intArray.size() - 1  && intArray.get(x) <  + intArray.get(x+1)){
							 outOS.println();	 
						 }

						if (intArray.get(x) < hmmFlank + hmmPad)
							costL += doubleArray.get(x);
						if (intArray.get(x) > hmmFlank + hmmPad + str.getPeriod())
							costR += doubleArray.get(x);													
					}
					outOS.println();
					outOS.print ("L = " + (costL/(hmmFlank + hmmPad)) + " R = " + costR/(hmmSeq.length() - hmmFlank - hmmPad - str.getPeriod()) + "\n");
					outOS.print("==================================================================\n");
					/*********************************************************/
				}

				SequenceReader msaReader
				= FastaReader.getReader(prefix	+ "o.fasta");
				ArrayList<Sequence> seqList = new ArrayList<Sequence>();
				Sequence nSeq = null;
				while ((nSeq = msaReader.nextSequence(Alphabet.DNA16())) != null) {
					seqList.add(nSeq);
				}				 
				//str.getChr()+"_"+str.getStart()+"_"+str.getEnd();

				if (np >= 2) {
					trVar.setVar2(var);
					trVar.addEvidence(seqList.size());
				} else {// nploidy ==1
					int llength = seqList.get(0).length();

					int gaps = call(seqList);

					var = (llength - gaps - end + start) * 1.0
							/ str.getPeriod();

					trVar.setVar(var);
					trVar.addEvidence(seqList.size());
				}
			}// if

			//Process process = 
			//		Runtime.getRuntime().exec("rm -f " + prefix + _tIndex + "i.fasta");
			//process.waitFor();

			//process = 
			//	Runtime.getRuntime().exec("cp " + prefix + "i.fasta " + prefix + _tIndex + "i.fasta");
			//process.waitFor();			

			outOS.print(trVar.toString(headers));
			outOS.print('\n');
		}// for

		reader.close();
		outOS.close();
	}
	/**
	 * 
	 * @param seqList
	 * @param start
	 *            : the start index of the list (inclusive)
	 * @param end
	 *            : the end index of the list (exclusive)
	 */
	static int call(ArrayList<Sequence> seqList, int indexStart, int indexEnd) {
		if (indexEnd <= indexStart)
			return 0;	

		// Get consensus
		int gaps = 0;
		Sequence nSeq = new Sequence(Alphabet.DNA6(), seqList.get(0).length(),
				"consensus");
		int[] votes = new int[6];
		for (int i = 0; i < nSeq.length(); i++) {			
			Arrays.fill(votes, 0);
			for (int s = indexStart; s < indexEnd; s++) {
				votes[seqList.get(s).symbolAt(i)]++;
			}
			byte best = 0;
			for (byte b = 1; b < 6; b++)
				if (votes[b] > votes[best])
					best = b;

			nSeq.setBase(i, best);
			if (best == 5)
				gaps++;
		}// for
		return gaps;
	}

	static int call(ArrayList<Sequence> seqList) {
		return call(seqList,0,seqList.size());

	}
	/************************************************
	static int[][] distance(ArrayList<Sequence> seqList) {
		int[][] dis = new int[seqList.size()][seqList.size()];

		for (int i = seqList.size() - 1; i > 1; i--) {
			for (int j = i - 1; j >= 0; j--) {
				dis[i][j] = dis[j][i] = distance(seqList.get(i), seqList.get(j));
			}
		}

		return dis;
	}

	static int distance(Sequence s1, Sequence s2) {
		int dis = 0;
		for (int i = 0; i < s1.length(); i++) {
			if (s1.getBase(i) != s2.getBase(i)) {
				if (s1.getBase(i) == DNA.GAP || s2.getBase(i) == DNA.GAP)
					dis += 1;
				else
					dis += 0;
			}
		}
		return dis;
	}
/************************************************/
}
