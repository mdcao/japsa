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
 * 28/05/2014 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsa.tools.bio.sim;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Genome;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;


/**
 * @author minhduc
 *
 */
@Deployable(
		scriptName = "jsa.sim.capture", 
		scriptDesc = "Simulate the capture process"
		)
public class SimulateCaptureCmd extends CommandLine{
	//CommandLine cmdLine;
	public SimulateCaptureCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());		

		addString("probe", "probes.bam",  "File containing probes in bam format");				
		addString("output", "fragment.fasta",  "Output of fragment file");
		addString("logFile", "-",  "Log file");
		addString("reference", "hg19.fas", "Name of reference genome");

		addInt("mean", 800, "Fragment size mean");
		addInt("std", 80, "Fragment size standard deviation");
		addInt("num", 1000000, "Number of fragments ");
		addInt("seed", 0, "Random seed, 0 for a random seed");

		addString("miseq", null, "Name of read file if miseq is simulated");
		addString("pacbio", null, "Name of read file if pacbio is simulated");

		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		SimulateCaptureCmd cmdLine = new SimulateCaptureCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/

		String logFile = cmdLine.getStringVal("logFile");
		String probe       =  cmdLine.getStringVal("probe");
		String output       =  cmdLine.getStringVal("output");
		String referenceFile =  cmdLine.getStringVal("reference");

		int mean = cmdLine.getIntVal("mean");
		int std = cmdLine.getIntVal("std");		
		int seed = cmdLine.getIntVal("seed");		
		int num = cmdLine.getIntVal("num");

		String miseq       =  cmdLine.getStringVal("miseq");
		String pacbio       =  cmdLine.getStringVal("pacbio");

		SequenceOutputStream miSeq1Fq = null, miSeq2Fq = null, pacbioFq = null;

		if (miseq != null){
			miSeq1Fq = SequenceOutputStream.makeOutputStream(miseq + "_1.fastq");	
			miSeq2Fq = SequenceOutputStream.makeOutputStream(miseq + "_2.fastq");
		}

		if (pacbio != null){
			pacbioFq = SequenceOutputStream.makeOutputStream(pacbio + ".fastq");
		}	



		SequenceOutputStream 
		logOS 	=	logFile.equals("-")? (new SequenceOutputStream(System.err))	:(SequenceOutputStream.makeOutputStream(logFile));

		seed = SimulateGenomeCmd.seed(seed);
		Random rnd = new Random(seed);

		logOS.print("#Seed " + seed + "\n");


		Genome genome = new Genome();		
		genome.read(referenceFile);
		ArrayList<Sequence> chrList = genome.chrList();

		logOS.print("Read " + chrList.size() + " chr " + genome.getLength() + "bp\n" );				
		long [] accLen = new long[chrList.size()];		
		accLen[0] = chrList.get(0).length();		
		Logging.info("Acc 0 " + accLen[0]);		
		for (int i = 1; i < accLen.length;i++){
			accLen[i] = accLen[i-1] + chrList.get(i).length();
			Logging.info("Acc " +i + " " +  accLen[i]);			
		}		

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(probe));						


		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);	
		int numFragment = 0;

		//actual number of fragments generated, including the non probed
		long numGen = 0;

		while (numFragment < num){

			numGen ++;
			if (numGen % 1000000 ==0){
				Logging.info("Generated " + numGen + " selected " + numFragment);
			}

			//toss the coin
			double r = rnd.nextDouble();
			long p = (long) (r * genome.getLength());			
			int index = 0;
			while (p > accLen[index])
				index ++;
			//identify the chroms
			if (index > 0){
				p = p - accLen[index - 1];
			}

			int myP = (int) p;
			//Logging.info("Found " + index + " " + myP + "  " + p);

			if (p > chrList.get(index).length()){
				Logging.exit("Not expecting " + p + " vs " + index, 1);
			}

			int fragLength = (int) (rnd.nextGaussian() * std + mean);
			if (myP + fragLength > chrList.get(index).length()){
				Logging.warn("Whoops");
				continue;
			}

			SAMRecordIterator iter = samReader.query(chrList.get(index).getName(), myP, myP + fragLength, true);

			int countProbe = 0;
			int myEnd = 0;
			while (iter.hasNext()){				
				SAMRecord sam = iter.next();
				
				int start = sam.getAlignmentStart();
				if (start < myEnd)
					continue;
				
				int end = sam.getAlignmentEnd();				
				//probe can only bind if > 90%
				if ((end - start) < 0.9 * sam.getReadLength())
					continue;
				
				countProbe += (end - start + 1);
				myEnd = end;				
			}
			iter.close();

			if (countProbe <= 0)
				continue;							 

			r = rnd.nextDouble();
			if (r < countProbe * 1.0 / fragLength){
				Sequence seq = chrList.get(index).subSequence(myP, myP + fragLength);
				seq.setName(chrList.get(index).getName() + "_" + (myP + 1) + "_" +(myP + fragLength));
				seq.writeFasta(sos);
				numFragment ++;		

				if (miSeq1Fq != null){
					simulatePaired(seq, miSeq1Fq, miSeq2Fq, rnd);
				}
			}
		}
		sos.close();


		if (miSeq1Fq != null)
			miSeq1Fq.close();

		if (miSeq2Fq != null)
			miSeq2Fq.close();

		if (pacbioFq != null)
			pacbioFq.close();		

		logOS.close();


	}

	/**
	 * Simulate a read from the start of a fragment
	 * @param fragment
	 * @param len
	 * @param snp
	 * @param indel
	 * @param ext
	 * @param rnd
	 * @return
	 */
	static Sequence simulateRead(Sequence fragment, int len,  double snp, double indel, double ext, Random rnd){
		Sequence read = new Sequence(fragment.alphabet(), len);
		int mIndex = 0, fIndex = 0;
		for (;mIndex < len && fIndex < fragment.length();){
			byte base = fragment.getBase(fIndex);
			if (rnd.nextDouble() < snp){
				read.setBase(mIndex, (byte) ((base + rnd.nextInt(3)) % 4));
				mIndex ++;
				fIndex ++;				
			}else if (rnd.nextDouble() < indel){
				if (rnd.nextDouble() < 0.5){
					//deletion
					do{
						fIndex ++;
					}while (rnd.nextDouble() < ext);
				}else{
					//insertion
					do{
						mIndex ++;
					}while (rnd.nextDouble() < ext);
				}//else					
			}else{
				read.setBase(mIndex, base);
				mIndex ++;
				fIndex ++;
			}//else
		}//for
		for (;mIndex < len;mIndex ++){
			//pad in random to fill in
			read.setBase(mIndex, (byte) rnd.nextInt(4));
		}
		return read;
	}

	/**
	 * Simulate MiSeq
	 * @param fragment
	 * @param o1
	 * @param o2
	 * @param rnd
	 * @throws IOException
	 */
	static void simulatePaired(Sequence fragment, SequenceOutputStream o1 , SequenceOutputStream o2, Random rnd) throws IOException{		

		double snp = 0.01;
		double indel = 0.0001;				
		double ext = 0.2;

		//double r = 
		int len = Math.min(250, fragment.length());
		String name = fragment.getName();				

		Sequence read1 = simulateRead(fragment, len, snp, indel, ext, rnd);
		Sequence read2 = simulateRead(Alphabet.DNA.complement(fragment), len, snp, indel, ext, rnd);

		if (rnd.nextBoolean()){
			o1.print("@");
			o1.print(name);						
			o1.print("\n");
			for (int i = 0; i < read1.length();i++)
				o1.print(read1.charAt(i));
			o1.print("\n+\n");

			for (int i = 0; i < read1.length();i++)
				o1.print("I");
			o1.print("\n");

			o2.print("@");
			o2.print(name);						
			o2.print("\n");
			for (int i = 0; i < read2.length();i++)
				o2.print(read2.charAt(i));
			o2.print("\n+\n");

			for (int i = 0; i < read2.length();i++)
				o2.print("I");
			o2.print("\n");			
		}else{
			o2.print("@");
			o2.print(name);						
			o2.print("\n");
			for (int i = 0; i < read1.length();i++)
				o2.print(read1.charAt(i));
			o2.print("\n+\n");

			for (int i = 0; i < read1.length();i++)
				o2.print("I");
			o2.print("\n");

			o1.print("@");
			o1.print(name);						
			o1.print("\n");
			for (int i = 0; i < read2.length();i++)
				o1.print(read2.charAt(i));
			o1.print("\n+\n");

			for (int i = 0; i < read2.length();i++)
				o1.print("I");
			o1.print("\n");			
		}
	}

	static void simulatePacBio(Sequence fragment, SequenceOutputStream o, Random rnd) throws IOException{
		double snp = 0.01;
		double indel = 0.1;				
		double ext = 0.4;

		int len = (int) (fragment.length() * .9);
		String name = fragment.getName();		
		Sequence read	=  (rnd.nextBoolean())?simulateRead(fragment, len, snp, indel, ext, rnd): simulateRead(Alphabet.DNA.complement(fragment), len, snp, indel, ext, rnd);

		o.print("@");
		o.print(name);						
		o.print("\n");
		for (int i = 0; i < read.length();i++)
			o.print(read.charAt(i));
		o.print("\n+\n");

		for (int i = 0; i < read.length();i++)
			o.print("E");
		o.print("\n");

	}
}
