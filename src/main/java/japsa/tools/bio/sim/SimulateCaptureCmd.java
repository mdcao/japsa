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
 * 28/05/2016 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsa.tools.bio.sim;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;


import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.bio.sim.IlluminaSequencing;
import japsa.bio.sim.PacBioSequencing;
import japsa.seq.Genome;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.Simulation;
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

		//Input/output
		addString("reference", null, "Name of genome to be ",true);
		addString("probe", null,  "File containing probes in bam format");
		addString("logFile", "-",  "Log file");
		addString("ID", "",  "A unique ID for the data set");

		addString("fragment", null, "Output of fragment file");
		addString("miseq", null, "Name of read file if miseq is simulated");
		addString("pacbio", null, "Name of read file if pacbio is simulated");

		//Fragment size distribution
		addInt("median", 700 , "Median of fragment size distribution");
		addDouble("shape", 5, "Shape parameter of the fragment size distribution");

		addInt("num", 1000000, "Number of fragments ");		

		//addDouble("mismatch",0.01,"probability of mismatches");
		//addDouble("deletion",0.01,"probability of deletion");
		//addDouble("insertion",0.01,"probability of insertion");
		//addDouble("extension",0.01,"probability of indel extention");


		//Specific parameter for each sequencing technology
		addInt("pblen", 30000, "PacBio: Average (polymerase) read length");

		addInt("illen", 300, "Illumina: read length");
		addString("ilmode", "pe", "Illumina: Sequencing mode: pe = paired-end, mp=mate-paired and se=singled-end");

		addInt("seed", 0, "Random seed, 0 for a random seed");
		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		CommandLine cmdLine = new SimulateCaptureCmd ();
		args = cmdLine.stdParseLine(args);	


		/**********************************************************************/

		String logFile = cmdLine.getStringVal("logFile");
		String probe       =  cmdLine.getStringVal("probe");
		String fragment       =  cmdLine.getStringVal("fragment");
		String ID       =  cmdLine.getStringVal("ID");
		String referenceFile =  cmdLine.getStringVal("reference");

		int median =  cmdLine.getIntVal("median");
		double shape = cmdLine.getDoubleVal("shape");		
		int seed =  cmdLine.getIntVal("seed");		
		int num =   cmdLine.getIntVal("num");

		int pblen = cmdLine.getIntVal("pblen");

		String miseq       =  cmdLine.getStringVal("miseq");
		String pacbio       =  cmdLine.getStringVal("pacbio");

		if (miseq == null && pacbio == null && fragment==null){
			System.err.println("At least one of fragment, miseq and pacbio has to be set\n" + cmdLine.usageString());			
			System.exit(-1);
		}

		int median2 = 400;

		double [] dist2 = new double[median2*4];

		double max = 0.0;
		for (int i = 0; i < dist2.length;i++){
			dist2[i] = Simulation.logLogisticPDF(i + 1, median2, 4);
			if (dist2[i] > max)
				max = dist2[i];			
		}	

		//Normalise to 1
		for (int i = 0; i < dist2.length;i++){			
			//Logging.info("dist2 [" + i + "] = " + dist2[i] + " max = " + max + " after " + dist2[i] / max);
			dist2[i] = dist2[i] / max;
		}

		SequenceOutputStream miSeq1Fq = null, miSeq2Fq = null, pacbioFq = null, sos = null;

		if (miseq != null){
			miSeq1Fq = SequenceOutputStream.makeOutputStream(miseq + "_1.fastq.gz");	
			miSeq2Fq = SequenceOutputStream.makeOutputStream(miseq + "_2.fastq.gz");
		}

		if (pacbio != null){
			pacbioFq = SequenceOutputStream.makeOutputStream(pacbio + ".fastq.gz");
		}	

		//TODO: what is it?
		int flank = median * 4;

		double hybridizationRatio = 0.5;

		SequenceOutputStream 
		logOS =	logFile.equals("-")? (new SequenceOutputStream(System.err)):(SequenceOutputStream.makeOutputStream(logFile));
		logOS.print("Parameters for simulation \n" + cmdLine.optionValues());

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

		//double mu = Math.log(median) + shape * shape;		
		//LogNormalDistribution logNormalDist = new LogNormalDistribution(mu, shape);
		//Reseed the distribution to make sure the same results if the same seed given
		//logNormalDist.reseedRandomGenerator(rnd.nextInt());		

		BitSet [] bitSets = null;
		SamReader samReader =  null;

		if (probe != null){
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			samReader = SamReaderFactory.makeDefault().open(new File(probe));
			SAMRecordIterator samIter = samReader.iterator();
			Logging.info("Mark capturable regions");

			bitSets = new BitSet[chrList.size()];
			for (int i = 0; i < bitSets.length;i++)
				bitSets[i] = new BitSet();

			while (samIter.hasNext()){
				SAMRecord sam = samIter.next();
				if (sam.getReadUnmappedFlag())
					continue;

				int start = sam.getAlignmentStart();			
				int end = sam.getAlignmentEnd();	
				//probe can only bind if > 90%
				if ((end - start) < 0.9 * sam.getReadLength())
					continue;			

				int refIndex = sam.getReferenceIndex();
				int i = Math.max(start - flank,0);
				for (; i < end  && i < chrList.get(refIndex).length();i++){
					bitSets[refIndex].set(i);
				}
			}		
			samIter.close();
			Logging.info("Mark capturable regions -- done");
		}

		if (fragment != null)
			sos = SequenceOutputStream.makeOutputStream(fragment);

		int numFragment = 0;		
		//actual number of fragments generated, including the non probed
		long numGen = 0;

		int 	
		fragmentRej1 = 0,
		fragmentRej2 = 0,
		fragmentRej3 = 0;		

		while (numFragment < num){
			numGen ++;
			if (numGen % 1000000 ==0){
				Logging.info("Generated " + numGen + " selected " + numFragment
						//+ "app = " + numFragmentApp
						+ "; reject1 = " + fragmentRej1 
						+ "; reject2 = " + fragmentRej2
						+ "; reject3 = " + fragmentRej3);				
			}

			//1. Generate the length of the next fragment
			int fragLength = 
					Math.max((int) Simulation.logLogisticSample(median, shape, rnd), 50);	

			//Logging.info("Gen0 " + fragLength);

			//2. Generate the position of the fragment
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

			//Take the min of the frag length and the length to the end
			fragLength = Math.min(fragLength, chrList.get(index).length() - myP);
			if (fragLength < 50){
				Logging.warn("Whoops");
				continue;
			}

			if (samReader != null){
				//if probe is provided, see if the fragment is rejected
				if (!bitSets[index].get(myP)){
					//Logging.info("Reject0 " + fragLength);

					fragmentRej1 ++;
					continue;//while
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
					//probe can only bind if > 80%
					if ((end - start) < hybridizationRatio * sam.getReadLength())
						continue;//while iter probe

					countProbe += (end - start + 1);
					myEnd = end;
				}
				iter.close();

				if (countProbe <= 0){
					//Logging.info("Reject1 " + fragLength + " " + countProbe);
					fragmentRej1 ++;
					continue;
				}


				r = rnd.nextDouble();

				double myOdd = countProbe * 1.0 / fragLength  + 0.5;
				//myOdd = (myOdd - 0.4)/0.6;

				if (r > myOdd){
					//bad luck, rejected
					//Logging.info("Reject2 " + fragLength + " " + countProbe);
					fragmentRej2 ++;
					continue;//while
				}
			}//if

			//Logging.info("Gen1 " + fragLength);
			//here, the fragment is captured

			//another round	of selection
			double myOdd = 0;
			if (fragLength >= dist2.length)
				myOdd = dist2[dist2.length - 1];
			else myOdd = dist2[fragLength];




			if (rnd.nextDouble() > myOdd){
				fragmentRej3 ++;				
				continue;
			}
			numFragment ++;
			//numFragmentApp += count;

			//Logging.info("Gen2 " + fragLength);
			//now that the fragment is to be sequenced
			Sequence seq = chrList.get(index).subSequence(myP, myP + fragLength);
			seq.setName(ID + "_" + chrList.get(index).getName() + "_" + (myP + 1) + "_" +(myP + fragLength));


			if (sos != null)
				seq.writeFasta(sos);

			if (miSeq1Fq != null){
				IlluminaSequencing.simulatePaired(seq, miSeq1Fq, miSeq2Fq, rnd);
			}

			if (pacbioFq != null){
				PacBioSequencing.simulatePacBio(seq, pblen, pacbioFq, rnd);
			}

		}
		if (sos != null)
			sos.close();

		if (miSeq1Fq != null)
			miSeq1Fq.close();

		if (miSeq2Fq != null)
			miSeq2Fq.close();

		if (pacbioFq != null)
			pacbioFq.close();		

		logOS.close();
	}	
	static class GenomicRegion{
		int length, start, chrIndex;
	}	
	
}
