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
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.Simulation;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * @author minhduc
 *
 */
@Deployable(
		scriptName = "jsa.sim.capsim", 
		scriptDesc = "Simulate capture sequencing"
		)
public class SimulateCaptureCmd extends CommandLine{
    private static final Logger LOG = LoggerFactory.getLogger(SimulateCaptureCmd.class);

    //CommandLine cmdLine;
	public SimulateCaptureCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		//Input/output
		addString("reference", null, "Name of genome to be ",true);
		addString("probe", null,  "File containing probes mapped to the reference in bam format");
		addString("logFile", "-",  "Log file");
		addString("ID", "",  "A unique ID for the data set");

		//addString("fragment", null, "Output of fragment file");
		addString("miseq", null, "Name of read file if miseq is simulated");
		addString("pacbio", null, "Name of read file if pacbio is simulated");

		//Fragment size distribution
		addInt("fmedian", 2000 , "Median of fragment size at shearing");
		addDouble("fshape", 6, "Shape parameter of the fragment size distribution");

		//addInt("smedian", 1300 , "Median of fragment size distribution");
		//addDouble("sshape", 6, "Shape parameter of the fragment size distribution");

		//addInt("tmedian", 0 , "Median of target fragment size (the fragment size of the data).\n If specified, " +
		//		"will override fmedian and smedian.\n Othersise will be estimated");
		//addDouble("tshape", 0, "Shape parameter of the effective fragment size distribution");

		addInt("num", 1000000, "Number of fragments ");

		//addDouble("mismatch",0.01,"probability of mismatches");
		//addDouble("deletion",0.01,"probability of deletion");
		//addDouble("insertion",0.01,"probability of insertion");
		//addDouble("extension",0.01,"probability of indel extention");
		//Specific parameter for each sequencing technology

		addInt("pblen", 30000, "PacBio: Average (polymerase) read length");

		addInt("illen", 300, "Illumina: read length");
		//addString("ilmode", "pe", "Illumina: Sequencing mode: pe = paired-end, mp=mate-paired and se=singled-end");

		addInt("seed", 0, "Random seed, 0 for a random seed");
		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		CommandLine cmdLine = new SimulateCaptureCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String logFile = cmdLine.getStringVal("logFile");
		String probe       =  cmdLine.getStringVal("probe");
		String ID       =  cmdLine.getStringVal("ID");
		String referenceFile =  cmdLine.getStringVal("reference");

		//int smedian =  cmdLine.getIntVal("smedian");
		//double sshape = cmdLine.getDoubleVal("sshape");

		//int tmedian =  cmdLine.getIntVal("tmedian");
		//double tshape = cmdLine.getDoubleVal("tshape");

		int seed =  cmdLine.getIntVal("seed");		
		int num =   cmdLine.getIntVal("num");

		int pblen = cmdLine.getIntVal("pblen");
		int pbshape = 6;

		String miseq       =  cmdLine.getStringVal("miseq");
		String pacbio      =  cmdLine.getStringVal("pacbio");

		if (miseq == null && pacbio == null){
			System.err.println("One of miseq or pacbio must be set\n" + cmdLine.usageString());			
			System.exit(-1);
		}

		if (miseq != null && pacbio != null){
			System.err.println("Only one of miseq or pacbio can be set\n" + cmdLine.usageString());
			System.exit(-1);
		}

		int fmedian =  cmdLine.getIntVal("fmedian");
		double fshape = cmdLine.getDoubleVal("fshape");

		int flank = fmedian  + (fmedian / 4);
		double hybridizationRatio = 0.5;

		//double [] dist2 = null;
		//if (tmedian <=0){
		//	LOG.info("Estimating target distribution");
		//	dist2 = new double[fmedian*6];
		//}else{
		//	dist2 = new double[tmedian*6];
		//	double max = 0;
		//	for (int i = 0; i < dist2.length;i++){
		//		dist2[i] = Simulation.logLogisticPDF(i + 1, smedian, sshape);
		//		if (dist2[i] > max)
		//			max = dist2[i];
		//	}
		//}

		//double [] dist2 = new double[smedian*4];
		//double max = 0.0;
		//for (int i = 0; i < dist2.length;i++){
		//	dist2[i] = Simulation.logLogisticPDF(i + 1, smedian, sshape);
		//	if (dist2[i] > max)
		//		max = dist2[i];			
		//}

		//Normalise to 1
		//for (int i = 0; i < dist2.length;i++){			
		//	//LOG.info("dist2 [" + i + "] = " + dist2[i] + " max = " + max + " after " + dist2[i] / max);
		//	dist2[i] = dist2[i] / max;
		//}

		SequenceOutputStream miSeq1Fq = null, miSeq2Fq = null, pacbioFq = null;

		if (miseq != null){
			miSeq1Fq = SequenceOutputStream.makeOutputStream(miseq + "_1.fastq.gz");	
			miSeq2Fq = SequenceOutputStream.makeOutputStream(miseq + "_2.fastq.gz");
		}

		if (pacbio != null){
			pacbioFq = SequenceOutputStream.makeOutputStream(pacbio + ".fastq.gz");
		}	

		SequenceOutputStream 
		logOS =	logFile.equals("-")? (new SequenceOutputStream(System.err)):(SequenceOutputStream.makeOutputStream(logFile));
		logOS.print("Parameters for simulation \n" + cmdLine.optionValues());

		seed = Simulation.seed(seed);
		Random rnd = new Random(seed);

		logOS.print("#Seed " + seed + "\n");

		Genome genome = new Genome();		
		genome.read(referenceFile);
		ArrayList<Sequence> chrList = genome.chrList();

		logOS.print("Read " + chrList.size() + " chr " + genome.getLength() + "bp\n" );

		BitSet [] bitSets = null;
		GenomicRegion genRegion = null;
		SamReader samReader =  null;
		long [] accLen = null;

		if (probe != null){
			genRegion = new GenomicRegion();

			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			samReader = SamReaderFactory.makeDefault().open(new File(probe));
			SAMRecordIterator samIter = samReader.iterator();
			LOG.info("Mark capturable regions");

			bitSets = new BitSet[chrList.size()];
			for (int i = 0; i < bitSets.length;i++)
				bitSets[i] = new BitSet();

			//Mark regions from which fragments *may* be captured
			while (samIter.hasNext()){
				SAMRecord sam = samIter.next();
				if (sam.getReadUnmappedFlag())
					continue;

				int start = sam.getAlignmentStart();			
				int end = sam.getAlignmentEnd();

				if ((end - start) < hybridizationRatio * sam.getReadLength())
					continue;

				int refIndex = sam.getReferenceIndex();
				//TODO: this part can be improved
				bitSets[refIndex].set(Math.max(start - flank,0), end);
			}		
			samIter.close();
			//LOG.info("Mark capturable regions -- done");
			for (int x=0; x < genome.chrList().size();x++){
				Sequence chrom = genome.chrList().get(x);
				BitSet myBitSet = bitSets[x];
				int regionStart = -1;				
				for (int i = 0; i < chrom.length();i++){
					if (myBitSet.get(i) && (regionStart < 0)){
						//start of a new region
						regionStart = i;						
					}else if (!myBitSet.get(i) && (regionStart >= 0)){
						//end of a region
						genRegion.addRegion(x, regionStart, i - regionStart);
						regionStart = -1;
					}
				}
				if (regionStart >=0){
					genRegion.addRegion(x, regionStart, chrom.length() - regionStart);
				}
			}//for
			//LOG.info("Mark capturable regions 2 -- done");
		}else{
			accLen = new long[chrList.size()];		
			accLen[0] = chrList.get(0).length();		
			LOG.info("Acc 0 " + accLen[0]);
			for (int i = 1; i < accLen.length;i++){
				accLen[i] = accLen[i-1] + chrList.get(i).length();
				LOG.info("Acc " +i + " " +  accLen[i]);
			}	

		}

		//if (fragment != null)
		//	sos = SequenceOutputStream.makeOutputStream(fragment);

		long numFragment = 0;		
		//actual number of fragments generated, including the non probed
		long numGen = 0;

		long
		fragmentRej1 = 0,
		fragmentRej2 = 0,
		fragmentRej3 = 0,
		fragmentRej4 = 0;		

		while (numFragment < num){
			numGen ++;
			if (numGen % 1000000 == 0){
				LOG.info("Generated " + numGen + " selected " + numFragment
						+ "; reject1 = " + fragmentRej1 
						+ "; reject2 = " + fragmentRej2
						+ "; reject3 = " + fragmentRej3
						+ "; reject4 = " + fragmentRej4);				
			}

			//1. Generate the length of the next fragment
			int fragLength = 
					Math.max((int) Simulation.logLogisticSample(fmedian, fshape, rnd), 50);	

			//LOG.info("Gen0 " + fragLength);

			//2. Generate the position of the fragment
			//toss the coin
			double r = rnd.nextDouble();

			int chrIndex = 0, chrPos = 0;

			if (genRegion != null){				
				long p = (long) (r * genRegion.totLength);
				int index = 0;

				while (p > genRegion.regions.get(index).accuLength){
					index ++;
				}

				if (index > 0){
					p = p - genRegion.regions.get(index - 1).accuLength;
				}

				if (p > genRegion.regions.get(index).length){
					LOG.error("Not expecting2 " + p + " vs " + index);
					System.exit(1);
				}				

				chrPos = ((int) p) + genRegion.regions.get(index).position;
				chrIndex = genRegion.regions.get(index).chrIndex;				
			}else{

				long p = (long) (r * genome.getLength());			
				int index = 0;
				while (p > accLen[index])
					index ++;
				//identify the chroms
				if (index > 0){
					p = p - accLen[index - 1];
				}			

				chrIndex = index;
				chrPos = (int) p;
				//LOG.info("Found " + index + " " + myP + "  " + p);
			}
			if (chrPos > chrList.get(chrIndex).length()){
				LOG.error("Not expecting " + chrPos + " vs " + chrIndex);
				System.exit(1);
			}			

			//Take the min of the frag length and the length to the end
			fragLength = Math.min(fragLength, chrList.get(chrIndex).length() - chrPos);
			if (fragLength < 50){
				LOG.warn("Whoops");
				continue;
			}

			if (samReader != null){
				//if probe is provided, see if the fragment is rejected
				if (!bitSets[chrIndex].get(chrPos)){
					//LOG.info("Reject0 " + fragLength);

					fragmentRej1 ++;
					continue;//while
				}

				/*******************************************************************************
				SAMRecordIterator iter = samReader.query(chrList.get(chrIndex).getName(), chrPos, chrPos + fragLength, false);
				int countProbe = 0;
				int myEnd = 0;
				while (iter.hasNext()){				
					SAMRecord sam = iter.next();

					int start = sam.getAlignmentStart();
					if (start < myEnd)
						continue;

					int end = sam.getAlignmentEnd();				
					//probe can only bind if > 80%
					//if ((end - start) < hybridizationRatio * sam.getReadLength())
					//	continue;//while iter probe

					if (start < chrPos)
						start = chrPos;

					if (end > chrPos + fragLength)
						end = chrPos + fragLength;

					countProbe += (end - start + 1);
					myEnd = end;
				}
				iter.close();

				if (countProbe <= 0){
					//LOG.info("Reject1 " + fragLength + " " + countProbe);
					fragmentRej2 ++;
					continue;
				}
				double myOdd = countProbe * 4.0/ fragLength;
				//myOdd = 4 * (myOdd - 0.1) / 0.9;				
				/*******************************************************************************/

				/*******************************************************************************/				
				SAMRecordIterator iter = samReader.query(chrList.get(chrIndex).getName(), chrPos, chrPos + fragLength, false);
				int countProbe = 0;

				while (iter.hasNext()){				
					SAMRecord sam = iter.next();

					int start = sam.getAlignmentStart();
					int end = sam.getAlignmentEnd();

					if (start < chrPos)
						start = chrPos;

					if (end > chrPos + fragLength)
						end = chrPos + fragLength;

					countProbe += (end - start + 1);					
				}
				iter.close();

				if (countProbe <= 0){
					//LOG.info("Reject1 " + fragLength + " " + countProbe);
					fragmentRej2 ++;
					continue;
				}
				double myOdd = (countProbe * 0.5 / fragLength) - 0.2;								
				/*******************************************************************************/

				r = rnd.nextDouble();
				if (r > myOdd){
					//bad luck, rejected
					//LOG.info("Reject2 " + fragLength + " " + countProbe);
					fragmentRej3 ++;
					//System.out.println("Rejected 3 " + myOdd + " vs " + r);
					continue;//while
				}
			}//if

			//LOG.info("Gen1 " + fragLength);
			//here, the fragment is captured

			//another round	of selection
			//double myOdd = 0;
			//if (fragLength >= dist2.length)
			//	myOdd = dist2[dist2.length - 1];
			//else myOdd = dist2[fragLength];


			//if (rnd.nextDouble() > myOdd){
			//	fragmentRej4 ++;				
			//	continue;
			//}
			numFragment ++;
			//numFragmentApp += count;

			//LOG.info("Gen2 " + fragLength);
			//now that the fragment is to be sequenced
			Sequence seq = chrList.get(chrIndex).subSequence(chrPos, chrPos + fragLength);
			seq.setName(ID + "_" + chrList.get(chrIndex).getName() + "_" + (chrPos + 1) + "_" +(chrPos + fragLength));

			//if (sos != null)
			//	seq.writeFasta(sos);

			if (miSeq1Fq != null){
				IlluminaSequencing.simulatePaired(seq, miSeq1Fq, miSeq2Fq, rnd);
			}

			if (pacbioFq != null){
				int readLen = Math.max((int) Simulation.logLogisticSample(pblen, pbshape, rnd), 50);
				PacBioSequencing.simulatePacBio(seq, readLen, pacbioFq, rnd);
			}

		}

		LOG.info("Generated " + numGen + " selected " + numFragment
				+ "; reject1 = " + fragmentRej1 
				+ "; reject2 = " + fragmentRej2
				+ "; reject3 = " + fragmentRej3
				+ "; reject4 = " + fragmentRej4);				

		//if (sos != null)
		//	sos.close();

		if (miSeq1Fq != null)
			miSeq1Fq.close();

		if (miSeq2Fq != null)
			miSeq2Fq.close();

		if (pacbioFq != null)
			pacbioFq.close();		

		logOS.close();
	}

	/**
	 * Implement regions that may be capturable
	 * @author minhduc
	 *
	 */
	static class GenomicRegion{
		Genome genome;
		long totLength = 0;
		static class Region{
			int chrIndex;
			int position;
			int length;
			long accuLength;
		}

		ArrayList<Region> regions = new ArrayList<Region>();

		Region addRegion(int cIndex, int pos, int length){
			Region region = new Region();
			region.chrIndex = cIndex;
			region.position = pos;
			region.length = length;
			totLength += length;			
			region.accuLength = totLength;

			regions.add(region);			
			return region;
		}
	}

}


/*RST*
----------------------------------------------------------------------------
*capsim*: Simulating the Dynamics of Targeted Capture Sequencing with CapSim
----------------------------------------------------------------------------

*capsim* (jsa.sim.capsim) is a tool to simulate target capture sequencing. Its
simulates the dynamics of capture process

<usage>

~~~~~~~~~~~~~
Usage samples
~~~~~~~~~~~~~

*RST*/
