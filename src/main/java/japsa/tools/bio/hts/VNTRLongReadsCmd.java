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
package japsa.tools.bio.hts;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.clustering.Cluster;
import org.apache.commons.math3.stat.clustering.KMeansPlusPlusClusterer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.bio.alignment.MultipleAlignment;
import japsa.bio.alignment.ProfileDP;
import japsa.bio.alignment.ProfileDP.EmissionState;
import japsa.bio.tr.ReadAllele;
import japsa.bio.tr.RepeatCluster;
import japsa.bio.tr.TandemRepeat;
import japsa.bio.tr.TandemRepeatVariant;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.seq.XAFReader;
import japsa.util.ByteArray;
import japsa.util.CommandLine;
import japsa.util.DoubleArray;
import japsa.util.IntArray;
import japsa.util.JapsaMath;
import japsa.util.deploy.Deployable;
import japsa.xm.expert.Expert;
import japsa.xm.expert.MarkovExpert;



/**
 * VNTR typing using long reads
 * 
 */

@Deployable(scriptName = "jsa.tr.longreads", scriptDesc = "VNTR typing using long reads")
public class VNTRLongReadsCmd  extends CommandLine {
	private static final Logger LOG = LoggerFactory.getLogger(VNTRLongReadsCmd.class);
	public VNTRLongReadsCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		CommandLine.Option referenceOpt =
				addString("reference", null, "Name of the reference genome ", true);
		///addStdInputFile();
		CommandLine.Option bamFileOpt = 
				addString("bamFile", null, "Name of the bam file", true);

		CommandLine.Option outputOpt =
				addString("output", "-",
						"Name of the output file, -  for stdout");

		CommandLine.Option xafFileOpt =
				addString("xafFile", null, "Name of the regions file in xaf",
						true);

		CommandLine.Option flankingOpt =
				addInt("flanking", 30, "Size of the flanking regions");

		CommandLine.Option minQualOpt =
				addInt("qual", 0, "Minimum quality");

		addInt("iteration", 1, "Number of iteration");
		addInt("nploidy",2,
				"The ploidy of the genome 1 =  happloid, 2 = diploid. Currenly only support up to 2-ploidy");
		addString("prefix", "",
				"Prefix of temporary files, if not specified, will be automatically generated");

		///////////////Adding galaxy support/////////////
		flankingOpt.setGalaxySetting(new GalaxySetting("integer", null,false));
		minQualOpt.setGalaxySetting(new GalaxySetting("integer", null,false));
		xafFileOpt.setGalaxySetting(new GalaxySetting("data", "tabular",false));

		GalaxySetting galaxyOutput = new GalaxySetting("data", "text",true); 
		galaxyOutput.setLabel("countRead.txt");
		outputOpt.setGalaxySetting(galaxyOutput);
		bamFileOpt.setGalaxySetting(new GalaxySetting("data", "bam",false));
		referenceOpt.setGalaxySetting(new GalaxySetting("data", "fasta",false));		
		setGalaxy(annotation.scriptName());	


		addStdHelp();		
	} 

	static Alphabet dna = Alphabet.DNA16();
	static IntArray profilePositions = new IntArray();
	static IntArray seqPositions = new IntArray();		
	static DoubleArray costGeneration = new DoubleArray();
	static ByteArray byteArray = new ByteArray();


	public static void main(String[] args) throws Exception,
	InterruptedException {
		/*********************** Setting up script ****************************/
		CommandLine cmdLine = new VNTRLongReadsCmd();
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		int flanking = cmdLine.getIntVal("flanking");
		if (flanking < 10)
			flanking = 10;

		int qual = cmdLine.getIntVal("qual");

		int np = cmdLine.getIntVal("nploidy");
		if (np > 2) {
			System.err.println("The program currenly only support haploid and diployd. Enter nploidy of 1 or 2");
			System.exit(1);
		}

		String bamFile = cmdLine.getStringVal("bamFile");
		String prefix = cmdLine.getStringVal("prefix");

		if (prefix == null || prefix.length() == 0) {
			prefix = "p" + System.currentTimeMillis();
		}
		/**********************************************************************/

		SequenceOutputStream outOS = SequenceOutputStream
				.makeOutputStream(cmdLine.getStringVal("output"));

		String[] headers = TandemRepeatVariant.SIMPLE_HEADERS;
		if (np > 1) {
			headers = TandemRepeatVariant.SIMPLE_HEADERS2;
		}

		TandemRepeatVariant.printHeader(outOS, headers);

		String strFile = cmdLine.getStringVal("xafFile");

		LOG.info("Read genome begins");
		HashMap <String, Sequence> genome = new HashMap <String, Sequence>();
		SequenceReader seqReader = SequenceReader.getReader(cmdLine.getStringVal("reference"));
		Sequence seq;
		while ((seq = seqReader.nextSequence(dna)) != null){
			genome.put(seq.getName(), seq);
		}
		seqReader.close();
		LOG.info("Read genome done");

		/**********************************************************************/
		XAFReader xafReader = new XAFReader(strFile);

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));						

		Expert.setAlphabet(Alphabet.DNA4());

		ArrayList<Sequence> readSequences = new ArrayList<Sequence>(); 

		//Go through the list of repeats
		while (xafReader.next() != null){
			List<ReadAllele> accumulator = new ArrayList<ReadAllele>();
			List<ReadAllele> accumulator_batch = new ArrayList<ReadAllele>();
			TandemRepeat str = TandemRepeat.read(xafReader);

			//start,end = the start and end of the region (including flanks)
			int start = Integer.parseInt(xafReader.getField("start")) - flanking;			
			int end = Integer.parseInt(xafReader.getField("end")) + flanking;
			String chrom = xafReader.getField("chrom");					

			if (seq == null || !seq.getName().equals(chrom)){
				seq = genome.get(chrom);
			}
			if (seq == null){
				xafReader.close();				
				LOG.error("Chrom in line " + xafReader.lineNo() + " not found!!!");
				System.exit(1);
			}

			if (end > seq.length())
				end = seq.length();

			if (start < 1)
				start = 1;

			int hmmFlank = flanking;
			int period = str.getPeriod();
			double fraction = str.getUnitNo() - Math.floor(str.getUnitNo());			
			int hmmPad = (int)(fraction * period ) ;

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
				if (x == hmmFlank + hmmPad -1 || x ==  hmmFlank + hmmPad + period - 1)
					outOS.print("==");
			}
			outOS.println();			

			//run on the reference
			//if (1==0)
			{
				Sequence refRepeat = seq.subSequence(start, end);
				refRepeat.setName("reference");
				processRead(refRepeat, dp, fraction,  hmmFlank, hmmPad, period,  outOS , accumulator);

			}

			SAMRecordIterator iter = reader.query(str.getParent(), start, end, false);

			String fileName = prefix + "_" + str.getID() + "_i.fasta";
			SequenceOutputStream os = SequenceOutputStream.makeOutputStream(fileName);

			//			double var = 0;
			TandemRepeatVariant trVar = new TandemRepeatVariant();
			trVar.setTandemRepeat(str);

			int readIndex = 0;

			//Clear the list
			readSequences.clear();
			while (iter.hasNext()) {
				SAMRecord rec = iter.next();
				// Check qualilty
				if (rec.getMappingQuality() < qual) {
					continue;
				}

				// Only reads that fully span the repeat and flankings
				int currentRefPos = rec.getAlignmentStart();
				if (currentRefPos > start)
					continue;
				if (rec.getAlignmentEnd() < end)
					continue;

				readIndex ++;
				////////////////////////////////////////////////////////////////////
				//assert currentRefBase < start

				Sequence readSeq = getReadPosition(rec,start,end);
				if (readSeq == null)
					continue;


				String readName = readSeq.getName();
				String [] toks =  readName.split("/",4);

				String polymerageRead = (toks.length > 1) ? toks[1] : toks[0];
				String subRead = (toks.length > 2) ? toks[2] : "_";
				//String alignSubRead = (toks.length > 3) ? toks[3] : "_";
				readSeq.setName(polymerageRead + "_" + subRead);				
				readSeq.writeFasta(os);

				//processRead(readSeq, dp, fraction,  hmmFlank, hmmPad, period,  outOS );
				readSequences.add(readSeq);
			}// while
			iter.close();
			os.close();
			//readSequences: an array of reads

			ProfileDP dpBatch = new ProfileDP(hmmSeq, hmmFlank + hmmPad, hmmFlank + hmmPad + str.getPeriod() - 1);//-1 for 0-index, inclusive
			processBatch(readSequences, dpBatch, fraction,  hmmFlank, hmmPad, period,  outOS, accumulator_batch );

			outOS.print(trVar.toString(headers));
			outOS.print('\n');
			Number[][] genotypes = RepeatCluster.genotype(accumulator);
			Number[][] genotypes_batch = RepeatCluster.genotype(accumulator_batch);
			Number[][] genotypes1 = clustering(accumulator);
			Number[][] genotypes_batch1 = clustering(accumulator_batch);
			Number[][] genotypes2 = clustering(RepeatCluster.removeOutliers(accumulator,2));
			Number[][] genotypes_batch2 = clustering(RepeatCluster.removeOutliers(accumulator_batch,2));
		
			outOS.print("Genotypes|Counts_outliers\t"+Arrays.asList(genotypes[0])+"\t"+Arrays.asList(genotypes[1])+"\n");
			outOS.print("Genotypes_batch|Counts_outliers\t"+Arrays.asList(genotypes_batch[0])+"\t"+Arrays.asList(genotypes_batch[1])+"\n");
			outOS.print("Genotypes|Counts_elbow\t"+Arrays.asList(genotypes1[0])+"\t"+Arrays.asList(genotypes1[1])+"\n");
			outOS.print("Genotypes_batch|Counts_elbow\t"+Arrays.asList(genotypes_batch1[0])+"\t"+Arrays.asList(genotypes_batch1[1])+"\n");
			outOS.print("Genotypes|Counts_elbow_outliers\t"+Arrays.asList(genotypes2[0])+"\t"+Arrays.asList(genotypes2[1])+"\n");
			outOS.print("Genotypes_batch|Counts_elbouw_outliers\t"+Arrays.asList(genotypes_batch2[0])+"\t"+Arrays.asList(genotypes_batch2[1])+"\n");
		}// for

		reader.close();
		outOS.close();
	}

	static private void processBatch(ArrayList<Sequence> readBatch, ProfileDP dpBatch, double fraction, int hmmFlank, int hmmPad, int period, SequenceOutputStream outOS, List<ReadAllele> read_alleles ) throws IOException{


		for (int round = 0; round < 5;round ++){
			double myCost = 0;
			int countIns = 0, countDel = 0, countMG = 0, countMB = 0;
			//Estimate error profiles
			for (Sequence readSeq:readBatch){
				EmissionState bestState = dpBatch.align(readSeq);
				//TODO: make a filter here: select only eligible alignment
				double alignScore = bestState.getScore();
				countIns += bestState.getCountIns();
				countDel += bestState.getCountDel();
				countMG += bestState.getCountMG();
				countMB += bestState.getCountMB();				
				//System.out.println("Score " + alignScore + " vs " + readSeq.length()*2 + " (" + alignScore/readSeq.length() +")");
				double bestIter = bestState.getIter() + fraction;
				myCost += alignScore;
			}//for readSeq
			double sum = 3.0 + countMG + countMB + countIns + countDel;
			double insP = (countIns + 1.0) /sum;
			double delP = (countDel + 1.0) /sum;
			double matP = (countMG + countMB + 1.0) /sum;
			double matchP = (countMG + 1.0) / (countMG + countMB + 2.0);
			double misMatchP = 1 - matchP;
			System.out.printf("Total: %3d %3d %3d %3d %8.4f %8.4f %8.4f %8.4f\n", countMG, countMB, countIns, countDel,  insP, delP,  misMatchP,myCost);
			dpBatch.setTransitionProbability(matP, insP, delP);
			dpBatch.setMatchProbability(matchP);			
		}//round

		//actually calling
		int index = 0;
		
		for (Sequence readSeq:readBatch){
			MarkovExpert expert = new MarkovExpert(1);
			double costM = 0;
			for (int x = 0; x< readSeq.length();x++){
				int base = readSeq.getBase(x);					
				costM -= JapsaMath.log2(expert.update(base));
			}	

			double backGround = costM / readSeq.length() - 0.1;
			boolean pass = true;
			outOS.print("Markov " + costM + "\t" + (costM / readSeq.length()) + "\n");


			EmissionState bestState = dpBatch.align(readSeq);
			double alignScore = bestState.getScore();
			//System.out.println("Score " + alignScore + " vs " + readSeq.length()*2 + " (" + alignScore/readSeq.length() +")");
			double bestIter = bestState.getIter() + fraction;
			read_alleles.add(new ReadAllele(readSeq.getName(), bestIter));
			//accumulator.add(bestIter);
			profilePositions.clear();
			seqPositions.clear();
			costGeneration.clear();
			byteArray.clear();

			EmissionState lastState = bestState;				
			bestState = bestState.bwdState;

			while (bestState != null){
				profilePositions.add(bestState.profilePos);	
				seqPositions.add(bestState.profilePos);
				costGeneration.add(lastState.score - bestState.score);

				if (bestState.seqPos == lastState.seqPos)
					byteArray.add((byte)Alphabet.DNA.GAP);
				else
					byteArray.add(readSeq.getBase(lastState.seqPos));						

				lastState = bestState;
				bestState = bestState.bwdState;
			}					

			double costL = 0, costR = 0, costCurrentRep = 0, costRep = 0;
			int stateL = 0, stateR = 0, stateCurrentRep = 0, stateRep = 0;
			int baseL = 0, baseR = 0, baseCurrentRep = 0, baseRep = 0;
			int bSeqL = 0, bSeqR = 0, bSeqCurrentRep = 0, bSeqRep = 0;

			int lastProfilePos = -1, lastSeqPos = -1;

			for (int x = profilePositions.size() - 1; x >=0; x--){
				outOS.print(Alphabet.DNA().int2char(byteArray.get(x)));

				int profilePos = profilePositions.get(x);
				int seqPos = seqPositions.get(x);

				if (profilePos < hmmFlank + hmmPad){
					stateL ++;
					costL += costGeneration.get(x);

					if (lastProfilePos != profilePos)
						baseL ++;

					if (lastSeqPos != seqPos)
						bSeqL ++;

				}else if(profilePos > hmmFlank + hmmPad + period){
					stateR ++;
					costR += costGeneration.get(x);	

					if (lastProfilePos != profilePos)
						baseR ++;

					if (lastSeqPos != seqPos)
						bSeqR ++;
				}else{
					stateCurrentRep ++;
					costCurrentRep += costGeneration.get(x);

					stateRep ++;
					costRep += costGeneration.get(x);

					if (lastProfilePos != profilePos){
						baseRep ++;
						baseCurrentRep ++;
					}

					if (lastSeqPos != seqPos){
						bSeqRep ++;
						bSeqCurrentRep ++;
					}

				}

				//end of a repeat cycle
				if (profilePos < lastProfilePos){
					outOS.print("<-----------------REP " + costCurrentRep  
							+ " " + stateCurrentRep
							+ " " + (stateCurrentRep == 0?"inf": "" + (costCurrentRep/stateCurrentRep))
							+ " " + baseCurrentRep
							+ " " + (baseCurrentRep == 0?"inf": "" + (costCurrentRep/baseCurrentRep))
							+ " " + bSeqCurrentRep
							+ " " + (bSeqCurrentRep == 0?"inf": "" + (costCurrentRep/bSeqCurrentRep))
							);

					if (costCurrentRep/bSeqCurrentRep > backGround){
						pass = false;
						outOS.print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
					}

					outOS.println();	 
					costCurrentRep = 0;//restart
					stateCurrentRep = 0;//restart
					baseCurrentRep = 0;
					bSeqCurrentRep = 0;
				}

				//left 
				if (profilePos >= hmmFlank + hmmPad && lastProfilePos < hmmFlank + hmmPad){
					outOS.print("<-----------------LEFT " + costL 
							+  " " + stateL 
							+  " " + (stateL == 0?"inf": "" + (costL/stateL))
							+  " " + baseL
							+  " " + (baseL == 0?"inf": "" + (costL/baseL))
							+  " " + bSeqL
							+  " " + (bSeqL == 0?"inf": "" + (costL/bSeqL))
							);
					if (costL/bSeqL > backGround){
						pass = false;
						outOS.print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
					}

					outOS.println();

				}

				//right
				//if (profilePos < hmmFlank + hmmPad + period && lastProfilePos >= hmmFlank + hmmPad + period){				
				//	outOS.print("<-----------------RIGHT " + costR
				//			+  " " + stateR
				//			+  " " + (stateR == 0?"inf": "" + (costR/stateR))
				//			+  " " + baseR
				//			+  " " + (baseR == 0?"inf": "" + (costR/baseR))
				//			+  " " + bSeqR
				//			+  " " + (bSeqR == 0?"inf": "" + (costR/bSeqR))						
				//			);
				//	outOS.println();
				//}	
				lastProfilePos = profilePos;
				lastSeqPos = seqPos;	

			}//for x
			//move to out of the loop
			outOS.print("<-----------------RIGHT " + costR
					+  " " +  stateR
					+  " " + (stateR == 0?"inf": "" + (costR/stateR))
					+  " " +  baseR
					+  " " + (baseR == 0?"inf": "" + (costR/baseR))
					+  " " +  bSeqR
					+  " " + (bSeqR == 0?"inf": "" + (costR/bSeqR))						
					);

			if (costR/bSeqR > backGround){
				pass = false;
				outOS.print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
			}
			//outOS.println();

			outOS.println();
			outOS.print ("L = " + (costL/(hmmFlank + hmmPad)) + " R = " + costR/(dpBatch.getProfileLength() - hmmFlank - hmmPad - period) + "\n");

			/*****************************************************************/				
			outOS.print("##" + readSeq.getName() +"\t"+bestIter+"\t"+readSeq.length() +"\t" +alignScore+"\t" + alignScore/readSeq.length() + '\t' + costM + "\t" + costM / readSeq.length() + "\t"  + costL + "\t" + stateL + "\t" + costR + "\t" + stateR + "\t" + (alignScore - costL - costR) + "\t" + stateRep + "\t" + pass + '\n');			
			outOS.print("==================================================================\n");	
			
		}

		//List<Cluster<ReadAllele>> clusters = clustering(read_alleles);
		//outOS.print(" Number of alleles: " + clusters.size() + ":\n");
		//for (Cluster<ReadAllele> cluster: clusters){
		//	outOS.print("  Allele " + cluster.getCenter().copy_number + " with " + cluster.getPoints().size() + " reads \n");
		//}

	}
	/*******************************************************************/				


	static private Number[][] clustering(List<ReadAllele> accumulator){
		Number[] res = new Number[]{Double.NaN, Double.NaN};
		if(accumulator.size()<2) return new Number[][] {res};
		KMeansPlusPlusClusterer<ReadAllele> clusterer = new KMeansPlusPlusClusterer<ReadAllele>(new Random());
		List<Cluster<ReadAllele>> clusters2 = clusterer.cluster(accumulator,2,20);
		List<Cluster<ReadAllele>> clusters1 = clusterer.cluster(accumulator,1,20);
		double  elbow2 = ReadAllele.inertia(clusters2);
		double  elbow1 = ReadAllele.inertia(clusters1);
		List<Cluster<ReadAllele>> clusters = elbow2 < elbow1 ? clusters2 : clusters1;
		Number[] n1 = new Number[clusters.size()];
		Number[] n2 = new Number[clusters.size()];
		for(int i=0; i<n1.length; i++){
			n1[i] = clusters.get(i).getCenter().copy_number;
			n2[i] = clusters.get(i).getPoints().size();
		}
		return new Number[][]{n1,n2};
	}


	static private void processRead(Sequence readSeq, ProfileDP dp, double fraction, int hmmFlank, int hmmPad, int period, SequenceOutputStream outOS, List<ReadAllele> read_alleles ) throws IOException{


		MarkovExpert expert = new MarkovExpert(1);
		double costM = 0;
		for (int x = 0; x< readSeq.length();x++){
			int base = readSeq.getBase(x);					
			costM -= JapsaMath.log2(expert.update(base));
		}	

		double backGround = costM / readSeq.length() - 0.1;
		boolean pass = true;


		outOS.print("Markov " + costM + "\t" + (costM / readSeq.length()) + "\n");

		EmissionState bestState = dp.align(readSeq);
		double alignScore = bestState.getScore();
		//System.out.println("Score " + alignScore + " vs " + readSeq.length()*2 + " (" + alignScore/readSeq.length() +")");
		double bestIter = bestState.getIter() + fraction;

		read_alleles.add(new ReadAllele(readSeq.getName(), bestIter));
		//accumulator.add(bestIter);
		
		/*******************************************************************/				
		profilePositions.clear();
		seqPositions.clear();
		costGeneration.clear();
		byteArray.clear();

		//double oldCost = bestState.score;
		EmissionState lastState = bestState;				
		bestState = bestState.bwdState;

		while (bestState != null){
			profilePositions.add(bestState.profilePos);	
			seqPositions.add(bestState.profilePos);
			costGeneration.add(lastState.score - bestState.score);

			if (bestState.seqPos == lastState.seqPos)
				byteArray.add((byte)Alphabet.DNA.GAP);
			else
				byteArray.add(readSeq.getBase(lastState.seqPos));						

			lastState = bestState;
			bestState = bestState.bwdState;
		}					

		double costL = 0, costR = 0, costCurrentRep = 0, costRep = 0;
		int stateL = 0, stateR = 0, stateCurrentRep = 0, stateRep = 0;
		int baseL = 0, baseR = 0, baseCurrentRep = 0, baseRep = 0;
		int bSeqL = 0, bSeqR = 0, bSeqCurrentRep = 0, bSeqRep = 0;

		int lastProfilePos = -1, lastSeqPos = -1;

		for (int x = profilePositions.size() - 1; x >=0; x--){
			outOS.print(Alphabet.DNA().int2char(byteArray.get(x)));

			int profilePos = profilePositions.get(x);
			int seqPos = seqPositions.get(x);

			if (profilePos < hmmFlank + hmmPad){
				stateL ++;
				costL += costGeneration.get(x);

				if (lastProfilePos != profilePos)
					baseL ++;

				if (lastSeqPos != seqPos)
					bSeqL ++;

			}else if(profilePos > hmmFlank + hmmPad + period){
				stateR ++;
				costR += costGeneration.get(x);	

				if (lastProfilePos != profilePos)
					baseR ++;

				if (lastSeqPos != seqPos)
					bSeqR ++;
			}else{
				stateCurrentRep ++;
				costCurrentRep += costGeneration.get(x);

				stateRep ++;
				costRep += costGeneration.get(x);

				if (lastProfilePos != profilePos){
					baseRep ++;
					baseCurrentRep ++;
				}

				if (lastSeqPos != seqPos){
					bSeqRep ++;
					bSeqCurrentRep ++;
				}

			}

			//end of a repeat cycle
			if (profilePos < lastProfilePos){
				outOS.print("<-----------------REP " + costCurrentRep  
						+ " " + stateCurrentRep
						+ " " + (stateCurrentRep == 0?"inf": "" + (costCurrentRep/stateCurrentRep))
						+ " " + baseCurrentRep
						+ " " + (baseCurrentRep == 0?"inf": "" + (costCurrentRep/baseCurrentRep))
						+ " " + bSeqCurrentRep
						+ " " + (bSeqCurrentRep == 0?"inf": "" + (costCurrentRep/bSeqCurrentRep))
						);

				if (costCurrentRep/bSeqCurrentRep > backGround){
					pass = false;
					outOS.print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
				}

				outOS.println();	 
				costCurrentRep = 0;//restart
				stateCurrentRep = 0;//restart
				baseCurrentRep = 0;
				bSeqCurrentRep = 0;
			}

			//left 
			if (profilePos >= hmmFlank + hmmPad && lastProfilePos < hmmFlank + hmmPad){
				outOS.print("<-----------------LEFT " + costL 
						+  " " + stateL 
						+  " " + (stateL == 0?"inf": "" + (costL/stateL))
						+  " " + baseL
						+  " " + (baseL == 0?"inf": "" + (costL/baseL))
						+  " " + bSeqL
						+  " " + (bSeqL == 0?"inf": "" + (costL/bSeqL))
						);
				if (costL/bSeqL > backGround){
					pass = false;
					outOS.print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
				}

				outOS.println();

			}

			//right
			//if (profilePos < hmmFlank + hmmPad + period && lastProfilePos >= hmmFlank + hmmPad + period){				
			//	outOS.print("<-----------------RIGHT " + costR
			//			+  " " + stateR
			//			+  " " + (stateR == 0?"inf": "" + (costR/stateR))
			//			+  " " + baseR
			//			+  " " + (baseR == 0?"inf": "" + (costR/baseR))
			//			+  " " + bSeqR
			//			+  " " + (bSeqR == 0?"inf": "" + (costR/bSeqR))						
			//			);
			//	outOS.println();
			//}	
			lastProfilePos = profilePos;
			lastSeqPos = seqPos;	

		}//for x

		//move to out of the loop
		outOS.print("<-----------------RIGHT " + costR
				+  " " +  stateR
				+  " " + (stateR == 0?"inf": "" + (costR/stateR))
				+  " " +  baseR
				+  " " + (baseR == 0?"inf": "" + (costR/baseR))
				+  " " +  bSeqR
				+  " " + (bSeqR == 0?"inf": "" + (costR/bSeqR))						
				);

		if (costR/bSeqR > backGround){
			pass = false;
			outOS.print(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
		}
		//outOS.println();

		outOS.println();
		outOS.print ("L = " + (costL/(hmmFlank + hmmPad)) + " R = " + costR/(dp.getProfileLength() - hmmFlank - hmmPad - period) + "\n");

		/*****************************************************************/				
		outOS.print("##" + readSeq.getName() +"\t"+bestIter+"\t"+readSeq.length() +"\t" +alignScore+"\t" + alignScore/readSeq.length() + '\t' + costM + "\t" + costM / readSeq.length() + "\t"  + costL + "\t" + stateL + "\t" + costR + "\t" + stateR + "\t" + (alignScore - costL - costR) + "\t" + stateRep + "\t" + pass + '\n');			
		outOS.print("==================================================================\n");	
		
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
	 * 
	 * @param seqList
	 * @param indexStart
	 *            : the start index of the list (inclusive)
	 * @param indexEnd
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

	static void aaa(Sequence seq, TandemRepeat str, int flanking, SamReader reader, int qual, String prefix,  int np) throws IOException, InterruptedException{
		String cmd = "kalign -gpo 60 -gpe 10 -tgpe 0 -bonus 0 -q -i " + prefix
				+ "i.fasta -o " + prefix + "o.fasta";

		//cmd = "clustalo --force -i " + prefix + "i.fasta -o " + prefix
		//			+ "o.fasta";


		String chrom = str.getChr();
		int start = str.getStart() - flanking;
		int end = str.getEnd() + flanking;

		if (start < 0) start = 0;
		if (end > seq.length())
			end = seq.length();

		SAMRecordIterator iter 
		= reader.query(chrom, start, end, false);


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
			LOG.info("Running " + cmd);
			Process process = Runtime.getRuntime().exec(cmd);
			process.waitFor();
			LOG.info("Done " + cmd);

			SequenceReader hmmSeqReader 
			= FastaReader.getReader(prefix	+ "i.fasta");
			Sequence readSeq;	

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
	}
}



