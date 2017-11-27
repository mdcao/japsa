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

package japsadev.tools;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.seq.XAFReader;
import japsa.util.BetaBinomialModel;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import org.apache.commons.math3.distribution.NormalDistribution;


/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.dev.captureVNTR", 
scriptDesc = "VNTR typing using capture sequencing")
public class CaptureVNTR extends CommandLine{	
	
static	 String[][] map = new String[][] {new String[] {"Illumina_NA12878" ,  "hg19_SMGcore.LCDG.20140908.006_S6_bmem.rdepth"}, 
			 new String[] {"Illumina_NA12878_1", "hg19_NA12878_800bp_Post_Capture_S2_bmem.rdepth2"},
		 new String[] {"Illumina_NA12891" ,  "hg19_NA12891_800bp_Post_Capture_S4_bmem.rdepth2"},
			 new String[] {"Illumina_NA12877" ,  "hg19_SMGcore.LCDG.20140908.001_S1_bmem.rdepth"  },
				 new String[] {"Illumina_NA12879" ,  "hg19_SMGcore.LCDG.20140908.002_S2_bmem.rdepth"  },
					 new String[] {"Illumina_NA12889" ,  "hg19_SMGcore.LCDG.20140908.003_S3_bmem.rdepth"  },
						 new String[] {"Illumina_NA12890" ,  "hg19_SMGcore.LCDG.20140908.004_S4_bmem.rdepth"  },
							 new String[] {"Illumina_NA12892"  , "hg19_SMGcore.LCDG.20140908.005_S5_bmem.rdepth"} } ;

	static Map<String, String>hm =  new HashMap();
	static{
		for(int i=0; i<map.length; i++){
			hm.put(map[i][1], map[i][0]);
		}
	}
	
	public CaptureVNTR(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("xafFile", "VNTR.xaf",  "Name of repeat file");
		addString("reference", "hg19.fas", "Name of reference genome");
		addString("target", "target.fa", "Where to write the target");		
		addString("technology", "pacbio", "Technology: pacbio or illumina");
		addString("directory", "./", "Directory with depth files");
		addInt("stat", 2, "0,1,2");
		addInt("readLength", 250, "Read length");

		addInt("stage", 6, "Stage of processing:\n"
			+ "0: Generate hmm profile, technology parameter is required\n"
			+ "1: Extract target sequences, prepare for alignment\n"
			+ "2: Look for reads spanning any of the repeats\n"
			+ "3: Read depth information extraction\n"
			+ "4:  ---\n"
			+ "5: Read depth analysis\n"
			+ "6: Read depth analysis with likelihood\n"
			);
		addString("output", "-", "Name of output file, - for standard out");
		addInt("pad", 10, "Gaps");
		addString("resample", null, "reference sample");
		addString("resAllele", null, "reference alleles");
		addString("CI", "95", "Confidence interval between 0 and 100");

		//addInt("CI", 95, "Confidence Interval");

		
		//addBoolean("reverse",false,"Reverse sort order");

		addStdHelp();		
	}  


   /*note the output file also includes information on downsampling */
	public static void main(String [] args) throws IOException, InterruptedException{
		CommandLine cmdLine = new CaptureVNTR();		
		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/	


		String xafFile       =  cmdLine.getStringVal("xafFile");	
		String referenceFile =  cmdLine.getStringVal("reference");
		String targetFile    =  cmdLine.getStringVal("target");
		String technology    =  cmdLine.getStringVal("technology");
		String output      =  cmdLine.getStringVal("output");
		String resample    =  cmdLine.getStringVal("resample");
		String resAllele    =  cmdLine.getStringVal("resAllele");
		int    stat   =  cmdLine.getIntVal("stat");
		int    readLength   =  cmdLine.getIntVal("readLength");

		int stage = cmdLine.getIntVal("stage");	
		String CI = cmdLine.getStringVal("CI");	
		//stage = 6;
		
		//int pad = cmdLine.getIntVal("pad");		

		if (stage == 0)
			generateProfile(xafFile, technology);

		if (stage == 1)
			stage1_extractTargetSequence(referenceFile, xafFile, targetFile);

		//Post stage 1:
		//bwa index target.fa
		//bwa mem
		//samtools view/sort/index/depth				
		if (stage == 2){
			for (int i= 0; i< args.length; i++)
				stage2_spanRead(args[i], xafFile);
		}
		if (stage == 3){
			System.out.println("This funciton no longer here, please use jsa.tr.trdepth");		
		}
		if (stage == 4){
			System.out.println("This funciton no longer here, please use stage 5");
		}

		if (stage == 5){
			stage5_readDepthAnalysis(xafFile,resample, args, output,stat, readLength);
		}

		if (stage == 6){
			String dir = cmdLine.getStringVal("directory");
			/*for(int i=0; i<args.length; i++){
				args[i] = dir+"/"+args[i];
			}*/
			String[] resamples = resample.split(":");
			/*for(int i=0; i<resamples.length; i++){
				resamples[i] = dir+"/"+resamples[i];
			}*/
			stage6_readDepthAnalysis(xafFile,resamples, resAllele, new File(dir+"/"+args[0]), output,stat, readLength, CI.split(":"));
		}

	}
	/**
	 * Analyse data resulted from stage 3
	 * @param xafFile: information of repeats
	 * @param resample: reference sample
	 * @param sFiles : analysied samples
	 * @param outputFile
	 * @throws IOException
	 * @throws InterruptedException
	 */
	static void stage5_readDepthAnalysis(String xafFile, String rData, String[] sFiles, String outputFile, int stat, int readLength) throws IOException, InterruptedException{
		if (sFiles.length ==0)
			return;	

		String [] sampleID = new String[sFiles.length];

		XAFReader [] sReaders = new XAFReader[sFiles.length];
		for (int i = 0; i < sFiles.length;i++){
			File file = new File(sFiles[i]);
			sampleID[i] = file.getName();
			sampleID[i] = sampleID[i].replaceAll(".xaf", "");
			sReaders[i] = new XAFReader(sFiles[i]);
		}

		XAFReader xafReader = new XAFReader(xafFile);
		XAFReader rReader = new XAFReader(rData);		

		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(outputFile);	
		//sos.print("#H:ID\tchrom\tstart\tend\trepLen\tseqLen");
		sos.print("#H:ID\tchrom\tstart\tend");
		for (int i = 0; i < sFiles.length;i++){
			sos.print("\t" + sampleID[i]);
		}		
		sos.print('\n');

		int repCol = stat * 2 + 6;
		int seqCol = repCol + 1;
		if (stat == 3){
			repCol = 14;
			seqCol = 12;			
		}
		if (stat == 4){
			repCol = 14;
			seqCol = 13;			
		}
		if (stat == 5){
			repCol = 14;
			seqCol = 15;			
		}
		if (stat == 6){
			repCol = 14;
			seqCol = 16;
		}

		while (xafReader.next() != null){
			//Extract			
			String ID = xafReader.getField("ID");
			String chrom = xafReader.getField("chrom");
			int startRep  = Integer.parseInt(xafReader.getField("start"));
			int endRep = Integer.parseInt(xafReader.getField("end"));

			//sum of lengths of the flanks
			int flankings = Integer.parseInt(xafReader.getField("lflank")) + Integer.parseInt(xafReader.getField("rflank")) + readLength;
			double period = Double.parseDouble(xafReader.getField("period"));
			double lengthR = period * Double.parseDouble(xafReader.getField("unitNo")) - readLength;

			rReader.next();
			if( !ID.equals(rReader.getField("ID"))){
				Logging.exit("Wrong ID at line " + rReader.lineNo(),1);
			}

			for (int i = 0; i < sFiles.length;i++){
				sReaders[i].next();
			}

			int samplingSize = 1000;

			//int recNo = xafReader.recordNo();

			//if (recNo == 8 || recNo == 32 || recNo == 38 || recNo == 57 || recNo == 64|| recNo == 65 || recNo == 86 || recNo == 87 || recNo == 93 || recNo == 100||recNo ==109 || recNo == 112 || recNo == 114 || recNo == 120){
			//}else
			//	continue;

			sos.print(ID);
			sos.print('\t');
			sos.print(chrom);
			sos.print('\t');
			sos.print(startRep);
			sos.print('\t');
			sos.print(endRep);
			//}				

			double countR = Double.parseDouble(rReader.getField(repCol));
			double totR   = Double.parseDouble(rReader.getField(seqCol));
			double ratioR = countR / totR;

			for (int i = 0; i < sFiles.length;i++){
				if( !ID.equals(sReaders[i].getField("ID"))){
					Logging.exit("Wrong ID at line " + rReader.lineNo(),1);
				}

				double countS = Double.parseDouble(sReaders[i].getField(repCol));
				double totS   = Double.parseDouble(sReaders[i].getField(seqCol));

				double ratioS = countS / totS;
				sos.print('\t');
				sos.print(ratioS/ratioR);

				NormalDistribution dist = BetaBinomialModel.ratioDistribution(totR - countR, totR, totS - countS, totS, samplingSize);

				double lengthShigh   = (flankings + lengthR) / (dist.getMean() - 2 * dist.getStandardDeviation()) - flankings; 
				double lengthSlow  = (flankings + lengthR) / (dist.getMean() + 2 * dist.getStandardDeviation()) - flankings;	
				double lengthS = (flankings + lengthR) / (dist.getMean()) - flankings;
				/*
				 * note:
				 * 
				 * r = (2f+l_s)/(2f+l_r)
				 * r_low = mean-2sd
				 * r_high=mean+2sd
				 * 
				 * lS =  (2f + lR)  * r - 2f
				 * 
				 */				

				sos.print('\t');
				sos.print((lengthR + readLength)/period);

				sos.print('\t');
				sos.print((lengthS + readLength)/period);

				sos.print('\t');				
				sos.print((lengthSlow + readLength)/period);
				sos.print('\t');
				sos.print((lengthShigh + readLength)/period);


				//sos.print('\t');
				//sos.print(dist.getMean() + "\t" + dist.getStandardDeviation());				
			}//for i
			sos.println();			
		}//while iter
		sos.close();
		xafReader.close();

		rReader.close();
		for (int i = 0; i < sFiles.length;i++){
			sReaders[i].close();
		}
	}

	 static double downsample = 1; // this will need to be modified mannually

	/**
	 * Analyse data resulted from stage 3
	 * @param xafFile: information of repeats
	 * @param resample: reference sample
	 * @param sDir:  file with samples to analyse
	 * @param outputFile
	 * @throws IOException
	 * @throws InterruptedException
	 */
	static void stage6_readDepthAnalysis(String xafFile, String[] rData, String resAllele, File sDir, String outputFile, int stat, int readLength1, String[] cistring) throws IOException, InterruptedException{
		File[] sFiles = sDir.listFiles();
		if (sFiles.length ==0)
			return;	
		double[] CI = new double[cistring.length]; //NOTE :  THIS IS A TERRIBLE WAY TO SET CI, SHOULD ADD AS COMMAND LINE PARAMETER
		for(int i=0; i<CI.length; i++){
			CI[i] = Double.parseDouble(cistring[i])/100.0;
			if(CI[i]<0 || CI[i]>1) throw new RuntimeException("CI needs to be between 0 and 100");
		}
		double readLength = (double) readLength1;
		String [] sampleID = new String[sFiles.length];

		XAFReader [] sReaders = new XAFReader[sFiles.length];
		for (int i = 0; i < sFiles.length;i++){
			File file = sFiles[i];
			sampleID[i] = file.getName();
			sampleID[i] = sampleID[i].replaceAll(".xaf", "");
			sReaders[i] = new XAFReader(file.getAbsolutePath());
		}

		XAFReader xafReader = new XAFReader(xafFile);
		XAFReader[] rReader = new XAFReader[rData.length];
		for(int i=0; i<rData.length; i++){
			rReader[i] = new XAFReader(sDir.getAbsolutePath()+"/"+rData[i]);	
		}	
		XAFReader rAlleleReader = new XAFReader(resAllele);

		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(outputFile);	
		VNTRGenotyper vg = new VNTRGenotyper(downsample);

		//sos.print("#H:ID\tchrom\tstart\tend\trepLen\tseqLen");
		sos.print("#H:ID\tchrom\tstart\tend");
		for (int i = 0; i < sFiles.length;i++){
			sos.print("\t" + sampleID[i]);
		}
		sos.print("\tsumCountS;sumTotS;countR;totR");
		sos.print('\n');

		int repCol = stat * 2 + 6;
		int seqCol = repCol + 1;
		if (stat == 3){
			repCol = 14;
			seqCol = 12;			
		}
		if (stat == 4){
			repCol = 14;
			seqCol = 13;			
		}
		if (stat == 5){
			repCol = 14;
			seqCol = 15;			
		}
		if (stat == 6){
			repCol = 14;
			seqCol = 16;			
		}		
		double sumCountS, sumTotS;
		
		for (; xafReader.next() != null;){
			
			//Extract			
			String ID = xafReader.getField("ID");
			
			String chrom = xafReader.getField("chrom");
			int startRep  = Integer.parseInt(xafReader.getField("start"));
			int endRep = Integer.parseInt(xafReader.getField("end"));

			//sum of lengths of the flanks
			//int    flankings = Integer.parseInt(xafReader.getField("lflank")) + Integer.parseInt(xafReader.getField("rflank")) + readLength;
			//double period = Double.parseDouble(xafReader.getField("period"));
			//double lengthR = period * Double.parseDouble(xafReader.getField("unitNo"));// - readLength;
			for(int i=0; i<rReader.length; i++){
				rReader[i].next();
			}
			if( !ID.equals(rReader[0].getField("ID"))){
				Logging.exit("Wrong ID at line a " + rReader[0].lineNo(),1);
			}
			
			rAlleleReader.next();
			if( !ID.equals(rAlleleReader.getField("ID"))){
				Logging.exit("Wrong ID at line b " + rAlleleReader.lineNo(),1);
			}
			
			for (int i = 0; i < sFiles.length;i++){
				sReaders[i].next();
			}		
			
			String refAlleleStr =  rAlleleReader.getField("refAllele");
			if (refAlleleStr == null){
				throw new RuntimeException("is null");
//				continue;//while
			}
			
			double refAllele = 0;
			try{
				refAllele = Double.parseDouble(refAlleleStr);
			}catch(Exception e){
//				continue;//while
				refAllele = Double.NaN;
			}			
//			if (refAllele == 0)
	//			
	//			continue;//while
				

			sos.print(ID);
			sos.print('\t');
			sos.print(chrom);
			sos.print('\t');
			sos.print(startRep);
			sos.print('\t');
			sos.print(endRep);
		//	sos.print('\t');
	//		sos.print(refAllele);
			//}				
			//repCol = 11; seqCol = 12;
			double countR = 0;
			
			double totR   =0;
			for(int i=0; i<rReader.length; i++){
				countR +=Double.parseDouble(rReader[i].getField(repCol));
				totR += Double.parseDouble(rReader[i].getField(seqCol));
			}
			//double ratioR = countR / totR;
			countR  = countR / readLength;
			totR = totR / readLength;
			
			
			
			sumCountS =0;
			sumTotS =0;
			
			for (int i = 0; i < sFiles.length;i++){
				
					//	String[] name = sFiles[i].split("/");
						String nme = sFiles[i].getName();
					//	System.err.println(nme);
				if( !ID.equals(sReaders[i].getField("ID"))){
					Logging.exit("Wrong ID at line " + rReader[0].lineNo(),1);
				}

				double countS = Double.parseDouble(sReaders[i].getField(repCol));
				double totS   = Double.parseDouble(sReaders[i].getField(seqCol));
				sumCountS+=countS;
				sumTotS+=totS;
				countS  = countS / readLength;
				totS = totS / readLength;
				
//				sos.print('\t');
				
			

				//double ref = 10;
				//double sample = refAllele;//unknown
				//
				//double est = (countS/(totS - countS))/(countR)
				vg.setRef(refAllele, countR, totR - countR) ;
				vg.setSample(countS, totS - countS);
				
				vg.bdw = 1;
			
				String result = vg.getConfs(CI);

								sos.print('\t');
								sos.print(result);//+"("+String.format("%5.2g", mass[2]).trim()+")");
									
			}//for i
			double doc = countR*readLength/(endRep-startRep);
			Double[] top = new Double[] {sumCountS, sumTotS,countR*readLength, totR*readLength};
			for(int i=0; i<top.length; i++) top[i] = top[i]/downsample;
			sos.print("\t");
			sos.print(String.format("%6.3e;%6.3e;%6.3e;%6.3e", top).replaceAll("\\s+", ""));
			sos.println();
		}//while iter
		sos.close();
		xafReader.close();
		for (int i = 0; i < rReader.length;i++){
			rReader[i].close();
		}
		for (int i = 0; i < sFiles.length;i++){
			sReaders[i].close();
		}
	}
	/**
	 * Stage 0
	 * @param xafFile
	 * @param technology
	 * @throws IOException
	 * @throws InterruptedException
	 */
	static void generateProfile(String xafFile, String technology) throws IOException, InterruptedException{
		Profile profile = illuminaProfile;
		if (technology.startsWith("pac"))
			profile = pacbioProfile;

		XAFReader xafReader = new XAFReader(xafFile);

		int nSeqs = 50;
		Random rd = new Random(1);		

		//int index = 0;
		while (xafReader.next() != null){				
			TandemRepeat str = TandemRepeat.read(xafReader);
			//TODO:
			if (str.getPeriod() <= 8){
				throw new RuntimeException("period too small");
				//continue;
			}

			String repUnit = xafReader.getField("repUnit");//chromosome
			String target = technology + "_" + str.getID();

			Sequence repSeq = 
				new Sequence(dna, repUnit,	target);

			SequenceOutputStream sos = SequenceOutputStream.makeOutputStream("f.fasta"); 
			byte [] nBytes = new byte[repSeq.length() * 2];			

			for (int i =0; i< nSeqs; i++){
				int seqIndex = 0;
				int j = 0;
				while (j < repSeq.length()){
					double toss = rd.nextDouble();
					if (toss < profile.subs){//subs
						int current = repSeq.getBase(j);
						if (current > 3)
							current = rd.nextInt(4);
						current += rd.nextInt(3);
						current  = (current + 1) % 4;
						nBytes[seqIndex] = (byte) current;
						seqIndex ++;
						j ++;
					}else if (toss < profile.subs + profile.del){
						//deletion
						j ++;
					}else if (toss < profile.subs + profile.del + profile.ins){
						//insertion
						nBytes[seqIndex] = (byte) rd.nextInt(4);
						seqIndex ++;					
					}else{//copy
						nBytes[seqIndex] = repSeq.getBase(j);
						seqIndex ++;
						j ++;					
					}
				}
				Sequence newSeq = new Sequence(dna, nBytes, seqIndex, "seq"+i);
				//TODO: 
				if(newSeq.length() > 6)
					newSeq.writeFasta(sos);
			}//for
			sos.close();

			Process process =  Runtime.getRuntime().exec("rm -f " + target + "o.fasta");
			process.waitFor();

			process =  Runtime.getRuntime().exec("cp f.fasta " + target + "o.fasta");
			process.waitFor();			

			String cmd = "kalign -q -i f.fasta -o " + target + "o.fasta";
			Logging.info("Running " + cmd);


			process = Runtime.getRuntime().exec(cmd);
			process.waitFor();
			Logging.info("Done " + cmd);

			process = Runtime.getRuntime().exec("hmmbuild --dna " + target +".hmm " + target + "o.fasta");
			process.waitFor();			
		}

		xafReader.close();		
	}

	static Alphabet dna = Alphabet.DNA();


	/**
	 * Read in the repeat information from xaf file and construct capture sequences 
	 * from the reference genome.
	 * @param referenceFile
	 * @param xafFile
	 * @param targetFile
	 * @throws IOException
	 */
	static void stage1_extractTargetSequence(String referenceFile, String xafFile, String targetFile) throws IOException{
		//Read in the genome
		Logging.info("Read genome begins");
		HashMap <String, Sequence> genome = new HashMap <String, Sequence>();
		SequenceReader reader = SequenceReader.getReader(referenceFile);
		Sequence seq;
		while ((seq = reader.nextSequence(dna)) != null){
			genome.put(seq.getName(), seq);
		}
		reader.close();
		Logging.info("Read genome done");

		XAFReader xafReader = new XAFReader(xafFile);
		SequenceOutputStream sos =  SequenceOutputStream.makeOutputStream(targetFile);

		System.out.println("#H:ID\tchrom\tstart\tend\tlflank\trflank");
		while (xafReader.next() != null){
			//Extract 
			String chrom = xafReader.getField("chrom");//chromosome
			int start  = Integer.parseInt(xafReader.getField("start"));
			int end    = Integer.parseInt(xafReader.getField("end"));			 		
			int rflank = Integer.parseInt(xafReader.getField("rflank"));
			int lflank = Integer.parseInt(xafReader.getField("lflank"));

			int mrflank = rflank + 300;
			int mlflank = lflank + 300;

			if (seq == null || !seq.getName().equals(chrom)){
				seq = genome.get(chrom);
			}
			if (seq == null){
				xafReader.close();
				sos.close();
				Logging.exit("Chrom in line " + xafReader.lineNo() + " not found!!!", 1);
			}
			Sequence s = seq.subSequence(start - mlflank - 1, end + mrflank);
			s.setName(chrom+"_"+start+"_"+end+"_"+mlflank);
			s.writeFasta(sos);

			System.out.println(xafReader.getField("ID") 
				+ "\t" + s.getName()
				+ "\t" + (mlflank + 1)
				+ "\t" + (mlflank + end -start)
				+ "\t" + lflank
				+ "\t" + rflank
				);

		}		

		xafReader.close();
		sos.close();
	}
	/**
	 * Stage 2
	 * @param bamFile
	 * @param pad
	 * @throws IOException 
	 */

	static void stage2_spanRead(String bamFile, String xafFile) throws IOException{		
//		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
//		SAMFileReader samReader = new  SAMFileReader(new File(bamFile));
		//SAMFileHeader samHeader = samReader.getFileHeader();
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFile));

		//String xafFile = "";
		XAFReader xafReader = new XAFReader(xafFile);


		while (xafReader.next() != null){
			String ID = xafReader.getField("ID");
			String chrom = xafReader.getField("chrom");
			int startRep  = Integer.parseInt(xafReader.getField("start"));
			int endRep = Integer.parseInt(xafReader.getField("end"));
			int period = Integer.parseInt(xafReader.getField("period"));

			SAMRecordIterator iter = samReader.query(chrom, startRep, endRep, false);
			int count = 0, countProper = 0;
			while (iter.hasNext()){
				SAMRecord sam = iter.next();
				if (sam.getAlignmentStart() > startRep)
					break;//already got inside the repeat
				if (sam.getAlignmentEnd() > endRep){
					count ++;
					if (sam.getReadPairedFlag()){
						if (!sam.getMateUnmappedFlag() && chrom.equals(sam.getMateReferenceName()))
							countProper ++;
					}

					Cigar cigar = sam.getCigar();
					int readPos = 0;
					int refPos = sam.getAlignmentStart();
					int indel = 0;

					for (final CigarElement e : cigar.getCigarElements()) {
						final int  length = e.getLength();
						switch (e.getOperator()) {
						case H :								                	
							break; // ignore hard clips
						case P : //pad is a kind of clipped							
							break; // ignore pads	                
						case S ://advance on the reference							
							readPos += length;
							break; // soft clip read bases	                	
						case N : 
							refPos += length; 
							break;  // reference skip

						case D :      	
							if (refPos >= startRep && refPos <= endRep){
								indel -= length;//need to fix this
							}
							refPos += length;
							break;	

						case I :
							if (refPos >= startRep && refPos <= endRep){
								indel += length;//need to fix this
							}
							readPos += length;
							break;
						case M :							
							readPos += length;
							refPos  += length;
							break;
						case EQ :
							readPos += length;
							refPos  += length;
							break;

						case X :
							readPos += length;
							refPos  += length;
							break;
						default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
						}//casse
					}//for
					System.out.println(ID + "\t" + (indel * 1.0/period)+ "\t"+sam.getCigarString()+"#"+(1 + sam.getAlignmentEnd() - sam.getAlignmentStart())+"#"+indel+"#"+period);					
				}
			}//while
			iter.close();
			if (count > 0)
				System.out.println(ID + "\t" + count+"\t"+countProper+"=====================================");

		}

		/*******************************************************************
		int nSeq = samHeader.getSequenceDictionary().size();
		for (int seqIndex = 0; seqIndex < nSeq; seqIndex++){
			String seqName = samHeader.getSequence(seqIndex).getSequenceName();

			//TODO: make this more genetic
			String [] toks = seqName.split("_");
			int start = Integer.parseInt(toks[3]);
			int end = Integer.parseInt(toks[2]) - Integer.parseInt(toks[1]) + start;			

			start -= pad;
			end += pad;		

			SAMRecordIterator iter = samReader.query(seqName, start, end, false);
			int count = 0, countProper = 0;
			while (iter.hasNext()){
				SAMRecord sam = iter.next();
				if (sam.getAlignmentStart() > start)
					break;//already got inside the repeat
				if (sam.getAlignmentEnd() > end){
					count ++;
					if (!sam.getMateUnmappedFlag() && seqName.equals(sam.getMateReferenceName()))
						countProper ++;
					System.out.println("#"+sam.getCigarString()+"#"+(1 + sam.getAlignmentEnd() - sam.getAlignmentStart()));


				}
			}//while
			iter.close();			
			System.out.println(seqName + "\t" + (end-start-pad-pad) + "\t"+ count+"\t"+countProper);
		}
		/*******************************************************************/		
		samReader.close();
		xafReader.close();
	}




	static class Profile{
		double subs = 0.1;
		double del =  0.05;
		double ins =  0.05;		

		Profile (double s, double d, double i){
			subs = s; del = d; ins = i;
		}	
	}

	//Subs = 1%, del = 0.01%, ins = 0.01%
	static Profile illuminaProfile = new Profile(0.01, 0.001, 0.001);	
	//Subs = 1%, del = 10%, ins = 5%
	static Profile pacbioProfile =   new Profile(0.01, 0.075, 0.075);


}
