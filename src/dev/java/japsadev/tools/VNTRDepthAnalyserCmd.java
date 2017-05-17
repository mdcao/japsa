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

import japsa.seq.SequenceOutputStream;
import japsa.seq.XAFReader;
import japsa.util.BetaBinomialModel;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.dev.vntrDepthAnalyser", 
scriptDesc = "VNTR typing using coverage depth information of Illumina capture sequencing")
public class VNTRDepthAnalyserCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(VNTRDepthAnalyserCmd.class);

	public VNTRDepthAnalyserCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("xafFile", "VNTR.xaf",  "Name of repeat file");		
		addInt("readLength", 250, "Average read length");
		addString("output", "-", "Name of output file, - for standard out");		
		addString("resample", null, "reference sample");
		addString("resAllele", null, "reference alleles");
		addString("model", "bn", "bn (beta-binomial) or nm (normal)");
		//addBoolean("reverse",false,"Reverse sort order");
		
		addStdHelp();		
	}  



	public static void main(String [] args) throws IOException, InterruptedException{
		CommandLine cmdLine = new VNTRDepthAnalyserCmd();		
		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/	


		String xafFile       =  cmdLine.getStringVal("xafFile");	
		String output      =  cmdLine.getStringVal("output");
		String resample    =  cmdLine.getStringVal("resample");
		String resAllele    =  cmdLine.getStringVal("resAllele");
		String model    =  cmdLine.getStringVal("model");
		int    readLength   =  cmdLine.getIntVal("readLength");

		depthAnalysis(xafFile,resample, resAllele, args, output, readLength,  model);
	}

	static void depthAnalysis(String xafFile, String rData, String resAllele, String[] sFiles, String outputFile, int readLength, String model) throws IOException, InterruptedException{
		if (sFiles.length ==0)
			return;	

		int stat = 2;
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
		XAFReader rAlleleReader = new XAFReader(resAllele);

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
			int    flankings = Integer.parseInt(xafReader.getField("lflank")) + Integer.parseInt(xafReader.getField("rflank"));
			double period = Double.parseDouble(xafReader.getField("period"));
			//double lengthR = period * Double.parseDouble(xafReader.getField("unitNo"));// - readLength;

			rReader.next();
			if( !ID.equals(rReader.getField("ID"))){
				LOG.error("Wrong ID at line a " + rReader.lineNo());
				System.exit(1);
			}

			rAlleleReader.next();
			if( !ID.equals(rAlleleReader.getField("ID"))){
				LOG.error("Wrong ID at line b " + rAlleleReader.lineNo());
				System.exit(1);
			}

			for (int i = 0; i < sFiles.length;i++){
				sReaders[i].next();
			}		

			String refAlleleStr =  rAlleleReader.getField("refAllele");
			if (refAlleleStr == null){
				continue;//while
			}

			double refAllele = 0;
			try{
				refAllele = Double.parseDouble(refAlleleStr);
			}catch(Exception e){
				continue;//while
			}			
			if (refAllele == 0)
				continue;//while

			double lengthR = period * refAllele;//


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
			//double ratioR = countR / totR;
			countR  = countR / readLength;
			totR = totR / readLength;

			for (int i = 0; i < sFiles.length;i++){
				if( !ID.equals(sReaders[i].getField("ID"))){
					LOG.error("Wrong ID at line " + rReader.lineNo());
					System.exit(1);
				}

				double countS = Double.parseDouble(sReaders[i].getField(repCol));
				double totS   = Double.parseDouble(sReaders[i].getField(seqCol));

				countS  = countS / readLength;
				totS = totS / readLength;

				sos.print('\t');

				/////////////////////////////////////////////////////////////////
				//Actual analysis here
				if (model.equals("nm")){
					NormalDistribution dist = BetaBinomialModel.ratioDistribution(totR - countR, totR, totS - countS, totS, 1000);
					double lengthShigh   = (flankings + lengthR) / (dist.getMean() - 2 * dist.getStandardDeviation()) - flankings; 
					double lengthSlow  = 1;//(flankings + lengthR) / (dist.getMean() + 2 * dist.getStandardDeviation()) - flankings;
				}else{
					double conf = 0.95;
					VNTRGenotyper vg = new VNTRGenotyper();

					vg.setRef(refAllele, countR, totR - countR) ;
					vg.setSample(countS, totS - countS);
					vg.bdw = 1;
					//				String result =  vg.getConf1(0.2);			
					String result = vg.getConf(0.95);
					sos.print(result);


					//			double[] reports = bn(refAllele, totR, countR, totS,  countS, conf);
					//					sos.print("("+reports[0] + "," + reports[1] + "," + reports[2] + "," + reports[3] + ")");






				}


				//System.err.println(genos[range[0]] +" to "+genos[range[1]] +" "+mass[2]+"\n max:"+mass[0]+ " maxp:"+mass[1]);
				//for(int x=range[0]; x<=range[1]; x++){
				//	System.err.println(genos[i]+" => "+prob[x]);
				//}

				//				sos.print('\t');
				//				sos.print((lengthR)/period);

				//				sos.print('\t');
				//				sos.print((lengthS + readLength)/period);

				//				sos.print('\t');				
				//				sos.print((lengthSlow + readLength)/period);
				//				sos.print('\t');
				//				sos.print((lengthShigh + readLength)/period);


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

	static double [] bn(double refAllele, double totR, double countR, double totS, double countS, double conf){
		/********************************************************************
		VNTRGenotyper vg = new VNTRGenotyper();

		vg.setRef(refAllele, countR, totR - countR) ;
		vg.setSample(countS, totS - countS);

		vg.bdw = 1;
		//				String result =  vg.getConf1(0.2);			
		String result = vg.getConf(0.95);



		double[] genos = new double[20];
		int half = genos.length/2;

		double startGeno = Math.max(0, (refAllele - half));  
		for (int x = 0; x < genos.length; x++){
			genos[x] = startGeno +x;
		}

		double[] prob = new double[genos.length];
		vg.probability(prob, genos);
		//System.err.println("Simulated  "+sample+" depth:"+depth);
		int[] range = new int[2];

		//double conf = 0.5;
		double[] mass = VNTRGenotyper.getconf(prob, genos, conf, range);

		//geno[maxi],  prob[maxi], sum
		//low, high, support,max likelihood genotype, prob of the max genotype
		return new double[] {genos[range[0]], mass[0], genos[range[1]], mass[2], mass[0], mass[1]};
		/********************************************************************/
		return null;
	}


	static double [] nm(double totR, double countR, double totS, double countS, int samplingSize){
		NormalDistribution dist = BetaBinomialModel.ratioDistribution(totR - countR, totR, totS - countS, totS, samplingSize);

		double lengthShigh   = 0;//(flankings + lengthR) / (dist.getMean() - 2 * dist.getStandardDeviation()) - flankings; 
		double lengthSlow  = 1;//(flankings + lengthR) / (dist.getMean() + 2 * dist.getStandardDeviation()) - flankings;
		return new double [] {lengthShigh, lengthSlow};
		//		double lengthS = (flankings + lengthR) / (dist.getMean()) - flankings;
		//		/*
		//		 * note:
		//		 * 
		//		 * r = (2f+l_s)/(2f+l_r)
		//		 * r_low = mean-2sd
		//		 * r_high=mean+2sd
		//		 * 
		//		 * lS =  (2f + lR)  * r - 2f
		//		 * 
		//		 */				
		//
		//		sos.print('\t');
		//		sos.print((lengthR + readLength)/period);/

		//		sos.print('\t');
		//		sos.print((lengthS + readLength)/period);

		//		sos.print('\t');				
		//		sos.print((lengthSlow + readLength)/period);
		//		sos.print('\t');
		//		sos.print((lengthShigh + readLength)/period);

	}

}
