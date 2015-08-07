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
 * 08/04/2012 - Minh Duc Cao: Revised                                        
 * 16/11/2013 - Minh Duc Cao 
 ****************************************************************************/

package japsa.bio.tr;

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.OutputStream;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;



import org.apache.commons.math3.distribution.NormalDistribution;
/**
 * FIXME: Need testing
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.trv.fragment2var",
            scriptDesc = "Analyse tandem repeat variation from fragment sizes")
public class Fragment2TRV {
	public static void main(String[] args) throws Exception {
		/*********************** Setting up script ****************************/		 
		String scriptName = "jsa.trv.fragment2var";
		String desc = "Analyse tandem repeat variation from fragment sizes\n";		
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [options] fragmentFile1 fragmentFile2 ...");
		/**********************************************************************/		

		cmdLine.addString("trFile", null,
				"Name of the tandem repeat position file",true);
		cmdLine.addString("output", "-",
				"Name of the output file (- for standard output)");

		cmdLine.addString("ansFile", null,	"Name of the answer file, ignored if null");		
		cmdLine.addString("bcfFile", null,	"Name of the BCF file, ignored if null");

		cmdLine.addInt("fm", 0,"Mean of insert size, enter <= 0 for estimating from data");
		cmdLine.addDouble("fd", 0,"Standard deviation of insert size, enter <= 0 for estimating from data");
		cmdLine.addInt("gap", 3,"Gap");

		cmdLine.addBoolean("sorted", true,
				"The insert file has been sorted");

		cmdLine.addStdHelp();	
		/**********************************************************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")){
			System.out.println(desc + cmdLine.usageMessage());			
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usageMessage());
			System.exit(-1);
		}	
		/**********************************************************************/

		String trFile = cmdLine.getStringVal("trFile");
		String outFile = cmdLine.getStringVal("output");
		String ansFile =  cmdLine.getStringVal("ansFile");
		String bcfFile =  cmdLine.getStringVal("bcfFile");

		int fm = cmdLine.getIntVal("fm");
		double fd = cmdLine.getDoubleVal("fd");		
		int gap = cmdLine.getIntVal("gap");



		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(outFile);
		//BufferedOutputStream out;
		//if (outFile.equals("-")){
		//	out = new  BufferedOutputStream(System.out);			
		//}else
		//	out = new BufferedOutputStream(new FileOutputStream(outFile));

		ArrayList<TandemRepeatVariant> ansList = null;
		if (ansFile != null)
			ansList = TandemRepeatVariant.readFromFile(ansFile);

		ArrayList<TandemRepeatVariant> trvList = TandemRepeatVariant.readFromFile(trFile);

		ArrayList<Indel> indelList = null;
		if (bcfFile != null){
			indelList = filterBCF(trvList, bcfFile);	
		}	

		//runAnalysisAAA(trvList,args, fm, fd, gap,out, ansList);

		if (cmdLine.getBooleanVal("sorted"))
			runAnalysisSort(trvList,args, fm, fd, gap,out, indelList, ansList);//new
		//runAnalysisAAA(trvList,args, fm, fd, gap,out, ansList);
		else
			runAnalysis(trvList,args, fm, fd, gap,out, ansList);//new

		out.close();

		// getInsertSizeSlow(samFile, output);		
		//System.out.println("##There are " + args.length + " libraries");
		//runAnalysis(anno,val,args, fm, fd, gap);//old
	}	
	
	/**
	 * Assume both inserts and annotations are sorted by the start position
	 * 
	 * @param iFile
	 * @param aFile
	 * @throws IOException
	 */

	static int bsearch(JapsaAnnotation anno, int start, int end, int gap){
		int low = 0, high = anno.numFeatures();


		while (low < high){
			int mid = (low + high) / 2;			
			//assert 0 <= mid < anno.size()
			JapsaFeature feature = anno.getFeature(mid);
			if (start >= feature.getStart() - gap){
				low = mid +1;
			}else if (end < feature.getEnd() + gap){
				high = mid;
			}else{
				return mid;
			}	

		}		
		return -1;
	}

	static int blsearch(ArrayList<TandemRepeatVariant> trvList, String chr, int start, int end, int gap){
		int low = 0, high = trvList.size();

		while (low < high){
			int mid = (low + high) / 2;			
			//assert 0 <= mid < anno.size()
			TandemRepeatVariant trv = trvList.get(mid);

			if (chr.compareTo(trv.getChr()) > 0)
				low = mid + 1;
			else if (chr.compareTo(trv.getChr()) < 0)
				high = mid;
			else if (start >= trv.getStart() - gap){
				low = mid +1;
			}else if (end <= trv.getEnd() + gap){
				high = mid;
			}else{
				return mid;
			}	

		}		
		return -1;
	}

	static void runAnalysisAAA(ArrayList<TandemRepeatVariant>  trvList, String[] insertFiles, double fm, double fd, int gap, OutputStream out, ArrayList<TandemRepeatVariant>  ansList) throws IOException {	
		//Perform an initial comparision (only if answer is ready)
		if (ansList != null)
			CompareTRV.compareStr(ansList, trvList,0,10000);

		double [][] checkVar = new double[insertFiles.length + 1][trvList.size()];

		for (int xx=0; xx < trvList.size(); xx++)
			checkVar[0][xx] =  trvList.get(xx).getVar();		

		for (int fIdx = 0; fIdx < insertFiles.length; fIdx ++){
			int [] counts = new int[trvList.size()];
			double [] sums = new double[trvList.size()];
			double [] sqSums = new double[trvList.size()];	

			//Estimated mean and variance
			//double [] [] eM = new double[numChr] [];
			//double [] [] eV = new double[numChr] [];
			BufferedReader in = SequenceReader.openFile(insertFiles[fIdx]);
			System.out.println("##Processing file " + insertFiles[fIdx]);

			//total (for global stats)
			int tCount = 0;
			double tSum = 0, tSqSum = 0;

			String line = "";
			//int noPair = 0;
			//Read the insert file to get all stats
			while ((line = in.readLine()) != null){
				if (line.startsWith("#"))
					continue;

				//noPair ++;
				String [] toks = line.split("\\t");
				int start = Integer.parseInt(toks[3]);
				int end = Integer.parseInt(toks[4]);
				tCount ++;
				double v = (end + 1 - start); 
				tSum += v;
				tSqSum += v * v;

				int mid = blsearch(trvList, toks[0], start, end, gap);
				if (mid >=0){				
					counts[mid] ++;
					sums[mid] += v;
					sqSums[mid] += v * v;
					//TandemRepeat tr = trList.get(mid);
					//System.out.println("In " + tr.getChr() + " : " + start + "  " + end + " include " + tr.getStart() + " - " + (tr.getEnd() - 1) + "  " + v + " -> " + sums[mid]);
				}
			}//while
			double gMean = tSum / tCount;//global mean
			double gVar = tSqSum / tCount  - gMean * gMean;//global variance

			if (fd > 0) gMean = fm;
			if (fd > 0) gVar = fd * fd;

			System.out.println("##Lib mean = " + gMean + "  std = " + Math.sqrt(gVar));

			//now doing inference
			for (int i = 0; i < trvList.size(); i++){
				TandemRepeatVariant tr = trvList.get(i);
				double mean = 0;//, var = 0;

				if (counts[i] != 0){
					mean = sums[i] / counts[i]; //sample mean
					//var = sqSums[i]  / counts[i] - mean * mean;//sample variance
					double meanP = tr.getMean(), varP = tr.getStd() * tr.getStd();

					double eVar = gVar * varP / (counts[i] * tr.getPeriod() * tr.getPeriod() * varP + gVar);
					double eMean = ((meanP / varP)  -   (counts[i] * tr.getPeriod() * (mean - gMean) / gVar)) * eVar;

					tr.setStd(Math.sqrt(eVar));
					tr.setMean(eMean);
					double predicted = eMean;
					//int predicted =  (int) Math.round(eMean);				
					tr.setVar(predicted);

					checkVar[fIdx + 1][i] = predicted;;

				}			
				NormalDistribution d = new NormalDistribution(tr.getMean(), tr.getStd());
				tr.setConfidence(d.cumulativeProbability(tr.getVar() - 0.5, tr.getVar() + 0.5));
			}//for i
			if (ansList != null)
				CompareTRV.compareStr(ansList, trvList,0,10000);
		}//for iFdx

		String[] headers = {TandemRepeat.chrHd, TandemRepeat.startHd, TandemRepeat.endHd, 
				TandemRepeat.periodHd, TandemRepeat.unitNoHd, TandemRepeatVariant.varHd, TandemRepeatVariant.confidenceHd, TandemRepeatVariant.meanHd, TandemRepeatVariant.stdHd};

		out.write(("#H:" + headers[0]).getBytes());
		for (int i=1; i < headers.length; i++)
			out.write(("\t"+headers[i]).getBytes());
		out.write('\n');
		for (int i = 0; i < trvList.size(); i++){
			out.write((trvList.get(i).toString(headers)).getBytes());
			if (ansList != null){
				out.write(("\t"+ansList.get(i).getVar()).getBytes());
				for (int xx = 0; xx < checkVar.length; xx++)
					out.write(("\t"+checkVar[xx][i]).getBytes());
			}
			//out.print();

			//	if (counts[i] > 0){
			//		double mean = (sums[i]/counts[i]);
			//		double var  = sqSums[i]/counts[i] - mean * mean;
			//		out.print("\t#"+counts[i]+"\t"+mean +"\t"+Math.sqrt(var)+"\t"+ var);
			//	}
			out.write('\n');
		}
		//ShortTandemRepeat.write(trvList, System.out, headers);
	}

	/**
	 * This method is obsoleted because of the better runAnalysisSort method 
	 * 
	 * @param trfFile
	 * @param insertFiles
	 * @param fm
	 * @param fd
	 * @param gap
	 * @param out
	 * @param ansList
	 * @throws IOException
	 */
	static void runAnalysis(ArrayList<TandemRepeatVariant>  trvList, String[] insertFiles, double fm, double fd, int gap, SequenceOutputStream out, ArrayList<TandemRepeatVariant>  ansList) throws IOException {		

		//Perform an initial comparision (only if answer is ready)
		if (ansList != null)
			CompareTRV.compareStr(ansList, trvList,0,10000);

		double [][] checkVar = new double[insertFiles.length + 1][trvList.size()];

		for (int xx=0; xx < trvList.size(); xx++)
			checkVar[0][xx] =  trvList.get(xx).getVar();		

		for (int fIdx = 0; fIdx < insertFiles.length; fIdx ++){
			int [] counts = new int[trvList.size()];
			double [] sums = new double[trvList.size()];
			double [] sqSums = new double[trvList.size()];	

			//Estimated mean and variance
			//double [] [] eM = new double[numChr] [];
			//double [] [] eV = new double[numChr] [];
			BufferedReader in = SequenceReader.openFile(insertFiles[fIdx]);	
			System.out.println("##Processing file " + insertFiles[fIdx]);

			//total (for global stats)
			int tCount = 0;
			double tSum = 0, tSqSum = 0;

			String line = "";
			//int noPair = 0;
			//Read the insert file to get all stats
			while ((line = in.readLine()) != null){
				if (line.startsWith("#"))
					continue;

				//noPair ++;
				String [] toks = line.split("\\t");
				int start = Integer.parseInt(toks[3]);
				int end = Integer.parseInt(toks[4]);
				tCount ++;

				double v = end + 1 - start;//Integer.parseInt(toks[5]);//(end + 1 - start); 
				tSum += v;
				tSqSum += v * v;

				int mid = blsearch(trvList, toks[0], start, end, gap);
				if (mid >=0){				
					counts[mid] ++;
					sums[mid] += v;
					sqSums[mid] += v * v;					
				}
			}//while
			double gMean = tSum / tCount;//global mean
			double gVar = tSqSum / tCount  - gMean * gMean;//global variance

			if (fd > 0) gMean = fm;
			if (fd > 0) gVar = fd * fd;

			System.out.println("##Lib mean = " + gMean + "  std = " + Math.sqrt(gVar) + " tCount = " + tCount);

			//now doing inference
			for (int i = 0; i < trvList.size(); i++){
				TandemRepeatVariant trv = trvList.get(i);
				double mean = 0;//, var = 0;
				System.out.println("#=#" + i + "\t" + counts[i]);
				if (counts[i] != 0){
					mean = sums[i] / counts[i]; //sample mean
					//var = sqSums[i]  / counts[i] - mean * mean;//sample variance
					double meanP = trv.getMean(), varP = trv.getStd() * trv.getStd();

					double eVar = gVar * varP / (counts[i] * trv.getPeriod() * trv.getPeriod() * varP + gVar);
					double eMean = ((meanP / varP) - (counts[i] * trv.getPeriod() * (mean - gMean) / gVar)) * eVar;
					//double eMean = ((meanP / varP) - (counts[i] * (mean - gMean) / gVar)) * eVar;

					trv.setStd(Math.sqrt(eVar));
					trv.setMean(eMean);
					double predicted = eMean;
					//int predicted =  (int) Math.round(eMean);				
					trv.setVar(predicted);

					checkVar[fIdx + 1][i] = predicted;

				}			
				NormalDistribution d = new NormalDistribution(trv.getMean(), trv.getStd());
				trv.setConfidence(d.cumulativeProbability(trv.getVar() - 0.5, trv.getVar() + 0.5));
			}//for i
			if (ansList != null)
				CompareTRV.compareStr(ansList, trvList,0,10000);
		}//for iFdx

		String[] headers = {TandemRepeat.chrHd, TandemRepeat.startHd, TandemRepeat.endHd, 
				TandemRepeat.periodHd, TandemRepeat.unitNoHd, TandemRepeatVariant.varHd, TandemRepeatVariant.confidenceHd, TandemRepeatVariant.meanHd, TandemRepeatVariant.stdHd};

		out.print("#H:" + headers[0]);
		for (int i=1; i < headers.length; i++)
			out.print("\t"+headers[i]);
		out.print('\n');
		for (int i = 0; i < trvList.size(); i++){
			out.print(trvList.get(i).toString(headers));
			//if (ansList != null){
			//	out.print("\t"+ansList.get(i).getVar());
			//	for (int xx = 0; xx < checkVar.length; xx++)
			//		out.print("\t"+checkVar[xx][i]);
			//}
			out.print('\n');
		}
	}

	/**
	 * Assume every thing is sorted, and the order of reference sequences
	 * in the fragment list (which is from the bam/sam file) is the same with
	 * in indel/tr list. Every sequence in the tr list should be included in 
	 * list of fragment list.
	 * 
	 * @param trfFile
	 * @param insertFiles
	 * @param fm
	 * @param fd
	 * @param gap
	 * @param out
	 * @param ansList
	 * @throws IOException
	 */
	static void runAnalysisSort(ArrayList<TandemRepeatVariant>  trvList, String[] insertFiles, double fm, double fd, int gap, SequenceOutputStream out, ArrayList<Indel> indelList, ArrayList<TandemRepeatVariant>  ansList) throws IOException {

		int indelIdx = 0;
		Indel indel = null;
		if (indelList != null && indelList.size() > 0){
			indel = indelList.get(0);
		}		
		//indel != null iff indelList != null

		//Perform an initial comparision (only if answer is ready)
		if (ansList != null)
			CompareTRV.compareStr(ansList, trvList,0,10000);

		for (int fIdx = 0; fIdx < insertFiles.length; fIdx ++){
			int [] counts = new int[trvList.size()];
			int countEff = 0;//count the number of effective fragments
			double [] sums = new double[trvList.size()];
			double [] sqSums = new double[trvList.size()];	

			int trvIndex = 0;//every tr prior to this will not overlap with the insert

			BufferedReader in = SequenceReader.openFile(insertFiles[fIdx]);	
			System.out.println("##Processing file " + insertFiles[fIdx]);

			//total (for global stats)
			int tCount = 0;
			double tSum = 0, tSqSum = 0;

			String line = "";
			//int noPair = 0;
			//Read the insert file to get all stats

			ArrayList<String> refIDs = new ArrayList<String>(10);

			int refIndexIndel = 0;
			int refIndexTR = 0;

			while ((line = in.readLine()) != null){
				if (line.startsWith("#"))
					continue;

				if (line.startsWith("@")){
					//parse to get the reference IDs
					if (line.startsWith("@SQ")){
						String [] toks = line.split("\t");
						for (int x = 0; x < toks.length;x++){
							if (toks[x].startsWith("SN:"))
								refIDs.add(toks[x].substring(3));									
						}
					} 

					continue;

				}


				if (line.length() < 10)
					continue;
				//if (tCount >= 7300)
				//	System.out.println("here");

				PEFragment fragment = PEFragment.readLine(line);
				//noPair ++;
				//String [] toks = line.split("\\t");
				int start = fragment.getStart();
				int end = fragment.getEnd();				
				////////////////////////////////////////////////////////////////////
				//////////////////////// Calibrate insert size /////////////////////

				//while (indel != null && 
				//		( (indel.chr.compareTo(fragment.getSeqID()) < 0) 
				//				|| (indel.chr.compareTo(fragment.getSeqID()) == 0 && indel.start < fragment.getStart()))){
				while (indel != null){ 
					//get the refIndex of the indel
					while (refIndexIndel < refIDs.size() && !indel.chr.equals(refIDs.get(refIndexIndel))){
						refIndexIndel ++;
					}
					//assert refIndex >= refID.size || indel.chr.equals(refIDs.get(refIndex) based on the
					if (refIndexIndel >= refIDs.size()){
						throw new RuntimeException("Sequence " + indel.chr + " in indel list not present in the reference genome");
					}

					if (refIndexIndel > fragment.getReferenceIndex())
						break;

					if (refIndexIndel == fragment.getReferenceIndex() 
							&& indel.start >= start)
						break;

					indelIdx ++;
					if (indelIdx > indelList.size()) 
						indel = null;
					else 
						indel = indelList.get(indelIdx);					
				}//while

				//the current indel is after the start of the fragment				
				if (indel !=null && refIndexIndel == fragment.getReferenceIndex()){
					Indel eIndel = indel;
					int eIdx = indelIdx;
					while (eIndel.chr.equals(indel.chr) && eIndel.start < fragment.getEnd()){
						fragment.iSize += eIndel.length;
						eIdx ++;
						if (eIdx < indelList.size()) 
							eIndel = indelList.get(eIdx);
						else 
							break;
					}
				}


				double v = fragment.getISize();
				//to review if added to global				

				int cover = 0;
				//move to the next tr
				//while (trvIndex < trvList.size() &&						
				//		(trvList.get(trvIndex).getTandemRepeat().compareTo(fragment ) < 0))
				while (trvIndex < trvList.size()){

					while (refIndexTR < refIDs.size() && !trvList.get(trvIndex).getChr().equals(refIDs.get(refIndexTR))){
						refIndexTR ++;
						//System.out.println("advance refIndex Indel " + refIndexStr);
					}
					//assert refIndexStr >= refIDs.size() or trvList.get(trvIndex).getChr().equals(refIDs.get(refIndexStr)
					if (refIndexTR >= refIDs.size()){
						throw new RuntimeException("Sequence " + trvList.get(trvIndex).getChr() + " in STR list not present in the reference genome");
					}
					//assert trvList.get(trvIndex).getChr().equals(refIDs.get(refIndexTR))
					if (refIndexTR > fragment.getReferenceIndex())
						break;

					if (refIndexTR == fragment.getReferenceIndex() &&
							trvList.get(trvIndex).getStart() >= fragment.getStart() )
						break;


					trvIndex ++;
				}
				//assert trv at trvIndex is the first one after fragment.start

				//search for the tr that is encompassed by this insert
				if (refIndexTR == fragment.getReferenceIndex()){
					int tIndex = trvIndex;				
					while (tIndex < trvList.size()){
						TandemRepeat tr = trvList.get(tIndex).getTandemRepeat();

						if (!tr.getChr().equals(refIDs.get(refIndexTR)))
							break;

						if (tr.getStart() > end)
							break;

						if (start < tr.getStart() - gap && end > tr.getEnd() + gap){
							countEff ++;
							counts[tIndex] ++;
							sums[tIndex] += v;
							sqSums[tIndex] += v * v;
							cover ++;
							//if (cover >= 2){
							//	System.out.println(insert + "    " + cover);
							//}
						}		
						tIndex ++;	
					}//while
				}//if
				if (cover == 0){
					tSum += v;
					tSqSum += v * v;
					tCount ++;
				}
			}//while
			
			double gMean = tSum / tCount;//global mean
			double gVar = tSqSum / tCount  - gMean * gMean;//global variance

			if (fd > 0) gMean = fm;
			if (fd > 0) gVar = fd * fd;

			System.out.println("##Lib mean = " + gMean + "  std = " + Math.sqrt(gVar) + " effective " + countEff);

			//now doing inference
			for (int i = 0; i < trvList.size(); i++){
				TandemRepeatVariant trv = trvList.get(i);
				double mean = 0;//, var = 0;
				//		System.out.println("#=#" + i + "\t" + counts[i]);

				if (counts[i] != 0){
					mean = sums[i] / counts[i]; //sample mean
					//var = sqSums[i]  / counts[i] - mean * mean;//sample variance
					double meanP = trv.getMean(), varP = trv.getStd() * trv.getStd();

					double eVar = gVar * varP / (counts[i] * trv.getPeriod() * trv.getPeriod() * varP + gVar);
					double eMean = ((meanP / varP) - (counts[i] * trv.getPeriod() * (mean - gMean) / gVar)) * eVar;


					//	System.out.print("\t" + counts[i]);
					trv.setStd(Math.sqrt(eVar));
					trv.setMean(eMean);
					double predicted = eMean;									
					trv.setVar(predicted);			
					trv.addEvidence(counts[i]);
				}			
				NormalDistribution d = new NormalDistribution(trv.getMean(), trv.getStd());
				trv.setConfidence(d.cumulativeProbability(trv.getVar() - 0.5, trv.getVar() + 0.5));
			}//for i
			if (ansList != null)
				CompareTRV.compareStr(ansList, trvList,0,10000);
		}//for iFdx

		String[] headers = TandemRepeatVariant.STANDARD_HEADERS;
		//{ShortTandemRepeat.chrHd, ShortTandemRepeat.startHd, ShortTandemRepeat.endHd, ShortTandemRepeat.periodHd, ShortTandemRepeat.unitNoHd, STRVariation.varHd, STRVariation.confidenceHd, STRVariation.meanHd, STRVariation.stdHd};

		TandemRepeatVariant.print(trvList, out, headers);

		//out.write(("#H:" + headers[0]).getBytes());
		//for (int i=1; i < headers.length; i++)
		//	out.write(("\t"+headers[i]).getBytes());
		//out.write('\n');
		//for (int i = 0; i < trvList.size(); i++){
		//	out.write((trvList.get(i).toString(headers)).getBytes());			
		//	//if (ansList != null){
		//	//	out.print("\t"+ansList.get(i).getVar());
		//	//	for (int xx = 0; xx < checkVar.length; xx++)
		//	//		out.print("\t"+checkVar[xx][i]);
		//	//}			
		//	out.write('\n');
		//}
	}	



	/**
	 * Return a list of small indels not included in any STR
	 * Note: entries in bcf file are sorted, and include indels only
	 * @param trvList
	 * @param bcfFile
	 * @throws IOException
	 */
	static ArrayList<Indel> filterBCF (ArrayList<TandemRepeatVariant>  trvList, String bcfFile) throws IOException{		
		ArrayList<Indel> listIndels = new ArrayList<Indel>();

		Hashtable<String,String> chrSet = new Hashtable<String,String>();

		int trIndx = 0;

		TandemRepeatVariant trv = null;
		if (trvList.size() > 0)
			trv = trvList.get(trIndx);		

		BufferedReader in = SequenceReader.openFile(bcfFile);
		String line = "";
		while ((line = in.readLine()) != null){
			line = line.trim();
			if (line.startsWith("#")) continue;
			if (line.length() == 0 ) continue;

			String[] toks = line.split("\\t");
			//ref = toks[3]
			//alt = toks[4]


			int refLength = toks[3].length();

			String [] alts = toks[4].split(",");


			int sum = 0, count = 0;

			for (int i = 0; i < alts.length; i++){
				sum += alts[i].length();
				count ++;
			}

			double altLength = sum * 1.0 / count;

			double indelLength = altLength - refLength;
			int position = Integer.parseInt(toks[1]) +  refLength;			

			String chr = chrSet.get(toks[0]);
			if (chr == null){
				chrSet.put(toks[0], toks[0]);
				chr = toks[0];
			}

			Indel indel = new Indel(chr, position, indelLength);

			while ( trv != null && ((trv.getChr().compareTo(chr) < 0) || (trv.getEnd() < position && trv.getChr().equals(chr))) ) {
				trIndx ++;
				if ( trIndx < trvList.size())				
					trv = trvList.get(trIndx);
				else 
					trv = null;
			}
			//
			if (trv != null && (trv.getStart() < position) && (trv.getChr().equals(chr)) ){
				//System.out.println((++aaaa) + "    " + chr+" - " + position + "  " + indelLength + " : " + trv);

			}else{
				listIndels.add(indel);
			}					
		}

		System.out.println(listIndels.size());

		return listIndels;		
	}	
}

class Indel{	
	String chr;
	int start;
	double length;// negative for deletion from ref, positive for insertion to ref.

	public Indel(String c, int s, double indelLength){
		chr = c;
		start = s;
		length = indelLength;
	}
}
