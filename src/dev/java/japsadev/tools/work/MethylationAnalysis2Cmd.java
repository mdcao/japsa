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

package japsadev.tools.work;


import java.io.BufferedReader;
import java.io.IOException;

import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * Filter a bam filem based on some criteria. Input file in bam format assumed
 * to be sorted and indexed
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.dev.methyC2", 
	scriptDesc = "Put the counts together"
	)
public class MethylationAnalysis2Cmd extends CommandLine{
	//CommandLine cmdLine;
	public MethylationAnalysis2Cmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("sample", null, "Sample ID",true);		
		addString("output", null, "name of the output file",true);
		addInt("window", 1000, "window");
		//addInt("filterBits", 0, "Filter reads based on flag. Common values:\n 0    no filter\n 256  exclude secondary alignment \n 1024 exclude PCR/optical duplicates\n 2048 exclude supplementary alignments");



		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		MethylationAnalysis2Cmd cmdLine = new MethylationAnalysis2Cmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String sample  = cmdLine.getStringVal("sample");
		String output = cmdLine.getStringVal("output");
		int window = cmdLine.getIntVal("window");
		

		analyse(sample, window, output);
	}

	static void analyse(String sample, int windowSize, String output) throws IOException{


		String [] chroms = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",  
			"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
			"chr21", "chr22", "chrX", "chrY"};

		SequenceOutputStream fCount = SequenceOutputStream.makeOutputStream(output);

		fCount.print("#CHROM\tID\tSTART\tEND\tFORMAT\tNORMAL\tTUMOUR\n");		

		for (String chrom:chroms){		
			BufferedReader ncBR =  SequenceReader.openFile("count_" + sample + "N-CNTL_S1_" + chrom + ".dat");
			BufferedReader tcBR =  SequenceReader.openFile("count_" + sample + "T-CNTL_S1_" + chrom + ".dat");
			BufferedReader nbBR =  SequenceReader.openFile("count_" + sample + "N-BSC_S1_" + chrom + ".dat");
			BufferedReader tbBR =  SequenceReader.openFile("count_" + sample + "T-BSC_S1_" + chrom + ".dat");

			int start = 0;
			double normalCount = 0, tumourCount = 0;

			int pos = 0;

			String ncLine = "";
			String nbLine, tbLine, tcLine;
			while ((ncLine = ncBR.readLine())!=null){
				nbLine = nbBR.readLine();			
				tcLine = tcBR.readLine();
				tbLine = tbBR.readLine();

				if (ncLine.startsWith("#"))
					continue;

				//double value = 0;

				//tumour control 
				String [] toks = tcLine.split("\t");
				pos = Integer.parseInt(toks[0]);				

				while (pos > start + windowSize){					
					fCount.print(chrom+"\t"+(start+1)+"_"+(start+windowSize)+"\t"+(start+1) + "\t" + (start+windowSize) + "\tDP\t"+((int) (normalCount)) + '\t' + ((int) (tumourCount)) + '\n');
					start += windowSize;
					normalCount = tumourCount = 0;					
				}				

				//if tumourcontrol base is NOT C
				if (Integer.parseInt(toks[1]) == 0){
					//fCount.print(chrom +"\t" + (pos -1 ) + "\t" + pos +"\t0\n");
					continue;
				}

				//normal control
				toks = ncLine.split("\t");
				if (Integer.parseInt(toks[1])==0){
					//					fCount.print(chrom +"\t" + (pos -1 ) + "\t" + pos +"\t0\n");
					continue;
				}


				toks = nbLine.split("\t");
				if (Integer.parseInt(toks[2])==0){
					//					fCount.print(chrom +"\t" + (pos -1 ) + "\t" + pos +"\t0\n");
					continue;
				}

				double normalRatio = Double.parseDouble(toks[1]) / Double.parseDouble(toks[2]);

				toks = tbLine.split("\t");
				if (Integer.parseInt(toks[2])==0){
					//					fCount.print(chrom +"\t" + (pos -1 ) + "\t" + pos +"\t0\n");
					continue;
				}
				double tumourRatio = Double.parseDouble(toks[1]) / Double.parseDouble(toks[2]);

				normalCount += normalRatio;
				tumourCount += tumourRatio;
				//				fCount.print(chrom +"\t" + (pos -1 ) + "\t" + pos +"\t" + (tumourRatio - normalRatio) + '\n');		
			} 
			tbBR.close();
			tcBR.close();
			nbBR.close();
			ncBR.close();

			while (pos > start){				
				fCount.print(chrom+"\t"+(start+1)+"_"+Math.min(pos, start+windowSize)+"\t"+(start+1) + "\t" + (start+windowSize) + "\tDP\t"+((int) (normalCount)) + '\t' + ((int) (tumourCount)) + '\n');
				start += windowSize;
				normalCount = tumourCount = 0;					
			}				
		}


		fCount.close();
	}
}
