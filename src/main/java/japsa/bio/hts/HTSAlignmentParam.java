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

package japsa.bio.hts;


import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.HTSUtilities;
import japsa.util.JapsaMath;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.hts.alignOpt", 
scriptDesc = "Parameter estimation for alignment of erronenous read data")
public class HTSAlignmentParam {
	public static void main(String [] args) throws IOException, InterruptedException{
		/*********************** Setting up script ****************************/
		Deployable annotation = HTSAlignmentParam.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options] + 1.fq + [2.fq]" , annotation.scriptDesc());		
		/**********************************************************************/	
		//cmdLine.addString("bamFile", null,  "Name of bam file", true);
		cmdLine.addString("reference", null, "Name of reference genome",true);

		cmdLine.addString("aligner", "bwa", "Alignment program, currently support bwa (bwa mem) and blasr");
		cmdLine.addString("index", null, "Index of the aligner",true);
		cmdLine.addString("prefix", "tmp", "Prefix");

		cmdLine.addInt("seedLength", 14, "Seed length");
		cmdLine.addInt("thread", 8, "Number of thread");

		cmdLine.addInt("qual", 0, "Minimum quality required");

		args = cmdLine.stdParseLine_old(args);			
		/**********************************************************************/
		//String bamFile = cmdLine.getStringVal("bamFile");
		reference = cmdLine.getStringVal("reference");		
		indexFile = cmdLine.getStringVal("index");
		prefix = cmdLine.getStringVal("prefix");

		seedLength   =  cmdLine.getIntVal("seedLength");
		thread   =  cmdLine.getIntVal("thread");		
		qual = cmdLine.getIntVal("qual");

		String aligner = cmdLine.getStringVal("aligner");

		if (args.length > 0)
			fq1 = args[0];

		if (args.length > 1)
			fq2 = args[1];

		if (fq1 == null){
			System.err.println("Specify a file for reads" + 
					cmdLine.errors());
			System.exit(1);
		}

		if (aligner.equals("bwa")){
			int
			copyScore = 1, //A		
			mmPen = 1,//B 2?

			delPen = 1,
			insPen = 1,//O 2?

			insExt = 1,//E
			delExt = 1;

			for (int i = 0; i < 10;i++){
				double [] params = runBWA(copyScore,mmPen, insPen, delPen,insExt,delExt,i);

				copyScore = (int) Math.round(params[0] * 3.0/ params[0]); //essentially =3.0
				mmPen = (int) Math.round(params[1]* 3.0/ params[0]);
				insPen = (int) Math.round(params[2]* 3.0/ params[0]);
				delPen = (int) Math.round(params[3]* 3.0/ params[0]);
				insExt = (int) Math.round(params[4]* 3.0/ params[0]);
				delExt = (int) Math.round(params[5]* 3.0/ params[0]);
			}		
		}else if (aligner.equals("bwaf")){
			double
			copyScore = 1, //A		
			mmPen = 2,//B

			delPen = 2,
			insPen = 2,//O

			insExt = 1,//E
			delExt = 1;			

			for (int i = 0; i < 10;i++){
				double [] params = runBWAF(copyScore,mmPen, insPen, delPen,insExt,delExt,i);
				//{-costCopy, costChange, costIns, costDel, costIE, costDE};
				copyScore = params[0];
				mmPen   =  params[1];
				insPen  =  params[2];
				delPen  =  params[3];
				insExt  =  params[4];
				delExt  =  params[5];			
			}			
		}else {
			System.err.println("Unknown program " + aligner);
			System.exit(1);
		}
		//paramEst(bamFile, reference, qual);
	}
	static String fq1 = null, fq2 = null, reference = null, indexFile = null;
	static String prefix;

	static int thread = 4, seedLength = 14, qual=0;

	static double [] runBWA(int copyScore, int mmPen, int insPen, int delPen, int insExt, int delExt, int iter) 
			throws IOException, InterruptedException{		

		String scriptFile = prefix + "runBWA.sh";
		String resultFile = prefix + "result"+iter+".sam";
		String resultLog = prefix + "result"+iter+".log";

		PrintStream ps = new PrintStream(new FileOutputStream(scriptFile));
		String bwaCommand = "/sw/bwa/current/bin/bwa mem -t " + thread 
				+ " -k " + seedLength 
				+ " -A " + copyScore 
				+ " -B " +  mmPen
				+ " -O " +  delPen+","+insPen
				+ " -E " +  delExt+","+insExt
				+ " -W20 -r10 -L0"
				+ " " + indexFile 
				+ " " + fq1 + ((fq2==null)?"":(" " +fq2))				
				+ " > " + resultFile + " 2> " + resultLog;
		ps.println(bwaCommand);
		ps.close();		

		Logging.info("Running " + bwaCommand);
		Process process = Runtime.getRuntime().exec("bash ./" + scriptFile);

		process.waitFor();
		double [] ret = paramEst(resultFile, reference, qual); 
		return ret;
	}

	static double [] runBWAF(double copyScore, double mmPen, double insPen, double delPen, double insExt, double delExt, int iter) 
			throws IOException, InterruptedException{		

		DecimalFormat df = new DecimalFormat("#.##");

		String scriptFile = prefix + "runBWA.sh";
		String resultFile = prefix + "result"+iter+".sam";
		String resultLog = prefix + "result"+iter+".log";

		PrintStream ps = new PrintStream(new FileOutputStream(scriptFile));
		String bwaCommand = "/home/minhduc/.usr/local/bin/bwaf mem -L0 -t " + thread 
				+ " -k " + seedLength 
				+ " -A " +  df.format(copyScore) 
				+ " -B " +  df.format(mmPen)
				+ " -O " +  df.format(delPen)+","+df.format(insPen)
				+ " -E " +  df.format(delExt)+","+df.format(insExt)
				+ " " + indexFile 
				+ " " + fq1 + ((fq2==null)?"":(" " +fq2))				
				+ " > " + resultFile + " 2> " + resultLog;
		ps.println(bwaCommand);
		ps.close();		

		Logging.info("Running " + bwaCommand);
		Process process = Runtime.getRuntime().exec("bash ./" + scriptFile);

		process.waitFor();

		double [] ret = paramEst(resultFile, reference, qual); 
		return ret;
	}

	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 * @param bamFile
	 * @param pad
	 * @throws IOException 
	 */
	static double[] paramEst(String bamFile, String refFile, int qual) throws IOException{	

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFile));	

		SAMRecordIterator samIter = samReader.iterator();
		//Read the reference genome
		ArrayList<Sequence> genomes = SequenceReader.readAll(refFile, Alphabet.DNA());

		//get the first chrom
		int currentIndex = 0;
		Sequence chr = genomes.get(currentIndex);

		long    totBaseIns = 0,
				totBaseDel = 0,
				totNumIns = 0,
				totNumDel = 0,
				totMisMatch = 0,
				totMatch = 0;
		long totReadBase = 0, totRefBase = 0;
		int  numReads = 0, numReadsConsidered = 0;


		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();
			numReads ++;
			if (sam.getReadUnmappedFlag())
				continue;

			if (sam.getMappingQuality() < qual)
				continue;			

			
			//int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index
			int refIndex = sam.getReferenceIndex();

			//if move to another chrom, get that chrom
			if (refIndex != currentIndex){
				currentIndex = refIndex;
				chr = genomes.get(currentIndex);
			}
			
			//make the read seq			
			Sequence readSeq = new Sequence(Alphabet.DNA(), sam.getReadString(), sam.getReadName());
			if (readSeq.length() <= 1){
				Logging.warn(sam.getReadName() +" ignored");
				continue;
			}
			
			japsa.util.HTSUtilities.IdentityProfile profile = 
					HTSUtilities.identity(chr, readSeq, sam);			

			
			totBaseIns  += profile.baseIns;
			totBaseDel  += profile.baseDel;
			totNumIns   += profile.numIns;
			totNumDel   += profile.numDel;			
			totMisMatch += profile.mismatch;
			totMatch    += profile.match;

			totReadBase += profile.readBase;
			totRefBase  += profile.refBase;

			numReadsConsidered ++;
		}		
		samReader.close();
		//Done
		System.out.println("===================================================");
		System.out.println(numReads + "  " + numReadsConsidered);		

		System.out.println("Deletion " + totBaseDel + " " + totNumDel +" " + totBaseDel*1.0/totRefBase);
		System.out.println("Insertion " + totBaseIns + " " + totNumIns+" " + totBaseIns*1.0/totRefBase);
		System.out.println("MisMatch " + totMisMatch +" " + totMisMatch*1.0/totRefBase);
		System.out.println("Match " + totMatch);

		System.out.println("ReadBase " + totReadBase);
		System.out.println("ReferenceBase " + totRefBase);	

		double totState0 = totMatch + totMisMatch;//
		double probDel = (totNumDel + 1.0) / (totState0 + 3.0);
		double probIns = (totNumIns + 1.0) / (totState0 + 3.0);

		double probMatch = 1.0 - probDel - probIns;
		double probCopy =  (totMatch + 1.0) / (totState0 + 2.0);
		double probChange = (totMisMatch + 1.0) / (totState0 + 2.0);

		double probDE = (1.0 + totBaseDel - totNumDel) / (2.0 +totBaseDel);
		double probIE = (1.0 + totBaseIns - totNumIns) / (2.0 +totBaseIns);

		double costMatch = - JapsaMath.log2(probMatch);		
		double costCopy =  - costMatch + JapsaMath.log2(probCopy) + 2;//minus 2
		double costChange = costMatch - JapsaMath.log2(probChange) + JapsaMath.log2(3.0) - 2;

		double costDel = - JapsaMath.log2(probDel);//plus 2 minus 2
		double costIns = - JapsaMath.log2(probIns);//plus 2 minus 2

		double costDE = - JapsaMath.log2(probDE);//plus 2 minus 2
		double costIE = - JapsaMath.log2(probIE);//plus 2 minus 2		


		System.out.printf("Indentity %f %f %f %f\n",1.0 *totMatch/(totMatch + totMisMatch + totBaseDel +totBaseIns),
				1.0 *totMisMatch/(totMatch + totMisMatch + totBaseDel +totBaseIns),
				1.0 *totBaseIns/(totMatch + totMisMatch + totBaseDel +totBaseIns),
				1.0 *totBaseDel/(totMatch + totMisMatch + totBaseDel +totBaseIns ));

		System.out.printf("Probs %f %f %f %f %f %f\n",probCopy, probChange, probIns, probDel, probIE, probDE);
		System.out.printf("Probs %f %f %f %f %f %f\n",probCopy * probMatch, probChange* probMatch, probIns, probDel, probIE, probDE);
		System.out.printf("Costs %f %f %f %f %f %f\n",costCopy, costChange, costIns, costDel, costIE, costDE);
		return new double[] {costCopy,  costChange, costIns , costDel, costIE , costDE, totReadBase, totMatch, totMisMatch, totNumIns, totNumDel, totBaseIns, totBaseDel};
	}

}
