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

package japsa.tools.bio.hts;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTextWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.hts.aareads", 
	scriptDesc = "Filter reads supporting alternative alleles"
	)
public class AlternativeAllelesCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(AlternativeAllelesCmd.class);

	//CommandLine cmdLine;
	public AlternativeAllelesCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 


		addString("input", null, "Name of the input bam file",true);
		addString("reference", null, "Reference file",true);
		addString("vcf", null, "Name of the vcf file",true);
		addString("output", "-", "Name of the output file");
		addInt("threshold", 500, "Maximum concordant insert size");

		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		AlternativeAllelesCmd cmdTool = new AlternativeAllelesCmd ();
		args = cmdTool.stdParseLine(args);

		/**********************************************************************/
		String input = cmdTool.getStringVal("input");
		String ref = cmdTool.getStringVal("reference");
		String output = cmdTool.getStringVal("output");
		String vcf = cmdTool.getStringVal("vcf");
		int threshold = cmdTool.getIntVal("threshold");

		addSequence(input, vcf, ref, output,threshold);

	}

	public static class VarRecord{
		String chrom;
		int pos;
		int base; //Alphabet.DNA.A, Alphabet.DNA.C, Alphabet.DNA.G, Alphabet.DNA.T;

		VarRecord(String s, int p, int b){
			chrom = s;
			pos = p;
			base = b;			
		}

		static VarRecord parseLine(String line){
			String [] toks = line.split("\t");
			if (toks.length < 4)
				throw new RuntimeException("Line " + line + " unexpected!!");

			if (toks[3].length() > 1)
				throw new RuntimeException("Field " + toks[3] + " unexpected in line " + line +   "!!");

			int b = -1;
			switch (toks[3].charAt(0)){
			case 'A':
			case 'a': 
				b = Alphabet.DNA.A;
				break;

			case 'C':
			case 'c': 
				b = Alphabet.DNA.C;
				break;

			case 'G':
			case 'g': 
				b = Alphabet.DNA.G;
				break;
			case 'T':
			case 't': 
				b = Alphabet.DNA.T;
				break;

			default:
				throw new RuntimeException("Field " + toks[3] + " unexpected in line " + line +   "!!");

			}

			int p = Integer.parseInt(toks[1]) - 1;
			return new VarRecord(toks[0], p, b);

		}
	}

	/*This class is used to make a single buffered reader of a multichrom VCF look like a series of single chrom VCF files */
	public static class MultiChromVCFReader{
		BufferedReader bf;
		String myChrom;
		boolean finishedChrom = false;
		VarRecord fVar;// = nextRecord(bf);
		MultiChromVCFReader(String vcfFile) throws IOException{
			bf =  SequenceReader.openFile(vcfFile);
			bf.readLine();//dont care the first line
			fVar = AlternativeAllelesCmd.nextRecord(bf);
			
		}
		String updateChrom(){
			myChrom = fVar.chrom;
			finishedChrom = false;
			return myChrom;
		}
		public VarRecord nextRecord()  throws IOException{
			fVar = AlternativeAllelesCmd.nextRecord(bf);
			if(fVar ==null) return null;
			if(!fVar.chrom.equals(myChrom)) finishedChrom = true; // this means we have moved to next chrom
			if(finishedChrom) return null;
			// TODO Auto-generated method stub
			else return fVar;
		}
		public void close() throws IOException{
			bf.close();
			
		}
	}
	
	static VarRecord nextRecord(BufferedReader br) throws IOException{
		String line = br.readLine();
		if (line == null)
			return null;

		return  VarRecord.parseLine(line);

	}

	static void addSequence(String inFile, String vcfFile, String reference,String outFile,  int threshold) throws IOException{
		///////////////////////////////////////////////////////////
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(inFile));						
		SAMFileHeader samHeader = samReader.getFileHeader();		
		SAMTextWriter samWriter = outFile.equals("-")?	(new SAMTextWriter(System.out))	:(new SAMTextWriter(new File(outFile)));
		samWriter.setSortOrder(SortOrder.unsorted, false);		
		samWriter.writeHeader( samHeader.getTextHeader());
		MultiChromVCFReader bf =  new MultiChromVCFReader(vcfFile);
		for(int i=0; bf.fVar!=null; i++){
			if(i > 0){
				 // need to refresh Sam Reader for every chromosome.  I am sure this could be more efficient - LC
				 samReader = SamReaderFactory.makeDefault().open(new File(inFile));						
				 samHeader = samReader.getFileHeader();	
			}
			String myChrom = bf.updateChrom();
			addSequence(inFile, bf, reference,  samReader,samHeader, samWriter, threshold);
			samReader.close();
		}
		bf.close();
		samWriter.close();
		
		///////////////////////////////////////////////////////////
	}
	
	static void addSequence(String inFile, MultiChromVCFReader bf, String reference,   SamReader samReader, SAMFileHeader samHeader, SAMTextWriter samWriter, int threshold) throws IOException{		
		//double sumIZ = 0, sumSq = 0;
		//int countGood = 0, countBad = 0, countUgly = 0;
		//int countALLGood = 0, countALLBad = 0, countALLUgly = 0;
		//double sumALLIZ = 0, sumALLSq = 0;
		//Good: 0 < insert size <= SIZE_THRESHOLD
		//Bad:  insertSize >SIZE_THRESHOLD
		//Ugly: insertSize=0
		
		HashSet<String> somaticSet = new HashSet<String>(); 
		
		

		LinkedList<VarRecord> varList = new  LinkedList<VarRecord>();

		String myChrom = bf.myChrom;

		varList.add(bf.fVar);

	//	String myChrom = fVar.chrom;
		Sequence refSeq = null;
		{
			LOG.info("Read reference started");
			ArrayList<Sequence> seqs = SequenceReader.readAll(reference, Alphabet.DNA());
			for (Sequence seq:seqs){
				if (seq.getName().equals(myChrom)){
					refSeq = seq;
					break;
				}
			}
			LOG.info("Read reference done");
		}

		if(refSeq == null){
			bf.close();
			throw new RuntimeException("Chrom " + myChrom + " not found in the reference!!");
		}

		boolean hasVar = true;

		

		int myChromIndex = samHeader.getSequenceIndex(myChrom);
		if(myChromIndex < 0){
			samWriter.close();
			samReader.close();
			bf.close();
			throw new RuntimeException("Chrom " + myChrom + " not found in bam file!!");
		}

		SAMRecordIterator samIter = samReader.query(refSeq.getName(),0,0,false);
		LOG.info(" " + samIter.hasNext());

		LOG.info("Chrom = " + myChrom + " RefName = " + refSeq.getName() + " chromIndex = " + myChromIndex);

		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();
			if (sam.getReadUnmappedFlag())
				continue;

			int samRefIndex = sam.getReferenceIndex(); 
			if (samRefIndex < myChromIndex)
				continue;

			if (samRefIndex > myChromIndex)
				break;//while

			//assert samRefIndex == myChromIndex
			
			//int insertSize = Math.abs(sam.getInferredInsertSize());
			
			//if (insertSize == 0){
			//	countALLUgly ++;
			//}else if (insertSize <= threshold){
			//	countALLGood ++;
			//	sumALLIZ += insertSize;
			//	sumALLSq += insertSize * insertSize;					
			//}else{
			//	countALLBad ++;
			//}
			
			String readName = sam.getReadName();
			if (somaticSet.contains(readName))
				continue;

			Sequence readSeq = new Sequence(Alphabet.DNA(), sam.getReadString(), sam.getReadName());
			boolean support = false;

			int readPos = 0;//start from 0					
			int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index
			for (final CigarElement e : sam.getCigar().getCigarElements()) {
				final int  length = e.getLength();
				switch (e.getOperator()) {
				case H :					
					break; // ignore hard clips
				case P : 					
					break; // ignore pads	                
				case S :
					readPos += length;
					break; // soft clip read bases	                	
				case N : 
					refPos += length;					
					break;  // reference skip

				case D ://deletion      	
					refPos += length;
					break;

				case I :	                	
					readPos += length;
					break;
				case M :
					for (int i = 0; i < length; i++){
						int readBase = readSeq.getBase(readPos + i);
						if (refSeq.getBase(refPos + i) != readBase){
							//1. 
							while(varList.size() > 0){
								VarRecord first = varList.getFirst();

								if (first.pos < sam.getAlignmentStart()){
									varList.removeFirst();
									continue;
								}
								break;
							}

							//2. go through the list
							int currentVarPos = -1;
							for (VarRecord var:varList){

								if (var.pos == refPos + i && var.base == readBase){
									//yay
									support = true;
									break;//for									
								}

								currentVarPos = var.pos;								
								if (currentVarPos > refPos + i)
									break;//for
							}							
							if (support)
								break;//for i

							while (currentVarPos < refPos + i && hasVar){
								VarRecord var = bf.nextRecord();
								if (var == null){
									hasVar = false;
									break;
								}
								varList.add(var);

								if (var.pos == refPos + i && var.base == readBase){
									//yay
									support = true;
									break;//for									
								}

								currentVarPos = var.pos;
							}							
						}//if
						if (support)
							break;//for

					}//for

					readPos += length;
					refPos  += length;
					break;

				case EQ :
					readPos += length;
					refPos  += length;

					break;
				case X :
					//do some thing here
					LOG.error("Var X is not currently support, please let Minh know if you see this");
					readPos += length;
					refPos  += length;
					break;
				default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
				}//case
				if (support)
					break;
			}//for

			if (support){
				//samWriter.writeAlignment(sam);				
				//if (insertSize == 0){
				//	countUgly ++;
				//}else if (insertSize <= threshold){
				//	countGood ++;
				//	sumIZ += insertSize;
				//	sumSq += insertSize * insertSize;					
				//}else{
				//	countBad ++;
				//}
				somaticSet.add(readName);
			}

		}//while
		
		samIter.close();
		
		samIter = samReader.query(refSeq.getName(),0,0,false);
		while (samIter.hasNext()){			
			SAMRecord sam = samIter.next();
			String readName = sam.getReadName();
			if (somaticSet.contains(readName)){
				samWriter.writeAlignment(sam);
			}
		}
		
		
	//	bf.close();
		/**********************************************************************
		
		System.out.println("================ ALL DATA===================");
		System.out.printf("Good insert fragments  (0<insert<=%d): %d\n", threshold,countALLGood);
		if (countALLGood>0){
			double mean  = sumALLIZ / countALLGood;
			double stdev = Math.sqrt(sumALLSq/countALLGood - mean * mean );			
			System.out.printf("  mean = %f, std=%f\n",mean, stdev);	
		}
		System.out.printf("Bad  insert fragments  (insert>%d): %d\n", threshold,countALLBad);
		System.out.printf("Ungly insert fragments (insert=0): %d\n", countALLUgly);
		
		
		System.out.println("================ SELECTED DATA===================");
		System.out.printf("Good insert fragments  (0<insert<=%d): %d\n", threshold,countGood);
		if (countGood>0){
			double mean  = sumIZ / countGood;
			double stdev = Math.sqrt(sumSq/countGood - mean * mean );			
			System.out.printf("  mean = %f, std=%f\n",mean, stdev);	
		}
		System.out.printf("Bad  insert fragments  (insert>%d): %d\n", threshold,countBad);
		System.out.printf("Ungly insert fragments (insert=0): %d\n", countUgly);
		
		/**********************************************************************/
	}


}
