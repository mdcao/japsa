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
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
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

		double sumDiff=0;
		double sumSq=0;
		int count=0;
		
		String in;
		public String toString(){
			return in+"\t"+(double)sumDiff/(double) count;
		}
		VarRecord(String s, int p, int b, String line){
			this.in = line;
			chrom = s;
			pos = p;
			base = b;			
		}
		void print(PrintWriter pw, String left, String right){
			double mean = (double)sumDiff/(double) count;
			double sd = Math.sqrt(sumSq/(double)count - Math.pow(mean,2));
			pw.print(in);
			pw.print("\t");
			pw.print(mean);
			pw.print("\t");
			pw.print(sd);
			pw.print("\t");
			pw.print(count);
			pw.print("\t");
			pw.print(left+"\t"+right);
			pw.println();
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
			return new VarRecord(toks[0], p, b, line);

		}

		public void add(int diff) {
			this.sumDiff+=diff;
			this.sumSq+=Math.pow((double) diff, 2);
			this.count++;
			// TODO Auto-generated method stub
			
		}
	}

	static VarRecord nextRecord(BufferedReader br) throws IOException{
		String line = br.readLine();
		if (line == null)
			return null;

		return  VarRecord.parseLine(line);

	}

	static class FragInfo{
		public FragInfo(int alignmentStart) {
			this.start = alignmentStart;
		}
		int start;
		int end;
		List<VarRecord> snps  = new ArrayList<VarRecord>();
		public void setEnd(int alignmentEnd) {
			this.end = alignmentEnd;
		}
		public String toString(){
			String st = this.start+"";
			if(snps!=null) st = st+":"+snps.toString();
			return st;
		}
		public void addSNP(VarRecord var) {
			snps.add(var);
			
		}
		public void complete() {
			int diff = end - start;
			for(int i=0; i<snps.size(); i++){
				snps.get(i).add(diff);
			}
		}
	}
	static ConcurrentMap<String, FragInfo> readToPos = new ConcurrentHashMap<String, FragInfo>();// maps read pair to position
	
	
	/*static class Fragments{
		
		public Fragments(Integer diff, String readName) {
			// TODO Auto-generated constructor stub
		}
		int sumFrag;
		int sum;
	}*/
	
	static void clearUpTo( PrintWriter pw, LinkedList<VarRecord> varList, int pos ,String chrom, Sequence ref){
		int i=0; 
		for (i=0; i<varList.size(); i++){
			VarRecord var = varList.get(i);
			if(!var.chrom.equals(chrom)) throw new RuntimeException("!!");
			if(var.pos>pos) break;
			var.print(pw, ref.subSequence(var.pos-10, var.pos).toString(), ref.subSequence(var.pos, var.pos+10).toString());
		}
		pw.flush();
		for(int j=i-1; j>=0;j--){
			varList.remove(j);
		}
		
		
	}
	static void clearUpTo(ConcurrentMap<String, FragInfo> readToPos2,int pos ,String chrom, PrintWriter ls){
	//	Set<String> torem = new HashSet<String>();
		for(Iterator<String> it = readToPos2.keySet().iterator(); it.hasNext();){
			String key = it.next();
			FragInfo val = readToPos2.get(key);
			if(val.start < pos){
			//	System.err.println("removing "+key+" "+val);
				ls.println(key);
				readToPos2.remove(key);
			}
		}
	}
	
	static void addSequence(String inFile, String vcfFile, String reference, String outFile, int threshold) throws IOException{		
		//double sumIZ = 0, sumSq = 0;
		//int countGood = 0, countBad = 0, countUgly = 0;
		//int countALLGood = 0, countALLBad = 0, countALLUgly = 0;
		//double sumALLIZ = 0, sumALLSq = 0;
		//Good: 0 < insert size <= SIZE_THRESHOLD
		//Bad:  insertSize >SIZE_THRESHOLD
		//Ugly: insertSize=0
		
		Set<String> somaticSet = new HashSet<String>(); 


		BufferedReader bf =  SequenceReader.openFile(vcfFile);
		
		bf.readLine();//dont care the first line

		LinkedList<VarRecord> varList = new  LinkedList<VarRecord>();

		VarRecord fVar = nextRecord(bf);

		varList.add(fVar);

		final String myChrom = fVar.chrom;
		PrintWriter vcf_out = new PrintWriter(new FileWriter(new File(vcfFile+".out.vcf")));
		PrintWriter ls = new PrintWriter(new FileWriter(new File("unmatched_reads.txt")));

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

		///////////////////////////////////////////////////////////
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(inFile));						


		SAMFileHeader samHeader = samReader.getFileHeader();		
		SAMTextWriter samWriter = outFile.equals("-")?	(new SAMTextWriter(System.out))	:(new SAMTextWriter(new File(outFile)));

		samWriter.setSortOrder(SortOrder.unsorted, false);		
		samWriter.writeHeader( samHeader.getTextHeader());	
		///////////////////////////////////////////////////////////

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
			
			
			clearUpTo(vcf_out,varList,sam.getAlignmentStart()-10000, myChrom , refSeq); // clear up varlist up to 10,000 bases before current
			clearUpTo(readToPos, sam.getAlignmentStart()-10000, myChrom, ls);
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
			FragInfo pairedPos = readToPos.remove(readName);
		//	boolean firstPair = sam.getFirstOfPairFlag();
		//	boolean secondPair = sam.getSecondOfPairFlag();
		//	boolean neg = sam.getReadNegativeStrandFlag();
			if(sam.getPairedReadName()!=null){
				if(pairedPos==null){
					readToPos.put(readName, pairedPos=new FragInfo(sam.getAlignmentStart()));
				}
				else{
					pairedPos.setEnd(sam.getAlignmentEnd());
					pairedPos.complete();
					//System.err.println("DIFFF "+ diff);
				}
			}
			if(pairedPos==null) continue;
		//	if (somaticSet.contains(readName))
		//		continue;

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
							/*while(varList.size() > 0){
								VarRecord first = varList.getFirst();

								if (first.pos < sam.getAlignmentStart()){
									varList.removeFirst();
									continue;
								}
								break;
							}*/

							//2. go through the list
							int currentVarPos = -1;
							for (VarRecord var:varList){

								if (var.pos == refPos + i && var.base == readBase){
									//yay
									support = true;
									pairedPos.addSNP(var);
									break;//for									
								}

								currentVarPos = var.pos;								
								if (currentVarPos > refPos + i)
									break;//for
							}							
							if (support)
								break;//for i

							while (currentVarPos < refPos + i && hasVar){
								VarRecord var = nextRecord(bf);
								if (var == null){
									hasVar = false;
									break;
								}
								varList.add(var);

								if (var.pos == refPos + i && var.base == readBase){
									//yay
									support = true;
									pairedPos.addSNP(var);
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
			//	if (support)	break;
				
			}//for

			if (support && somaticSet!=null){
				//samWriter.writeAlignment(sca$240am);				
				//if (insertSize == 0){
				//	countUgly ++;
				//}else if (insertSize <= threshold){
				//	countGood ++;
				//	sumIZ += insertSize;
				//	sumSq += insertSize * insertSize;					
				//}else{
				//	countBad ++;
				//}
				//pairedPos.
				somaticSet.add(readName);
			}

		}//while
		
		samIter.close();
		
		if(somaticSet!=null){
			samIter = samReader.query(refSeq.getName(),0,0,false);
			while (samIter.hasNext()){			
				SAMRecord sam = samIter.next();
				String readName = sam.getReadName();
				if (somaticSet.contains(readName)){
					samWriter.writeAlignment(sam);
				}
			}
			
			
			samWriter.close();
			samReader.close();
		}
		bf.close();
		clearUpTo(vcf_out,varList,refSeq.length(), myChrom , refSeq); // clear up varlist up to 10,000 bases before current

		vcf_out.close();
		ls.close();
		
		
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
