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
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.zip.GZIPOutputStream;

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
static boolean writeSAM=false;
static int extraLength=20;
public static boolean partial=true;
	//CommandLine cmdLine;
	public AlternativeAllelesCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 


		addString("input", null, "Name of the input bam file",false);
		addString("reference", null, "Reference file",true);
		addBoolean("writeSAM",false, "write sam output",false);
		addString("vcf_plasma", null, "Name of the vcf file",true);
		addString("vcf_tumor", null, "Name of the vcf file tumor",true);
		addString("output", "results", "Name of the output directory");
		addString("chrom", "chr1", "which chromosome",true);
		addInt("threshold", 2000, "Maximum concordant insert size");
		//addInt("maxIns", 10000, "Maximum concordant insert size");

		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		AlternativeAllelesCmd cmdTool = new AlternativeAllelesCmd ();
		args = cmdTool.stdParseLine(args);

		/**********************************************************************/
		String input = cmdTool.getStringVal("input");
		String ref = cmdTool.getStringVal("reference");
		String output = cmdTool.getStringVal("output");
		String vcf = cmdTool.getStringVal("vcf_plasma");
		String vcf1 = cmdTool.getStringVal("vcf_tumor");
		String myChrom = cmdTool.getStringVal("chrom");
		writeSAM = cmdTool.getBooleanVal("writeSAM");
		 sizeThreshold = cmdTool.getIntVal("threshold");
		maxIns = sizeThreshold * 2;
		String percentiles_ = "0.1:0.25:0.5:0.75:0.9";
		String[] percs = percentiles_.split(":");
		percentiles = new double[percs.length];
		for(int i=0; i<percentiles.length; i++){
			percentiles[i] = Double.parseDouble(percs[i]);
		}
		File vcf1_ = vcf1==null ? null : new File(vcf1);
		if(!vcf1_.exists()) vcf1_ = null;
		addSequence(input, new File(vcf), vcf1_,ref, output, myChrom);

	}
	
	static int sizeThreshold = 2000;
	static int maxIns = 4000;
	static double[] percentiles =null;// new double[] {0.25, 0.5, 0.75};
static boolean noBAM = false;
	public static class VarRecord{
		String chrom;
		int pos;
		int base; //Alphabet.DNA.A, Alphabet.DNA.C, Alphabet.DNA.G, Alphabet.DNA.T;

		
		double minLength = Double.NaN;
		double maxLength = Double.NaN;
		double sumDiff=0;
		double sumSq=0;
		int count=0;
		int count1 =0;
		
		int[][] vals = new int[2][2];
		
		
		public void add(int[] b) {
			for(int i=0; i<b.length; i++){
				vals[i][0] += b[i];
				vals[i][1] += (int) Math.pow(b[i], 2);
			}
			count1++;
		}
		
		
		String in;
		public String toString(){
			return in+"\t"+(double)sumDiff/(double) count;
		}
		
		double[] percVals;
		VarRecord(String s, int p, int b, String line){
			this.in = line;
			this.percVals = new double[percentiles.length];
			chrom = s;
			pos = p;
			base = b;			
		}
		
		public void calcFeatures(){
			//Collections.sort(lens);
			this.count = lens.size();
			if(count>0){
				minLength = lens.firstKey();
				maxLength = lens.lastKey();
			}
			if(lens.size()>0){
				medianQ(percentiles, percVals, this.lens);
			}else{
				Arrays.fill(percVals, Double.NaN);
			}
		}
		
		
		static String[] extra;
		{
			String[] extra1 = 	new String[] {"tumor", 
				"meanLen", "sdLen","lenMin", "lenMax", "fragCount","left","right","meanBaseQ", "sdBaseQ", "meanMapQ","sdMapQ", "count","count_frac","diff_prev","diff_nxt"};
			String[] extra2 = new String[percentiles.length];
			for(int i=0; i<extra2.length; i++){
				extra2[i] = "lenQ_"+percentiles[i];
			}
			String extra3 = combine(extra1,":")+":"+combine(extra2,":");
			extra = extra3.split(":");
		}
		String print( String left, String right){
			this.calcFeatures();
			//if(noBAM) return in;
			StringBuffer sb = new StringBuffer();
			double mean = (double)sumDiff/(double) count;
			double sd = Math.sqrt(sumSq/(double)count - Math.pow(mean,2));
			sb.append(in);sb.append("\t");
			sb.append(this.tumor!=null);sb.append("\t");
			sb.append(mean);sb.append("\t");
			sb.append(sd);sb.append("\t");
			sb.append(minLength); sb.append("\t");
			sb.append(maxLength); sb.append("\t");
			sb.append(count);	sb.append("\t");
			sb.append(left+"\t"+right);
			for(int i=0; i<vals.length; i++){
				double mean1 = (double)vals[i][0]/(double) count1;
				double sd1 = Math.sqrt((double) vals[i][1]/(double)count1 - Math.pow(mean1,2));
				sb.append("\t");
				sb.append(mean1);
				sb.append("\t");
				sb.append(sd1);
			}
			sb.append("\t");
			sb.append(count1);sb.append("\t");
			sb.append((double)count/(double)count1);sb.append("\t");
			sb.append(this.diff_prev);	sb.append("\t");
			sb.append(this.diff_nxt);  
			for(int i=0; i<percentiles.length; i++){
			sb.append("\t");sb.append(this.percVals[i]);
			}

//			if(this.tumor!=null){
//				sb.append("\t");
//				sb.append(this.tumor.in);
//			}
			return sb.toString();
			//pw.println();
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
		SortedMap<Integer, Integer> lens = new TreeMap<Integer,Integer>();
		
		
		int zeroCount=0;
		int aboveThreshCnt=0;
		public void add(Integer diff) {
			if(diff==0){
				zeroCount++;
			}else if(diff>sizeThreshold){
				aboveThreshCnt++;
			}
			Integer v  = lens.get(diff);
			lens.put(diff, v==null ? 1 : v+1);
			this.sumDiff+=diff;
			this.sumSq+=Math.pow((double) diff, 2);
			this.count++;
			/*
			if(count==0){
				minLength = diff;
				maxLength=diff;
			}else{
				minLength = Math.min(diff, minLength);
				maxLength = Math.max(diff, minLength);
			}
			
			this.count++;
			*/
			// TODO Auto-generated method stub
			
		}
		private static double interpolate(int cov0, int cov, double perc0, double perc, double d) {
			// TODO Auto-generated method stub
			return (double) cov0 + ((d-perc0)/ (perc-perc0)) * ((double)cov-(double) cov0);
		}
		
		public static void medianQ(double[] percentiles, double[] vals, SortedMap<Integer, Integer>map){
			double cumul0=0;
			double cumul1 =0;
			double total = map.values().stream().mapToInt(Integer::intValue).sum();
			 int mapq0=0;
			Iterator<Integer> it = map.keySet().iterator();
			while(it.hasNext()){
				Integer key = it.next();
				Integer cnt = map.get(key);
				int mapq1= key;
				cumul1 = cumul0+cnt;
				double perc0 = cumul0/total;
				double perc = cumul1/total;
				for(int j=0; j<percentiles.length; j++){
					if(percentiles[j] >= perc0  && percentiles[j] <=perc){
						vals[j]  = interpolate(mapq0, mapq1, perc0, perc, percentiles[j]);
					}
				}
				mapq0 = mapq1;
				cumul0 = cumul1;
			}
		}
		
		
		int diff_nxt=-1;
		int diff_prev=-1;
		
		public void setNext(VarRecord nxt) {
			if(nxt.chrom.equals(chrom)){
				diff_nxt = nxt.pos-pos;
				//System.err.println(diff_nxt);
			}
			// TODO Auto-generated method stub
			
		}
		public void setPrevious(VarRecord prev) {
			if(prev.chrom.equals(chrom)){
				diff_prev = pos - prev.pos;
		//		System.err.println(diff_prev);
			}
			
		}
		VarRecord tumor = null;
		public void setMatchedTumour(VarRecord varRecord) {
			tumor = varRecord;
			
		}
		
		
	}

	static VarRecord nextRecord(BufferedReader br, VarRecord previous, String myChrom) throws IOException{
		String line = br.readLine();
		if (line == null)
			return null;
		
		VarRecord nxt =  VarRecord.parseLine(line);
		if(previous==null){
			while(!nxt.chrom.equals(myChrom)){
				line = br.readLine();
				nxt =  VarRecord.parseLine(line);
			}
		}
		if(!nxt.chrom.equals(myChrom)){
			return null;
		}
		if(previous!=null){
			previous.setNext(nxt);
			nxt.setPrevious(previous);
		}
		return nxt;

	}

	/*static class FragInfo{
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
		public int fraglength() {
			// TODO Auto-generated method stub
			return end-start;
		}
	}*/
	//static Map<String, FragInfo> readToPos = new HashMap<String, FragInfo>();// maps read pair to position
	
	
	/*static class Fragments{
		
		public Fragments(Integer diff, String readName) {
			// TODO Auto-generated constructor stub
		}
		int sumFrag;
		int sum;
	}*/
	
	static void clearUpTo( PrintWriter pw, List<VarRecord> varList, int pos ,String chrom, Sequence ref){
		int i=0; 
		for (i=0; i<varList.size(); i++){
			VarRecord var = varList.get(i);
			if(!var.chrom.equals(chrom)) throw new RuntimeException("!!");
			if(var.pos>pos) break;
			String st = var.print(ref.subSequence(var.pos-extraLength, var.pos).toString(), ref.subSequence(var.pos, var.pos+extraLength).toString());
			if(CHECK){
				String[] str = st.split("\t");
				if(str.length!=no_cols) throw new RuntimeException("wrong number cols "+str.length+" "+no_cols);
			}
			pw.println(st);
		}
		pw.flush();
		for(int j=i-1; j>=0;j--){
		//	System.err.println("removing up to "+pos+"  -> "+varList.get(j));
			varList.remove(j);
		}
		
		
	}
	static int unMatchedCount=0;
	static boolean printUnmatched=false;
	/*static void clearUpTo(Map<String, FragInfo> readToPos2,int pos ,String chrom, PrintWriter ls){
		//Set<String> torem = new HashSet<String>();
		for(Iterator<String> it = readToPos2.keySet().iterator(); it.hasNext();){
			String key = it.next();
			FragInfo val = readToPos2.get(key);
			if(val.start < pos){
			//	System.err.println("removing "+key+" "+val);
				if(printUnmatched) ls.println(key);
				unMatchedCount++;
				//torem.add(key);
				it.remove();
//				readToPos2.remove(key);
			}
		}
		//readToPos2
	}*/
	
	static List<String> header = new ArrayList<String>();
	static List<String> headerT = new ArrayList<String>();
	static int no_cols;
	static boolean CHECK=false;
	static void addSequence(String inFile, File vcfFile, File vcfFile1, String reference, String outFile, 
			String myChrom) throws IOException{		
		//double sumIZ = 0, sumSq = 0;
		//int countGood = 0, countBad = 0, countUgly = 0;
		//int countALLGood = 0, countALLBad = 0, countALLUgly = 0;
		//double sumALLIZ = 0, sumALLSq = 0;
		//Good: 0 < insert size <= SIZE_THRESHOLD
		//Bad:  insertSize >SIZE_THRESHOLD
		//Ugly: insertSize=0
		
		File outDir = new File(outFile);
		outDir.mkdir();
		
		Set<String> somaticSet = writeSAM ?  new HashSet<String>() : null;


		BufferedReader bf =  SequenceReader.openFile(vcfFile.getAbsolutePath());
		BufferedReader bf_T =  vcfFile1!=null && vcfFile1.exists() ? SequenceReader.openFile(vcfFile1.getAbsolutePath()) : null;

		
		
		header.addAll(Arrays.asList(bf.readLine().split("\t")));//dont care the first line
		if(bf_T!=null) 		headerT.addAll(Arrays.asList(bf_T.readLine().split("\t")));//dont care the first line

	
		
		List<VarRecord> varList = new  ArrayList<VarRecord>();
		Map<Integer, VarRecord> varList_tumour = new  HashMap<Integer, VarRecord>();
		
		if(bf_T!=null){
			VarRecord prevVar_T = null;
			VarRecord fVar_T = nextRecord(bf_T, prevVar_T, myChrom);
			while(fVar_T!=null){
	//		System.err.println("adding "+fVar_T);
				varList_tumour.put(fVar_T.pos, fVar_T);
				prevVar_T = fVar_T;
				fVar_T = nextRecord(bf_T, prevVar_T, myChrom);
			}
		}
		
		
		
	//	final String myChrom = fVar.chrom;
		OutputStream os = new GZIPOutputStream(new FileOutputStream(new File(outDir, vcfFile.getName()+"."+myChrom+".out.vcf.gz")));
		PrintWriter vcf_out = new PrintWriter(new OutputStreamWriter(os));
		PrintWriter ls = null;
		if(printUnmatched) ls = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(outDir, vcfFile.getName()+".unmatched_reads.txt.gz")))));

		
		
		
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
		
		File sFile = inFile==null ? null : new File(inFile);
		noBAM = sFile==null || !sFile.exists();
	
		header.addAll(Arrays.asList(VarRecord.extra));
		String header_st = combine(header,"\t");
		no_cols = header_st.split("\t").length;
		vcf_out.println(header_st);
			
		VarRecord prevVar = null;
		VarRecord fVar = nextRecord(bf, prevVar, myChrom);
		if(noBAM){
			while(fVar!=null){
	//		System.err.println("adding "+fVar_T);
				fVar.setMatchedTumour(varList_tumour.get(fVar.pos));
				varList.add(fVar);
				prevVar = fVar;
				fVar = nextRecord(bf, prevVar, myChrom);
			}
			
			
			clearUpTo(vcf_out,varList,refSeq.length()+10, myChrom , refSeq); // clear up varlist up to 10,000 bases before current
			if(varList.size()>0) throw new RuntimeException("did not clear all");
			vcf_out.close();
			if(ls!=null) ls.close();
			return;
		}else{
			fVar.setMatchedTumour(varList_tumour.get(fVar.pos));
			varList.add(fVar);
			prevVar = fVar;
		}
		SAMRecordIterator samIter = null;
		SamReader samReader = null;
		SAMTextWriter samWriter = null;
		int myChromIndex= -1;
		
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			 samReader = SamReaderFactory.makeDefault().open(new File(inFile));						
	
	
			SAMFileHeader samHeader = samReader.getFileHeader();		
			 samWriter = null;
			if(writeSAM){
			//	samWriter = outFile.equals("-")?	(new SAMTextWriter(System.out))	:(new SAMTextWriter(new File(outFile)));
				samWriter =new SAMTextWriter(new File(outDir, myChrom+".sam"));
				samWriter.setSortOrder(SortOrder.unsorted, false);		
				samWriter.writeHeader( samHeader.getTextHeader());	
			}
			///////////////////////////////////////////////////////////
	
			myChromIndex = samHeader.getSequenceIndex(myChrom);
			if(myChromIndex < 0){
				if(writeSAM) samWriter.close();
				samReader.close();
				bf.close();
				throw new RuntimeException("Chrom " + myChrom + " not found in bam file!!");
			}
			 samIter = samReader.query(refSeq.getName(),0,0,false);
				LOG.info(" " + samIter.hasNext());

		
		
		LOG.info("Chrom = " + myChrom + " RefName = " + refSeq.getName() + " chromIndex = " + myChromIndex);
		long time0=System.currentTimeMillis();
		int len = refSeq.length();
		long totalDiff=0;
		long totalDiffCount=0;
		int variantReads=0;
		int[] counts_shared = new int[sizeThreshold+1];
		int[] counts_unique = new int[sizeThreshold+1];
		int[] counts_filtered = new int[sizeThreshold+1];
		int[] counts = new int[sizeThreshold+1];
		
		
		
		for (int count=0; samIter.hasNext(); count++){
			SAMRecord sam = samIter.next();
		
			if (sam.getReadUnmappedFlag())
				continue;
			if(count % 1000000==0){
				
				long time1 = System.currentTimeMillis();
				double diff = time1 - time0;
				double time = ((double)diff)/(double)10000;
				time0 = time1;
				double avgLen = ((double)totalDiff/(double)totalDiffCount);
		//		System.err.println("processed "+count+" reads. Currently at "+sam.getAlignmentStart()/1e6+ " of "+len/1e6+" "+time+ "millis per read  "+diff/1000.0+"secs per million");
		//		System.err.println("percentage progress: "+ (double)sam.getAlignmentStart()/(double)len);
			//	System.err.println("total variant supporting reads "+variantReads+" of "+count);
			//	System.err.println("Unmatched: "+unMatchedCount+" matched: "+totalDiffCount+" frag length "+avgLen);
				 totalDiff=0;
				 totalDiffCount=0;
				 unMatchedCount=0;
			}
			boolean clear = count % 10000==0; // clear every 10000
			int samRefIndex = sam.getReferenceIndex(); 
			if (samRefIndex < myChromIndex){
			//	System.err.println("continue");
				continue;
			}

			if (samRefIndex > myChromIndex){
				//System.err.println("break");
				break;//while
			}
			
			
			if(clear && partial) clearUpTo(vcf_out,varList,sam.getAlignmentStart()-maxIns, myChrom , refSeq); // clear up varlist up to 10,000 bases before current
			
			String readName = sam.getReadName();	
			int insertSize = Math.abs(sam.getInferredInsertSize());
			if(insertSize==0 || insertSize > sizeThreshold){
				AlternativeAllelesCmd.unMatchedCount+=1;
				//continue;
			}else{
				counts[insertSize]++;
				totalDiff +=insertSize;
				totalDiffCount+=1;
			}
			Sequence readSeq = new Sequence(Alphabet.DNA(), sam.getReadString(), sam.getReadName());
			boolean support = false;
			boolean shared=false;
			boolean unique = false;
			
			int readPos = 0;//start from 0					
			int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index
			byte[] baseQ = sam.getBaseQualities();
			int mapQ = sam.getMappingQuality();
			int[] array = new int[2];
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
									variantReads+=1;
									//System.err.println("found read with variant "+refPos+" "+sam.getReadName());
									support = true;
									if(var.tumor!=null){
										shared=true;
									}else{
										unique= true;
									}
									var.add(insertSize);
								//	pairedPos.addSNP(var);
									array[0] = baseQ[readPos+i]; array[1]= mapQ; //array[2] = st; array[3] = end;
									var.add(array);
									
									break;//for									
								}

								currentVarPos = var.pos;								
								if (currentVarPos > refPos + i)
									break;//for
							}							
							if (support)
								break;//for i

							while (currentVarPos < refPos + i && hasVar){
								VarRecord var = nextRecord(bf, prevVar, myChrom);
								prevVar = var;
								if (var == null){
									hasVar = false;
									break;
								}
								var.setMatchedTumour(varList_tumour.get(var.pos));
							//	System.err.println("adding "+var.toString());
								varList.add(var);

								if (var.pos == refPos + i && var.base == readBase){
									//yay
								//	System.err.println("found read with variant "+refPos+" "+sam.getReadName());
									variantReads+=1;
									support = true;
									if(var.tumor!=null){
										shared=true;
									}else{
										unique= true;
									}
									var.add(insertSize);
//									pairedPos.addSNP(var);
									array[0] = baseQ[readPos+i]; array[1]= mapQ; //array[2] = st; array[3] = end;
									var.add(array);
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
			if(support && insertSize < sizeThreshold){
				if(shared) counts_shared[insertSize]++;
				if(unique)counts_unique[insertSize]++;
				counts_filtered[insertSize]++;
			}
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
		
		if(writeSAM && somaticSet!=null){
			samIter = samReader.query(refSeq.getName(),0,0,false);
			while (samIter.hasNext()){			
				SAMRecord sam = samIter.next();
				String readName = sam.getReadName();
				if (somaticSet.contains(readName)){
					samWriter.writeAlignment(sam);
				}
			}
			
			
			samWriter.close();
			
		}
		samReader.close();
		bf.close();
		clearUpTo(vcf_out,varList,refSeq.length()+10, myChrom , refSeq); // clear up varlist up to 10,000 bases before current
		if(varList.size()>0) throw new RuntimeException("did not clear all");
		vcf_out.close();
		if(ls!=null) ls.close();
		
		
		PrintWriter cnts_out = new PrintWriter(
				new GZIPOutputStream(new FileOutputStream(new File(outDir, vcfFile.getName()+"."+myChrom+".out.cnts.gz"))));
		cnts_out.println("len,filtered,shared,unique,all");
		for(int i=0; i<counts.length; i++){
			cnts_out.println(i+","+counts_filtered[i]+","+counts_shared[i]+","+counts_unique[i]+","+counts[i]);
		}
		cnts_out.close();
		
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
	 static String combine(String[] header2, String join) {
		 return combine(Arrays.asList(header2),join);
	 }
	 static String combine(List<String> header2, String join) {
		StringBuffer sb = new StringBuffer(header2.get(0));
		for(int i=1; i<header2.size(); i++){
			sb.append(join);
			sb.append(header2.get(i));
		}
		return sb.toString();
	}


}
