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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import japsa.seq.Alphabet;
import japsa.seq.FastqSequence;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.HTSUtilities.IdentityProfile;
import japsa.util.deploy.Deployable;


/**
 * @author lachlancoin
 *
 */
@Deployable(
	scriptName = "jsa.hts.repeatDetection",
	scriptDesc = "Detecting repeats in long read sequencing data")
public class RepeatDetectionCmd extends CommandLine{
	static String src_tag = "SC";
	static String pool_tag = "PT";
	
	
	
	
	
	private static final class CombinedIterator implements Iterator<SAMRecord> {
		private final SAMRecordIterator[] samIters;
		int current_sam_index =0;
		int currentIndex =0; //relates to chromosome
		private final SAMRecord[] currentVals;
		private final boolean[] returned;
		private final int[] cnts ;
		int max;
		Collection<String>[] readList ;
		Set<Integer> chrs;
		private CombinedIterator(SAMRecordIterator[] samIters, int max, Collection<String>[]readList, Set<Integer> chrs) {
			this.samIters = samIters;
			this.readList = readList;
			currentVals = new SAMRecord[samIters.length];
			returned = new boolean [samIters.length];
			this.max = max;
			cnts = new int[samIters.length];
			this.chrs = chrs;
		}

		@Override
		public boolean hasNext() {
			for(int i=0; i<samIters.length; i++){
				if(samIters[i].hasNext() && cnts[i]<max) return true;
			}
			return false;
		}
		int pool_ind=-1;
		private SAMRecord next(int i){
			SAMRecord sr;
			if(currentVals[i]!=null && !returned[i]){
				sr = currentVals[i];
			}else{
				if(cnts[i]<max){
					
					//sr  = samIters[i].next();
				/*	while(!(
							sr==null || 
							(readList==null || readList.contains(sr.getReadName())) ||
							(chrs==null || chrs.contains(sr.getReferenceIndex()))
							)
							)*/
					inner: while(true)
					{
						sr  = samIters[i].next();
						if(sr==null) break inner;
						pool_ind =-1;
						if(readList!=null){
							inner1: for(int i2=0; i2<readList.length; i2++){
								if(readList[i2].contains(sr.getReadName())){
									pool_ind = i2;
									break inner1;
								}
							}
						}
						if(readList!=null && pool_ind>=0) break inner;
						if(readList==null &&  chrs!=null && chrs.contains(sr.getReferenceIndex())) break inner;
						if(readList==null  && chrs==null) break inner;
					}
					if(sr!=null) {
						sr.setAttribute(pool_tag, pool_ind);
						sr.setAttribute(src_tag, i);
					}
					currentVals[i] = sr;
					returned[i] = false; 
				}else{
					sr = null;
				}
			}
			return sr;
		}

		@Override
		public SAMRecord next() {
			//int curr_index = current_sam_index;
			SAMRecord sr = next(current_sam_index);
			
			if(sr==null || sr.getReferenceIndex()>currentIndex ){
				int[] ref_inds = new int[samIters.length];
				int min_ind =-1;
				int minv = Integer.MAX_VALUE;
				for(int i=0; i<samIters.length; i++){
					SAMRecord sr_i = next(i);
					if(sr_i!=null) {
						ref_inds[i] = sr_i.getReferenceIndex(); //note this only advances if not returned;
						if(ref_inds[i]<minv){
							min_ind = i;
							minv = ref_inds[i];
						}
					}
				}
				current_sam_index = min_ind;
			}
			
			if(current_sam_index<0) return null;
			sr = this.currentVals[current_sam_index];
			if(sr!=null) this.currentIndex = sr.getReferenceIndex();
			cnts[current_sam_index]++;
			returned[current_sam_index] = true; 
			return sr;
			
		}
	}
	
	
//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

	 public static int insThresh = 200;
	 public static int flankThresh = 200;
	 static final FastqWriterFactory factory = new FastqWriterFactory();
	public RepeatDetectionCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("bamFile", null,  "Name of bam file", true);
		addString("bedFile", null,  "Name of bed file", true);
		addString("readList", null,  "reads to include", false);
		addString("reference", null, "Name of reference genome",true);
		addString("pattern", null, "Pattern of read name, used for filtering");
		addInt("qual", 0, "Minimum quality required");
		addInt("insThresh", 200, "Insertion threshold");
		addInt("flankThresh", 200, "Flank threshold");
		addStdHelp();		
	} 



	public static void main(String [] args) throws IOException, InterruptedException{		 		
		CommandLine cmdLine = new RepeatDetectionCmd();		
		args = cmdLine.stdParseLine(args);		

		String reference = cmdLine.getStringVal("reference");		
		int qual = cmdLine.getIntVal("qual");
		String pattern = cmdLine.getStringVal("pattern");
		String bamFile = cmdLine.getStringVal("bamFile");
		String bedFile = cmdLine.getStringVal("bedFile");
		String readL= cmdLine.getStringVal("readList");
		
		File resDir = new File("./results");
		resDir.mkdir();
		RepeatDetectionCmd.insThresh = cmdLine.getIntVal("insThresh");
		RepeatDetectionCmd.flankThresh = cmdLine.getIntVal("flankThresh");
		String[] bamFiles_ = bamFile.split(":");
		if(bamFile.equals("all") || bamFile.equals(".")){
			bamFiles_ = (new File("./")).list(new FilenameFilter(){

				@Override
				public boolean accept(File dir, String name) {
					return name.endsWith(".bam");
				}
				
			});
		}
		
		errorAnalysis(bamFiles_, reference, pattern, qual, resDir,readL==null ? null : readL.split(":"));		


		//paramEst(bamFile, reference, qual);
	}

	static FastqRecord makeRecord(Sequence readSeq, String baseQ, String chr, String suffix,  int st_read, int end_read, int st_ref, int end_ref){
		int len = end_read - st_read;
		return  new FastqRecord(
				readSeq.getName()+suffix+" "+chr+" "+len+" "+st_read+","+end_read+" "+st_ref+","+end_ref,
				new String(readSeq.subSequence(st_read, end_read).charSequence()),"",
				baseQ.substring(st_read, end_read)
				);
	}


	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 */
	static void errorAnalysis(String[] bamFiles_, String refFile, String pattern, int qual, File resDir, String[] readList) throws IOException{	
		FastqWriter fastq = factory.newWriter( new File(resDir,"insertion.fastq"));
		FastqWriter fastq_flank = factory.newWriter( new File(resDir,"flank.fastq"));
		FastqWriter fastq_none = factory.newWriter( new File(resDir,"no_insertion.fastq"));
		List<Insertions> insertions = new ArrayList<Insertions>();
		//Read the reference genome
		ArrayList<Sequence> genomes = SequenceReader.readAll(refFile, Alphabet.DNA());

	Integer max_reads = Integer.MAX_VALUE;
		Set<String>chrToInclude = null;
		int len = bamFiles_.length;
		final SAMRecordIterator[] samIters = new SAMRecordIterator[len];
		SamReader[] samReaders = new SamReader[len];
		for (int ii = 0; ii < len; ii++) {
			int source_index = ii;
			
			String bamFile = bamFiles_[ii];
			File bam = new File( bamFile);
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			if ("-".equals(bamFile))
				samReaders[ii] = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				samReaders[ii] = SamReaderFactory.makeDefault().open(bam);

			 samIters[ii] = samReaders[ii].iterator();
			// Read the reference genome
		}
		Set<Integer> chrom_indices_to_include = null;
		if(chrToInclude!=null){
			chrom_indices_to_include= new HashSet<Integer>();
			for(int i=0; i<genomes.size(); i++){
				if(chrToInclude.contains(genomes.get(i).getName())) chrom_indices_to_include.add(i);
			}
		}
		Collection[] reads = getReads(readList);
		Iterator<SAMRecord> samIter;
		if(reads==null && chrom_indices_to_include==null && samIters.length==1){
			samIter = samIters[0];
		}else{
			samIter = 		new CombinedIterator(samIters, max_reads,reads, chrom_indices_to_include);
		}
		
		
		
		
		//get the first chrom
		int currentIndex = 0;
		Sequence chr = genomes.get(currentIndex);

	
		long    totBaseIns = 0,
			totBaseDel = 0,
			totNumIns = 0,
			totNumDel = 0,
			totMisMatch = 0,
			totMatch = 0,
			totClipped = 0;

		long totReadBase = 0, totRefBase = 0;
		int  numReads = 0;

		int numNotAligned = 0;

		//String log = "###Read_name\tRead_length\tReference_length\tInsertions\tDeletions\tMismatches\n";
		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();

			if (pattern != null && (!sam.getReadName().contains(pattern)))
				continue;

			//make the read seq			
			Sequence readSeq = new Sequence(Alphabet.DNA(), sam.getReadString(), sam.getReadName());
			if (readSeq.length() <= 1){
				//LOG.warn(sam.getReadName() +" ignored");
				//TODO: This might be secondary alignment, need to do something about it
				continue;
			}			

			numReads ++;


			if (sam.getReadUnmappedFlag()){
				numNotAligned ++;
				continue;
			}

			int flag = sam.getFlags();


			if (sam.getMappingQuality() < qual) {
				numNotAligned ++;
				continue;
			}



			//int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index
			int refIndex = sam.getReferenceIndex();

			//if move to another chrom, get that chrom
			if (refIndex != currentIndex){
				currentIndex = refIndex;
				chr = genomes.get(currentIndex);
			}
			
			insertions.clear();
			japsa.util.HTSUtilities.IdentityProfile profile = 
				identity(chr, readSeq, sam, insertions);	
			
			String desc = chr.getName()+","+sam.getAlignmentStart()+","+sam.getAlignmentEnd()+","+readSeq.length();
			//sb.append(" ");
		if(insertions.size()>0 || true ) 	System.err.println(readSeq.getName()+" "+desc+" "+(readSeq.length()-(sam.getAlignmentEnd()-sam.getAlignmentStart()))+"  "+ insertions.size());
			String baseQ = sam.getBaseQualityString();
			if(insertions.size()>0){
				//Collections.sort(insertions);
				//StringBuffer sb = new StringBuffer(desc);
				//int prevp = 0;
				int readStart =0;
				int refStart = sam.getReferencePositionAtReadPosition(readStart);
				int readEnd=-1;
				int refEnd = -1;
				for(int i=0; i<insertions.size();i++){
					Insertions ins = insertions.get(i);
					readEnd = ins.readStart;
					refEnd = ins.refStart;
					
					if(readEnd-readStart > flankThresh){
						FastqRecord left = makeRecord(readSeq,baseQ, desc,".F"+i , readStart, readEnd, refStart, refEnd);
						fastq_flank.write(left);
					}
					FastqRecord repeat = makeRecord(readSeq, baseQ, desc, ".R"+i, ins.readStart, ins.readEnd(), ins.refStart, ins.refStart);
					fastq.write(repeat);
					refStart = ins.refStart;  //next refStart
					readStart = ins.readEnd(); //next readStart
				}
				readEnd = readSeq.length();
				refEnd = sam.getReferencePositionAtReadPosition(readEnd);
				int i = insertions.size();
				if(readEnd-readStart > flankThresh){
					FastqRecord left = makeRecord(readSeq,baseQ, desc,".F"+i , readStart, readEnd, refStart, refEnd);
					fastq_flank.write(left);
				}
			}else{
				if(readSeq.length()>flankThresh){
					fastq_none.write(new FastqRecord(readSeq.getName()+ " "+desc, new String(readSeq.charSequence()), "", baseQ));
				}
			}
			//log+=sam.getReadName() + "\t" + profile.readBase + "\t" + profile.refBase + "\t" + profile.baseIns + "\t" + profile.baseDel + "\t" + profile.mismatch+ "\n";

			totBaseIns  += profile.baseIns;
			totBaseDel  += profile.baseDel;
			totNumIns   += profile.numIns;
			totNumDel   += profile.numDel;			
			totMisMatch += profile.mismatch;
			totMatch    += profile.match;

			totReadBase += profile.readBase;
			totRefBase  += profile.refBase;			
			totClipped += profile.readClipped;
			//numReadsConsidered ++;
			
		}		
		for(int i=0; i<samReaders.length; i++){
		samReaders[i].close();
		}

		System.out.println("========================= TOTAL ============================");

		//Done

		fastq.close();fastq_flank.close();fastq_none.close();
		System.out.println("Deletion " + totBaseDel + " " + totNumDel +" " + totBaseDel*1.0/totRefBase);
		System.out.println("Insertion " + totBaseIns + " " + totNumIns+" " + totBaseIns*1.0/totRefBase);
		System.out.println("MisMatch " + totMisMatch +" " + totMisMatch*1.0/totRefBase);
		System.out.println("Match " + totMatch);
		System.out.println("Clipped " + totClipped);
		


		System.out.println("ReadBase " + totReadBase);
		System.out.println("ReferenceBase " + totRefBase);	

		double totState0 = totMatch + totMisMatch;//
		double probDel = (totNumDel + 1.0) / (totState0 + 3.0);
		double probIns = (totNumIns + 1.0) / (totState0 + 3.0);

		//double probMatch = 1.0 - probDel - probIns;
		double probCopy =  (totMatch + 1.0) / (totState0 + 2.0);
		double probChange = (totMisMatch + 1.0) / (totState0 + 2.0);

		double probDE = (1.0 + totBaseDel - totNumDel) / (2.0 +totBaseDel);
		double probIE = (1.0 + totBaseIns - totNumIns) / (2.0 +totBaseIns);		

		System.out.printf("Identity %f %f %f %f\n",1.0 *totMatch/(totMatch + totMisMatch + totBaseDel +totBaseIns),
			1.0 *totMisMatch/(totMatch + totMisMatch + totBaseDel +totBaseIns),
			1.0 *totBaseIns/(totMatch + totMisMatch + totBaseDel +totBaseIns),
			1.0 *totBaseDel/(totMatch + totMisMatch + totBaseDel +totBaseIns ));

		System.out.printf("Probs %f %f %f %f %f %f\n",probCopy, probChange, probIns, probDel, probIE, probDE);		

		System.out.println("========================= SUMMARY============================");

		System.out.printf("Total reads      :  %d\n", numReads);
		System.out.printf("Unaligned reads  :  %d\n", numNotAligned);

		System.out.printf("Deletion rate    : %.4f\n",totBaseDel*1.0/totRefBase);
		System.out.printf("Insertion rate   : %.4f\n",totBaseIns*1.0/totRefBase);
		System.out.printf("Mismatch rate    : %.4f\n",totMisMatch*1.0/totRefBase);
		System.out.printf("Identity rate    : %.4f\n",totMatch*1.0/totRefBase);
		System.out.println("=============================================================");
		System.out.println("=============================================================");

		//System.out.println(log);
	}
 static Collection[] getReads(String[] readList) throws IOException {
	Collection[] reads = null;
	if(readList!=null && readList.length>0){
	reads = new Collection[readList.length];
	 for(int i=0; i<reads.length; i++){
		reads[i] = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(readList[i])));
		String st;
		while((st = br.readLine())!=null){
			String st_ = st.split("\\s+")[0];
		//	System.err.println(st_);
			reads[i].add(st_);
		}
		br.close();
	 }
	}
	return reads;
	}
static class Insertions implements Comparable{
	Integer readStart=0;
	Integer refStart =0;
	Integer length =0;
	Insertions(int refStart, int readStart, int length){
		this.readStart = readStart;
		this.refStart = refStart;
		this.length = length;
	}
	public Integer readEnd() {
		// TODO Auto-generated method stub
		return readStart + length;
	}
	@Override
	public int compareTo(Object o) {
		// TODO Auto-generated method stub
		return length.compareTo(((Insertions)o).length);
	}
	public String toString(){
		return refStart+","+readStart+","+length;
	}
}
	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public static IdentityProfile identity(Sequence refSeq, Sequence readSeq,  SAMRecord sam, List<Insertions> insertions){
		IdentityProfile profile = new IdentityProfile();

		int readPos = 0;//start from 0					
		int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index				

		profile.readClipped = 0;
		profile.refClipped = sam.getAlignmentStart() + refSeq.length() - sam.getAlignmentEnd();
		profile.baseDel = 0;
		profile.baseIns = 0;
		profile.numDel = 0;
		profile.numIns = 0;
		profile.match = 0;
		profile.mismatch = 0;
		profile.refBase = 0;
		profile.readBase = 0;//the number of bases from ref and read

		for (final CigarElement e : sam.getCigar().getCigarElements()) {
			final int  length = e.getLength();
			switch (e.getOperator()) {
			case H :
				//nothing todo
				profile.readClipped += length;
				break; // ignore hard clips
			case P : 
				profile.readClipped += length;
				//pad is a kind of hard clipped ?? 					
				break; // ignore pads	                
			case S :
				//advance on the reference
				profile.readClipped += length;
				readPos += length;
				break; // soft clip read bases	                	
			case N : 
				refPos += length; 
				profile.refClipped += length;
				break;  // reference skip

			case D ://deletion      	
				refPos += length;
				profile.refBase += length;

				profile.baseDel += length;
				profile.numDel ++;
				break; 	

			case I :	    
				if(length>insThresh){
					insertions.add(new Insertions(refPos, readPos, length));
				//	System.err.println(length);
				}
				readPos += length;
				profile.readBase += length;
				
				
				profile.baseIns += length;
				profile.numIns ++;
				break;
			case M :
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++){
					if (refSeq.getBase(refPos + i) == readSeq.getBase(readPos + i))
						profile.match ++;
					else
						profile.mismatch ++;
				}
				profile.readBase += length;
				profile.refBase += length;

				readPos += length;
				refPos  += length;
				break;

			case EQ :
				readPos += length;
				refPos  += length;

				profile.readBase += length;
				profile.refBase += length;
				profile.match += length;
				break;

			case X :
				readPos += length;
				refPos  += length;

				profile.readBase += length;
				profile.refBase += length;

				profile.mismatch += length;
				break;
			default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}//case
		}//for			

		return profile;

	}
	
	
}

/*RST*
----------------------------------------------------------
*jsa.hts.errorAnalysis*: Error analysis of sequencing data
----------------------------------------------------------

*jsa.hts.errorAnalysis* assesses the error profile of sequencing data by getting the numbers
of errors (mismatches, indels etc) from a bam file. Obviously, it does not distinguish
sequencing errors from mutations, and hence consider mutations as errors. It is best to use
with the bam file from aligning sequencing reads to a reliable assembly of the sample.

<usage>

*RST*/
