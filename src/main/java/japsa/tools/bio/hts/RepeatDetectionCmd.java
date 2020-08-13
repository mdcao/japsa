/*****************************************************************************
 * Copyright (c) Lachlan Coin, University of Meblourne, All rights reserved.         *
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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqEncoder;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.SequenceUtil;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.tools.seq.SequenceUtils;
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
	
	
	static class Repeat{
		int st; int end; String label; String seq; int period; float nrep;
		public boolean equals(Object obj){
			return pos0 ==((Repeat)obj).pos0;
		}
		public int hashCode(){
			return pos0;
		}
		int pos0; // position in compressed reference
		Repeat(String[] str, int key){
			this.st =Integer.parseInt(str[3]);
			this.end = Integer.parseInt(str[4]);
			this.label = str[7];
			this.seq = str[8];
			this.period = Integer.parseInt(str[6]);
			this.nrep =1.0f +((float)end - (float) st)/(float) period;
		this.pos0 = key;
		}
		Repeat(){
			st = -1;
			end=-1;
			label="NA";
			seq = "NA";
			period=-1;
			nrep=Float.NaN;
		}
		int offset =0;
		public void setOffset(int i) {
			offset  = i;
		}
		public int pos(){
			return end+offset;
		}
		
		
	}
	static Repeat null_rep = new Repeat();
	
	static class Maps{
		SortedMap<Integer, Repeat> map = new TreeMap<Integer,Repeat>();
		List<String> header;
		 Maps(File bed) throws FileNotFoundException, IOException{
			 BufferedReader br = 	new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(bed))));
			header = Arrays.asList( br.readLine().split("\t"));
			 String st = "";
			 while((st = br.readLine())!=null){
				 String[] str = st.split("\t");
				 Integer key = Integer.parseInt(str[2]);

				 Repeat rep = new Repeat( str, key);
				// if(first) this.chr = str[0];
				 map.put(key, rep);//Integer.parseInt(str[4]));
			 }
			 br.close();
		 }
		 public int getPos(int pos){
			 if(map.containsKey(pos)) return map.get(pos).st;
			SortedMap <Integer, Repeat>hm = map.headMap(pos+1);
			if(hm.size()==0) return pos;
			Integer lastKey = hm.lastKey();
			Repeat rep = map.get(lastKey);
			return rep.end + (pos-lastKey);
		 }
		public Repeat find(Insertions ins) {
			// TODO Auto-generated method stub
			if(map.containsKey(ins.refStart)) return map.get(ins.refStart);
			SortedMap <Integer, Repeat>hm  =  map.tailMap(ins.refStart);
			Integer key1 = hm.size()==0 ? null : hm.firstKey();
			SortedMap <Integer, Repeat>hm1  =  map.headMap(ins.refStart);
			Integer key2 = hm1.size()==0 ? null : hm1.lastKey();
		//	ins.overlap(key1);
			int diff1 = ins.overlap(key1);
			int diff2 = ins.overlap(key2);
			if(Math.max(diff1, diff2) >=-50){
				Integer v = diff1 >=diff2 ? key1 : key2;
				//System.err.println(ins.refStart+" "+(v-ins.refStart));
				return map.get(v);
			}else{
				System.err.println("not found "+ins.refStart+" "+diff1+','+diff2);
				return null;
			}
		}
		
	 }
	
	
	
	
	
	
//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

	 public static int insThresh = 200;
	 public static int insThreshKnown = 200;
	 public static int insThreshUnknown = 500;
	 public static int flankThresh = 200;
	 public static boolean extractInsertion = false;
	 static final FastqWriterFactory factory = new FastqWriterFactory();
	public RepeatDetectionCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addBoolean("extractInsertion", false, "whether to include only insertion seq in fastq");
		addString("bamFile", null,  "Name of bam file", false);
		addString("fastq", null,  "Name of fastq file", false);
		addString("resDir", null,  "results dir", false);
		addString("chromsDir", null,  "Dir with compressed chroms", true);
		addString("readList", null,  "reads to include", false);
		addString("chromIndices", null,  "chrom indices to include", false);
		addString("reference", null, "Name of reference genome",true);
		addString("pattern", null, "Pattern of read name, used for filtering");
		addInt("qual", 0, "Minimum quality required");
		addInt("flank_req", 20, "Minimum flank required on the read either side of the insertion");
		addInt("bin", 1000, "For grouping of insertions");
		addString("chroms", null, "chrom indices to include",false);

	//	addBoolean("stageOne", true, "first stage");
		addString("insThresh", "200:200", "Insertion threshold known:Unknown");
		addInt("flankThresh", 200, "Flank threshold");
		addString("mm2_path", "/sw/minimap2/current/minimap2",  "minimap2 path", false);
		addString("mm2_memory", (Runtime.getRuntime().maxMemory()-1000000000)+"",  "minimap2 memory", false);
		addInt("mm2_threads", 4, "threads for mm2", false);
		addStdHelp();		
	} 

	static int mm2_threads=4;
	static String mm2_mem = "5g";
	static String mm2Preset="map-ont";
	static String mm2_path="/home/lachlan/github/minimap2/minimap2";
static int flank_req;

	public static void main(String [] args) throws IOException, InterruptedException{		 		
		CommandLine cmdLine = new RepeatDetectionCmd();		
		args = cmdLine.stdParseLine(args);		

		String reference = cmdLine.getStringVal("reference");		
		int qual = cmdLine.getIntVal("qual");
		String pattern = cmdLine.getStringVal("pattern");
		String bamFile = cmdLine.getStringVal("bamFile");
		String chromsDir = cmdLine.getStringVal("chromsDir");
		String readL= cmdLine.getStringVal("readList");
		flank_req = cmdLine.getIntVal("flank_req");
		mm2_threads = cmdLine.getIntVal("mm2_threads");
		mm2_mem = cmdLine.getStringVal("mm2_mem");
		mm2_path = cmdLine.getStringVal("mm2_path");
		bin = cmdLine.getIntVal("bin");
		String resD = (cmdLine.getStringVal("resDir"));
		File resDir = resD==null ? new File("./"+((new File(chromsDir))).getName()) : new File(resD);
		resDir.mkdir();
		String[] insT = cmdLine.getStringVal("insThresh").split(":");
		RepeatDetectionCmd.insThreshKnown = Integer.parseInt(insT[0]);
		RepeatDetectionCmd.insThreshUnknown = Integer.parseInt(insT[1]);
		RepeatDetectionCmd.insThresh = Math.min(insThreshKnown, insThreshUnknown);
		RepeatDetectionCmd.flankThresh = cmdLine.getIntVal("flankThresh");
		RepeatDetectionCmd.extractInsertion = cmdLine.getBooleanVal("extractInsertion");
		
	//	RepeatDetectionCmd.chromsDir = chromsDir;
		
		String[] bamFiles_ = bamFile.split(":");
		if(bamFile.equals("all") || bamFile.equals(".")){
			bamFiles_ = (new File("./")).list(new FilenameFilter(){

				@Override
				public boolean accept(File dir, String name) {
					return name.endsWith(".bam");
				}
				
			});
		}
		File insertionsDir = new File(resDir,"insertions");
		File insertionsDir1 = new File(resDir,"insertionsNoOverlap");
		delete(insertionsDir); delete(insertionsDir1);
	
		analysis(bamFiles_, new File(reference), pattern, cmdLine.getStringVal("chromIndices"), qual, resDir,new File(chromsDir), readL==null ? null : readL.split(":"));		


		//paramEst(bamFile, reference, qual);
	}

	private static void delete(File insertionsDir) {
		if(insertionsDir.exists()){
			File[] f = insertionsDir.listFiles();
			for(int i=0; i<f.length; i++) {
				f[i].delete();
			}
		}else{
			insertionsDir.mkdir();
		}
		
	}
/** null value of negStrand indicates do not look at strand - keep same as original */
	static FastqRecord makeRecord(Sequence readSeq, String baseQ, byte[] quals,  String suffix,  int st_read, int end_read, String desc, boolean negStrand){
	
		String sequence = new String(readSeq.subSequence(st_read, end_read).charSequence());
		String baseQL = baseQ.substring(st_read, end_read);
		String avgQ="";
		
			String	ch=negStrand ?  "_rev" : "";
			 double sum =0;
				for(int i=st_read; i<end_read; i++){
					sum+=quals[i];
				}
				double length = end_read-st_read;
				avgQ = " "+String.format("%5.3g", sum/length).trim();
			
		
		return  new FastqRecord(
				readSeq.getName()+ch+suffix+" "+st_read+"-"+end_read+" "+desc+avgQ,
				negStrand ? SequenceUtil.reverseComplement(sequence): sequence,
				"",
				negStrand ?(new StringBuilder(baseQL)).reverse().toString(): baseQL
				);
	}

	

	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 */
	static void analysis(String[] bamFiles_, File refFile, String pattern, String chroms, 
			int qual, File resDir, File chromsDir, String[] readList) throws IOException{	
		Integer max_reads = Integer.MAX_VALUE;
	//	Set<String>chrToInclude = null;
		int len = bamFiles_.length;
			final Iterator<SAMRecord>[] samIters = new Iterator[len];
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
		}
		
		
		Iterator<SAMRecord>  samIter = SequenceUtils.getCombined(samIters, getReads(readList), max_reads, chroms, true);
		
			processStageOne(samIter, samIters.length, resDir, chromsDir, pattern, qual, refFile.getName());
		for(int i=0; i<samReaders.length; i++){
			samReaders[i].close();
		}
		//for(int i=0; i<samReaders.length; i++){
		//	samReaders[i].close();
		//	}
		//get the first chrom
		
		//System.out.println(log);
	}
	
	//${mm2} -t ${threads} -ax map-ont ${reference} ${read_dir}/${file}.fastq
	
	 private static void processStageOne(Iterator<SAMRecord> samIter, int num_source, File resDir, File chromDir, String pattern, int qual, String reference) throws IOException {
		 FastqW fastq  =null;
		 if(!chromDir.exists()) throw new RuntimeException(chromDir.getAbsolutePath()+" does not exist");
		 SamReader[] outr = new SamReader[num_source];
		 int currentIndex = -1;
			//Sequence chr = genomes.get(currentIndex);
			int  numReads = 0;
			int numNotAligned = 0;

			while (samIter.hasNext()){
				SAMRecord sam = samIter.next();
			//	System.err.println(sam.getReferenceIndex()+" "+sam.getReferenceName());
				if(sam==null) break;
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
				String chrom = sam.getReferenceName();
				//if move to another chrom, get that chrom
				if (refIndex != currentIndex){
					currentIndex = refIndex;
					//sam.get
					//chr = genomes.get(currentIndex);
					if(fastq!=null){
							fastq.close();
						
						
					}
					String out = reference.replace(".fasta", "").replace(".fa", "").replace(".gz", "");
					String outf = chromDir.getAbsolutePath()+"/"+chrom+"."+out+".no_repeats.fa.gz";
					File bed = new File(chromDir.getAbsolutePath()+"/"+chrom+"."+out+".bed.gz");
					String  mm2Index =outf+".mmi";
					File seqF = new File(outf);
					if(seqF.exists()){
						ArrayList<Sequence> genomes = SequenceReader.readAll(outf, Alphabet.DNA());
						Maps maps = new Maps(bed);
						//int chk = maps.getPos(212738260);
						//System.err.println(212738260+"->"+chk);
						
							fastq  =new FastqW (resDir, currentIndex,"", mm2Index, genomes, maps);
					
					}else{
						System.err.println("does not exist "+seqF.getAbsolutePath());
						fastq = null;
					}
					
				}
				int source = (Integer) sam.getAttribute(SequenceUtils.src_tag);
				if(fastq!=null){
					String baseQ = sam.getBaseQualityString();
					readSeq.setName(readSeq.getName()+"."+source);
				//	String strand = sam.getReadNegativeStrandFlag() ? "-" : "+";
			//		String desc = "";//chr.getName();//+","+sam.getAlignmentStart()+","+sam.getAlignmentEnd()+","+readSeq.length();
					
					FastqRecord repeat =  new FastqRecord(readSeq.getName(),	readSeq.toString(),	"",baseQ);
													
													
														
						//	makeRecord(readSeq, baseQ,  sam.getBaseQualities(), "",0, sam.getReadLength(), "", sam.getReadNegativeStrandFlag());
					fastq.write(repeat);
				}
			}
			if(fastq!=null){
				fastq.close();
			}	
			System.out.println("========================= SUMMARY============================");
			System.out.printf("Total reads      :  %d\n", numReads);
			System.out.printf("Unaligned reads  :  %d\n", numNotAligned);
			System.out.println("=============================================================");

			
		}
	 
	 
	static ExecutorService executor = Executors.newFixedThreadPool(1);
	 static class FastqW{
			final FastqWriter fq;
			final  File inFile; 
			final String  mm2Index;
			final File resDir;
			final ArrayList<Sequence > genomes;
			final int chrom_index;
			
			final Maps maps;
			public FastqW(File resDir, int currentIndex, String prefix, String mm2Index, ArrayList<Sequence > genomes, Maps maps) throws FileNotFoundException, IOException {
				this.inFile = 	new File(resDir,currentIndex+prefix+".fastq");
				fq = factory.newWriter( inFile);
			//	fq = new FQWriter(new FileOutputStream(inFile));
				this.genomes = genomes;
				this.maps = maps;
				this.resDir = resDir;
				this.mm2Index = mm2Index;
				this.chrom_index = currentIndex;
			}
			public void close() throws IOException{
				fq.close();
				
				Iterator<SAMRecord> reader = SequenceUtils.getSAMIteratorFromFastq(inFile, mm2Index, mm2_path, mm2_threads, mm2Preset,  mm2_mem);
				
				Runnable run = new Runnable(){

					@Override
					public void run() {
							processStageTwo(genomes, chrom_index, reader, resDir,  null, 0, maps);
					}
					
				};
//				executor.execute(run);
				run.run();
			}
			
			public void write(FastqRecord repeat) {
				fq.write(repeat);
				
			}
		}
	 static float bin = 1000f;
	 
static class FQWriter implements FastqWriter{
	private  PrintStream writer;
	public FQWriter(String prefix, String suffix, boolean append) {
		try{
		writer = new PrintStream(new FileOutputStream(prefix+"."+suffix, append));
		 
		}catch(Exception exc){
			exc.printStackTrace();
			
			System.exit(0);
		}
		
	}

	public FQWriter(FileOutputStream os) {
		writer = new PrintStream(os);
	}

	@Override
	public void close() {
		writer.println();
		writer.close();
		
	}

	@Override
	public void write(FastqRecord rec) {
		// TODO Auto-generated method stub
		 FastqEncoder.write(writer, rec);
	}
	
}
static boolean writeNone =true;

 private static void processStageTwo(ArrayList<Sequence> genomes, int chrom_index, Iterator<SAMRecord> samIter, File resDir, String pattern, int qual, Maps maps) {
	File insertionsDir = new File(resDir, "insertions");
	File insertionsDir1 = new File(resDir, "insertionsNoOverlap");
	//SAMRecordIterator it = reader.iterator();
	//it.close();
//	reader.close();
	//CombinedIterator samIter = getCombined(samReaders, null,Integer.MAX_VALUE);
	String prefix = insertionsDir.getAbsolutePath()+"/"+chrom_index;
	String prefix1 = insertionsDir1.getAbsolutePath()+"/"+chrom_index;
	boolean append = true;
	 	
	 FastqWriter  fastq_none = writeNone ? new FQWriter(resDir.getAbsolutePath()+"/"+chrom_index,"no_insertion.fastq",append) : null;

		 List<Insertions>  insertions= new ArrayList<Insertions>();
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
		//int noInsertion=0;
		//int hasInsertion=0;
		int numNotAligned = 0;
		Map<Integer, Integer> num_insertions = new HashMap<Integer, Integer>();

		//String log = "###Read_name\tRead_length\tReference_length\tInsertions\tDeletions\tMismatches\n";
		outer1: while (samIter.hasNext()){
			SAMRecord sam = samIter.next();
			//boolean negStrand = sam.getReadNegativeStrandFlag();
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
				//	fastq  = factory.newWriter( new File(resDir,chr+".fastq"));
			}
			
			if(insertions !=null) insertions.clear();
			japsa.util.HTSUtilities.IdentityProfile profile = identity(chr, readSeq, sam, insertions);
			if(insertions.size()>1){
				for(int i=insertions.size()-1; i>=1; i--){
					if(insertions.get(i-1).merge(insertions.get(i),-20)){
						insertions.remove(i);
					}
				}
			}
			int st0 = sam.getAlignmentStart(); int end0 = sam.getAlignmentEnd(); 
			int st1 = maps.getPos(st0); int end1 = maps.getPos(end0);

			String desc = chr.getName()+","+st0+","+end0+","+(end0-st0)+" "+chr.getName()+","+st1+","+end1+","+(end1-st1)+" "+sam.getReadLength()+","+sam.getMappingQuality();
			//sb.append(" ");
		//	if(insertions!=null && 
			//		insertions.size()>0  ) 	System.err.println(readSeq.getName()+" "+desc+" "+(readSeq.length()-(sam.getAlignmentEnd()-sam.getAlignmentStart()))+"  "+ insertions.size());
			String baseQ = sam.getBaseQualityString();
				//FastqWriter  fastq_flank = null;
				int cnt=0;
				List<Repeat> reps = new ArrayList<Repeat>();
				List<Insertions> insertions1 = new ArrayList<Insertions>();
				
				int cnt_match=0;
				for(int i=0; i<insertions.size();i++){
					Insertions ins = insertions.get(i);
					 Repeat rep = maps.find(ins);
					 if(rep==null && ins.length >=insThreshUnknown){
						 insertions1.add(insertions.get(i));
						 reps.add(null);
					 }else if(rep!=null && ins.length >=insThreshKnown){
						 cnt_match++;
						 int ind =  reps.indexOf(rep);
						 if(ind>=0){
							if(ins.overlap(rep.pos0)>  insertions1.get(ind).overlap(rep.pos0)) insertions1.set(ind, ins);
						 }else{
							 insertions1.add(ins);
							 reps.add(rep);
						 }
					 }
				}
			//	putIfAbsent(insertions1.size(), 0);
				num_insertions.put(insertions1.size(), num_insertions.getOrDefault(insertions1.size(), 0)+1);
				
				if(insertions1.size()==0){
					System.err.println(readSeq.getName());
					if(writeNone){
						
						fastq_none.write(new FastqRecord(readSeq.getName()+ " "+desc, new String(readSeq.charSequence()), "", baseQ));
					}
					continue outer1;
				}
				Repeat rep; String nme = null;String seq = null;
				int readStart =0;
				outer: for(int i=0; i<insertions1.size();i++){
					Insertions ins = insertions1.get(i);
					rep = reps.get(i);
					int readEnd1 = ins.readEnd;
					
					if(i==insertions1.size()-1){
						readEnd1 = readSeq.length();
					}else{
						Insertions nxt = insertions1.get(i+1);
						Repeat nxt_rep = reps.get(i+1);
						double dist = Math.max(0,nxt.readStart - ins.readEnd);
						readEnd1 = readEnd1 + (int) Math.round(dist/2.0);
					}
					if(rep==null){
						nme = prefix1 + "."+round(ins.refStart);
						seq="";
					}else{
						nme = prefix+"."+rep.label;
						seq = rep.seq;
					}
					if(extractInsertion){
						if(rep==null){
							readStart = Math.max(readStart,ins.readStart-10);
							readEnd1 = Math.min(readEnd1, ins.readEnd+10);
						}else{ // include one repeat width either way
							readStart = Math.max(readStart,ins.readStart-rep.period);
							readEnd1 = Math.min(readEnd1, ins.readEnd+rep.period);
						}
					}
						int ins_pos = maps.getPos(ins.refStart);
					//	double len = ins.length;
						double nrep_ = rep==null ? Double.NaN : 1.0 + (double)ins.length/(double) rep.period;
						String nrep = rep==null ? "NA": String.format("%5.3g", nrep_).trim();
						String ratio = rep==null ? "NA": String.format("%5.3g", nrep_/(double) rep.nrep).trim();
						String period = rep==null ? "NA" : rep.period+"";
						int repStart = ins.readStart - readStart;
						int repEnd = ins.readEnd - readStart;
						String diff = rep==null ? "NA" : ""+(ins.refStart - rep.pos0);
						//System.err.println(sam.getReadNegativeStrandFlag());
						//System.err.println(readSeq.getName()+" "+nme);
						FastqRecord repeat =  makeRecord(readSeq, baseQ, sam.getBaseQualities(),".R"+i, readStart, readEnd1,
								repStart+"-"+repEnd+" "+ins.length+" "+ins_pos+" "+diff+" "+period
								+" "+nrep+" "+ratio+" "+seq+" "+ins.left_flank+","+ins.right_flank+ ","+sam.getMappingQuality(),
								sam.getReadNegativeStrandFlag());
						FastqWriter fastq  = new FQWriter(nme,"ins.fastq", append);
						fastq.write(repeat); fastq.close();
						readStart =readEnd1; //next readStart
				}

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
		
	//	samIter.close();
		System.out.println("========================= TOTAL "+chr.getName()+" ============================");

		//Done

		//fastq.close();fastq_flank.close();
		if(fastq_none!=null) fastq_none.close();
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

		System.out.println("========================= SUMMARY "+chr.getName()+"============================");

		System.out.printf("Total reads      :  %d\n", numReads);
		System.out.printf("Unaligned reads  :  %d\n", numNotAligned);
		System.out.println("Reads with insertion  "+ num_insertions);

		
		System.out.printf("Deletion rate    : %.4f\n",totBaseDel*1.0/totRefBase);
		System.out.printf("Insertion rate   : %.4f\n",totBaseIns*1.0/totRefBase);
		System.out.printf("Mismatch rate    : %.4f\n",totMisMatch*1.0/totRefBase);
		System.out.printf("Identity rate    : %.4f\n",totMatch*1.0/totRefBase);
		System.out.println("=============================================================");
		System.out.println("=============================================================");

		
	}





private static String round(int st0) {
	return ""+(int)Math.round(st0/bin);
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
	Integer refEnd =0;
	
	Integer length =0;
	Integer readEnd =0;
	Integer left_flank =0;
	Integer right_flank =0;
	//String key;
	Insertions(int refStart, int readStart, int length, int left_flank, int right_flank){
		this.readStart = readStart;
		this.refStart = refStart;
		this.refEnd = refStart; //need endd so we can merge adjacent insertions
		this.length = length;
		this.readEnd = readStart+length;
		this.left_flank = left_flank;
		this.right_flank =right_flank;
	//	int refEnd = refStart+length;
		
	}
	
	public int  overlap(Integer pos0){
		if(pos0==null) return Integer.MIN_VALUE;
		else return Math.min(refEnd-pos0, pos0-refStart );
	}
	public boolean merge(Insertions insertions, int threshold) {
		double overlap = Math.min(insertions.readEnd - readStart, readEnd  - insertions.readStart);
		if(overlap > threshold){
		this.readStart = Math.min(insertions.readStart, readStart);
		this.readEnd = Math.max(insertions.readEnd, readEnd);
		this.length = readEnd -readStart;
		this.refStart = Math.min(refStart, insertions.refStart);
		this.refEnd = Math.max(refEnd, insertions.refEnd);
		return true;
		}else{
			return false;
		}
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

static int flankReq=20;

/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public static IdentityProfile identity(Sequence refSeq, Sequence readSeq,  SAMRecord sam, List<Insertions> insertions){
		IdentityProfile profile = new IdentityProfile();
		int readLen = readSeq.length();
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
				if(length>insThresh ){
					insertions.add(new Insertions(refPos, readPos, length, readPos, readLen -(readPos+length)));
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
