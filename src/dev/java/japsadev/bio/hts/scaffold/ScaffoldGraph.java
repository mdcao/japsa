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

/**************************     REVISION HISTORY    **************************
 * 20/12/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.bio.hts.scaffold;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import japsa.seq.Alphabet;
import japsa.seq.JapsaFeature;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.Logging;
import japsa.util.Range;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.ProcessBuilder.Redirect;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;

public abstract class ScaffoldGraph{
	public static volatile int maxRepeatLength=7500; //for ribosomal repeat cluster in bacteria (Koren S et al 2013), it's 9.1kb for yeast.
	public static volatile int marginThres = 1000;
	public static volatile int minContigLength = 200;
	public static volatile int illuminaReadLength = 300; //Illumina MiSeq 2x300bp
	public static volatile int minSupportReads = 1;
	public static volatile boolean verbose = false;
	public static volatile boolean reportAll = true;
	public static volatile boolean updateGenome = true;
	public static volatile boolean eukaryotic = false;
	public static volatile boolean select = false;

	
	public volatile boolean annotation = false;
	public static volatile byte assembler =0b00; // 0 for SPAdes, 1 for ABySS
	
	public static HashMap<Integer,Integer> countOccurence=new HashMap<Integer,Integer>();
	
	public String prefix = "out";					
	public static double estimatedCov = 0;
	private static double estimatedLength = 0;

	//below maps contain avatar of contigs and bridges only, 
	//not the actual ones in used (because of the repeats that need to be cloned)

	ArrayList<Contig> contigs;
	HashMap<String, ContigBridge> bridgeMap= new HashMap<String, ContigBridge>();
	static HashMap<Integer, ArrayList<ContigBridge>> bridgesFromContig = new HashMap<Integer, ArrayList<ContigBridge>>();

	static HashMap<Integer, Range> contigsRange = new HashMap<Integer, Range>(); // for each contigs, save a pair of cooridnates (i,j)
	//for which, (0,i) and (j,contig.length()-1) are of low confident
	// so that an useful alignment must spanning these regions
	
	Scaffold [] scaffolds; // DNA translator, previous image of sequence is stored for real-time processing
	int scfNum, cirNum; // assembly statistics: number of contigs and circular ones.
	
	// Constructor for the graph with contigs FASTA file (contigs.fasta from SPAdes output)
	public ScaffoldGraph(String sequenceFile) throws IOException{
		//1. read in contigs
		SequenceReader reader = SequenceReader.getReader(sequenceFile);
		Sequence seq;
		contigs = new ArrayList<Contig>(); 

		int index = 0;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			Contig ctg = new Contig(index, seq);		

			String 	name = seq.getName(),
					desc = seq.getDesc();
			
			double mycov = 0;
			//SPAdes header: >%d_length_%d_cov_%f
			//SPAdes header: >MicroManage_ID %d_length_%d_cov_%f
			if(assembler==0b00){
				String [] toks = (name+desc).split("_");
				for (int i = 0; i < toks.length - 1;i++){
					if ("cov".equals(toks[i])){
						mycov = Double.parseDouble(toks[i+1]);
						break;
					}
				}
			} 
			//ABySS header: >%d %d %d, ID, length, kmer_sum (number of total mapped kmers)
			else if(assembler==0b01){
				String [] toks = desc.split("\\s");
				mycov = Double.parseDouble(toks[1])/(Double.parseDouble(toks[0])-Graph.getKmerSize());
			}
			
			/////////////////////////////////////////////////////////////////////
			/**
			*convert from k-mer coverage (read kmers per contig kmers) 
			*to read coverage (read bp per contig bp):
			*Cx=Ck*L/(L-k+1)
			*/
			System.out.printf("%s before=%.2f",ctg.getName(), mycov);
			mycov=mycov*illuminaReadLength/(illuminaReadLength-Graph.getKmerSize()+1);

			estimatedCov += mycov * seq.length();
			estimatedLength += seq.length();
			
			System.out.printf(" after=%.2f length=%d ave=%.2f\n", mycov, seq.length(), estimatedCov/estimatedLength);

			ctg.setCoverage(mycov);

			contigs.add(ctg);
			bridgesFromContig.put(ctg.getIndex(), new ArrayList<ContigBridge>());
			
			contigsRange.put(ctg.getIndex(), new Range(0,ctg.length()-1));
			
			index ++;
		}
		reader.close();
		
		estimatedCov /= estimatedLength;
		

		Logging.info("Read length:" + illuminaReadLength + "bp Kmer: " + Graph.getKmerSize() + "bp Average coverage:" + estimatedCov + " Length: " + estimatedLength);
		//turn off verbose mode if the genome is bigger than 100Mb. 
		if(estimatedLength > 100000000 || contigs.size() > 10000){
//			Logging.warn("Verbose mode and realtime updating are automatically disabled due to too large genome!");
//			verbose=false;
			updateGenome=false;
		}
				
		//2. Initialise scaffold graph
		scaffolds = new Scaffold[contigs.size()];		
		
		for (int i = 0; i < contigs.size();i++){				
			scaffolds[i] = new Scaffold(contigs.get(i));
			//point to the head of the scaffold
			contigs.get(i).head = i;
		}//for

		scfNum=contigs.size();
		cirNum=0;
		 
		
		////////////////////////////////////////////////////////////////////////
		/**
		 * A-statistics: A(delta,r)=log(e)*delta*n/G -r*log(2)
		 * delta: contig length, r: number of reads comprise this contig, 
		 * n: total number of reads, G: genome size
		 * Recalculated by Cx, contig_len, read_len, estimatedCov (average read coverage over the genome)
		 */
		for(Contig ctg:contigs){
			double astats = ctg.length()*estimatedCov/illuminaReadLength - Math.log(2)*ctg.getCoverage()*ctg.length()/illuminaReadLength;
			astats*=Math.log10(Math.E);
			ctg.setAstatistics(astats);
			//Logging.info(ctg.getName() + " A_stat=" + astats + " isRepeat=" + isRepeat(ctg));
			
			ctg.setRepeatFlag(astats>5?false:true);
		}

	}//constructor

	public String getAssemblerName(){
		if(assembler==0b01)
			return new String("ABySS");
		else
			return new String("SPAdes");
	}
	/* Read short-read assembly information from SPAdes output: assembly graph (assembly_graph.fastg) and 
	 ** traversed paths (contigs.pahth) to make up the contigs
	 */ 
	public void readMore(String assemblyGraph, String paths) throws IOException{	
		//1. Read assembly graph and store in a string graph
		Graph g = new Graph(assemblyGraph,assembler);

		//		for(Vertex v:g.getVertices()){
		//			System.out.println("Neighbors of vertex " + v.getLabel() + " (" + v.getNeighborCount() +"):");
		//			for(Edge e:v.getNeighbors())
		//				System.out.println(e + "; ");
		//			System.out.println();
		//		}

		Contig.setGraph(g);

		if(assembler==0b00){
			//2. read file contigs.paths from SPAdes
			BufferedReader pathReader = new BufferedReader(new FileReader(paths));
	
			String s;
			//Read contigs from contigs.paths and refer themselves to contigs.fasta
			Contig curContig = null;
			while((s=pathReader.readLine()) != null){
				if(s.contains("NODE"))
					curContig=getSPadesContig(s);
				else if(curContig!=null)
					curContig.setPath(new Path(g,s));
	
			}
			pathReader.close();
		} else if(assembler==0b01){ //for the case of ABySS: contig and vertex of assembly graph are the same
			for(Contig ctg:contigs)
				ctg.setPath(new Path(g,ctg.getName()+"+"));
		}
		
		if(ScaffoldGraph.verbose)
			Logging.info("Short read assembler " + (assembler==0b00?"SPAdes":"ABySS") + " kmer=" + Graph.getKmerSize());
	}

	public Contig getSPadesContig(String name){
		if(name.contains("'"))
			return null;


		Contig res = null;

		for(Contig ctg:contigs){
			// Extract to find contig named NODE_x_ 
			//because sometimes there are disagreement between contig name (_length_) in contigs.paths and contigs.fasta in SPAdes!!!
			// 
			if( ctg.getName().contains(name) || (ctg.getName()+ctg.getDesc()).contains("NODE_"+name.split("_")[1]+"_") ){
				res = ctg;
				break;
			}
		}

		if(res==null && verbose){
			System.out.println("Contig not found:" + name);
		}

		return res;
	}

	public synchronized int getN50(){
		int [] lengths = new int[scaffolds.length];
		int count=0;
		double sum = 0;
//		for (int i = 0; i < scaffolds.length;i++){
//			if(scaffolds[i].isEmpty()) continue;
//			int len = scaffolds[i].length();
//			
//			if(contigs.get(i).head == i)
//				if (	(!isRepeat(contigs.get(i)) && len > maxRepeatLength) //here are the big ones
//						|| scaffolds[i].closeBridge != null //circular plasmid contigs
//						|| (reportAll && needMore(contigs.get(i)) && contigs.get(i).coverage > .5*estimatedCov)) //short,repetitive sequences here if required	
//				{
//					lengths[count] = len;
//					sum+=len;
//					count++;
//				}
//		}
		for (int i = 0; i < scaffolds.length;i++){
			if(scaffolds[i].isEmpty()) continue;
			
//			if(select 
//					&& !contigs.get(i).isMapped()
//				)
//				continue;
			
			int len = scaffolds[i].length();
			
			if(contigs.get(i).head == i){
				if (!isMarker(contigs.get(i)) && !reportAll ) 
					continue;
			}
			else if(isMarker(contigs.get(i)) || !needMore(contigs.get(i)) )
				continue;
			
			lengths[count] = len;
			sum+=len;
			count++;		
		}
		Arrays.sort(lengths);

		int index = lengths.length;
		double contains = 0;
		while (contains < sum/2){
			index --;
			contains += lengths[index];
		}

		return lengths[index];	
	}

	public synchronized String getGapsInfo(){
		int 	gapCount=0,
				gapMaxLen=0;

		for (int i = 0; i < scaffolds.length;i++){
			if(scaffolds[i].isEmpty()) continue;
			int len = scaffolds[i].length();
			if ((contigs.get(i).head == i 
					&& isMarker(contigs.get(i))
					&& len > maxRepeatLength
					)
					|| scaffolds[i].closeBridge != null)
			{
				for(ContigBridge brg:scaffolds[i].bridges){
					if(brg.getBridgePath()==null){
						gapCount++;
						if(brg.getTransVector().distance(brg.firstContig, brg.secondContig) > gapMaxLen)
							gapMaxLen=brg.getTransVector().distance(brg.firstContig, brg.secondContig);
					}
				}
				if(scaffolds[i].closeBridge!=null){
					ContigBridge brg=scaffolds[i].closeBridge;
					if(brg.getBridgePath()==null){
						gapCount++;
						if(brg.getTransVector().distance(brg.firstContig, brg.secondContig) > gapMaxLen)
							gapMaxLen=brg.getTransVector().distance(brg.firstContig, brg.secondContig);
					}
				}

			}
		}


		return gapCount+" ("+gapMaxLen+")";	
	}
	/**
	 * MDC added second version that include bwa
	 * @param bamFile
	 * @param minCov
	 * @param qual
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void makeConnectionsWithBWA(String inFile, double minCov, int qual, String format, String bwaExe, int bwaThread, String bwaIndex) throws IOException, InterruptedException{

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);

		SamReader reader = null;
		Process bwaProcess = null;

		if (format.endsWith("am")){//bam or sam
			if ("-".equals(inFile))
				reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				reader = SamReaderFactory.makeDefault().open(new File(inFile));	
		}else{
			Logging.info("Starting bwa  at " + new Date());

			ProcessBuilder pb = null;
			if ("-".equals(inFile)){
				pb = new ProcessBuilder(bwaExe, 
						"mem",
						"-t",
						"" + bwaThread,
						"-k11",
						"-W20",
						"-r10",
						"-A1",
						"-B1",
						"-O1",
						"-E1",
						"-L0",
						"-a",
						"-Y",
//						"-K",
//						"20000",
						bwaIndex,
						"-"
						).
						redirectInput(Redirect.INHERIT);
			}else{
				pb = new ProcessBuilder(bwaExe, 
						"mem",
						"-t",
						"" + bwaThread,
						"-k11",
						"-W20",
						"-r10",
						"-A1",
						"-B1",
						"-O1",
						"-E1",
						"-L0",
						"-a",
						"-Y",
//						"-K",
//						"20000",
						bwaIndex,
						inFile
						);
			}
			bwaProcess  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();

			//Logging.info("bwa started x");	
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bwaProcess.getInputStream()));
		}


		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		ReadFilling readFilling = null;
		ArrayList<AlignmentRecord> samList = null;// alignment record of the same read;	
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();			
			
			if (rec.getReadUnmappedFlag())
				continue;
			if (rec.getMappingQuality() < qual)
				continue;
			
			Contig tmp = contigs.get(rec.getReferenceIndex());
			if(tmp==null){
				Logging.error("Contig " + rec.getReferenceIndex() + " doesn't exist!");
				System.exit(1);
			}
				
			AlignmentRecord myRec = new AlignmentRecord(rec, tmp);
//			Arrays.fill(tmp.isMapped, myRec.refStart, myRec.refEnd, 1);
			
//			System.out.println("Processing record of read " + rec.getReadName() + " and ref " + rec.getReferenceName() + (myRec.useful?": useful ":": useless ") + myRec);


			//////////////////////////////////////////////////////////////////
			// make bridge of contigs that align to the same (Nanopore) read. 
			// Note that SAM file MUST be sorted based on readID (samtools sort -n)
			//	which is natural if it is the output from an aligner (bwa, minimap2)

			//not the first occurrance				
			if (readID.equals(myRec.readID)) {				
				if (myRec.useful){				
					for (AlignmentRecord s : samList) {
						if (s.useful){
							this.addBridge(readFilling, s, myRec, minCov); //stt(s) < stt(myRec) -> (s,myRec) appear once only!
						}
					}
				}
			} else{
			
//			if (!readID.equals(myRec.readID)){
//				//process samlist here
//				processAlignments(samList);
					
					
				samList = new ArrayList<AlignmentRecord>();
				readID = myRec.readID;	
				readFilling = new ReadFilling(new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + readID), samList);	
			}			
			samList.add(myRec);

		}// while
		iter.close();

		//outOS.close();
		reader.close();		
		if (bwaProcess != null){
			bwaProcess.waitFor();
		}

	}
	/**
	 * SHN added minimap2 as the default aligner
	 * @param bamFile
	 * @param minCov
	 * @param qual
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void makeConnectionsWithMinimap2(String inFile, double minCov, int qual, String format, String mm2Preset, int mm2Threads, String mm2Index) throws IOException, InterruptedException{

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);

		SamReader reader = null;
		Process mm2Process = null;

		if (format.endsWith("am")){//bam or sam
			if ("-".equals(inFile))
				reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				reader = SamReaderFactory.makeDefault().open(new File(inFile));	
		}else{
			Logging.info("Starting minimap2 at " + new Date());

			ProcessBuilder pb = null;
			if ("-".equals(inFile)){
				pb = new ProcessBuilder("minimap2", 
						"-t",
						"" + mm2Threads,
						"-ax",
						mm2Preset,
						"-I",
						"40g",
						mm2Index,
						"-"
						).
						redirectInput(Redirect.INHERIT);
			}else{
				pb = new ProcessBuilder("minimap2", 
						"-t",
						"" + mm2Threads,
						"-ax",
						mm2Preset,
						"-I",
						"40g",
						mm2Index,
						inFile
						);
			}
			mm2Process  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();

			//Logging.info("bwa started x");	
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(mm2Process.getInputStream()));
		}


		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		ReadFilling readFilling = null;
		ArrayList<AlignmentRecord> samList = null;// alignment record of the same read;	
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();			
			
			if (rec.getReadUnmappedFlag())
				continue;
			if (rec.getMappingQuality() < qual)
				continue;
			
			Contig tmp = contigs.get(rec.getReferenceIndex());
			if(tmp==null){
				Logging.error("Contig " + rec.getReferenceIndex() + " doesn't exist!");
				System.exit(1);
			}
				
			AlignmentRecord myRec = new AlignmentRecord(rec, tmp);
//			Arrays.fill(tmp.isMapped, myRec.refStart, myRec.refEnd, 1);
			
//			System.out.println("Processing record of read " + rec.getReadName() + " and ref " + rec.getReferenceName() + (myRec.useful?": useful ":": useless ") + myRec);


			//////////////////////////////////////////////////////////////////
			// make bridge of contigs that align to the same (Nanopore) read. 
			// Note that SAM file MUST be sorted based on readID (samtools sort -n)
			//	which is natural if it is the output from an aligner (bwa, minimap2)

			//not the first occurrance				
			if (readID.equals(myRec.readID)) {				
				if (myRec.useful){				
					for (AlignmentRecord s : samList) {
						if (s.useful){
							this.addBridge(readFilling, s, myRec, minCov); //stt(s) < stt(myRec) -> (s,myRec) appear once only!
						}
					}
				}
			} else{
			
//			if (!readID.equals(myRec.readID)){
//				//process samlist here
//				processAlignments(samList);
					
					
				samList = new ArrayList<AlignmentRecord>();
				readID = myRec.readID;	
				readFilling = new ReadFilling(new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + readID), samList);	
			}			
			samList.add(myRec);

		}// while
		iter.close();

		//outOS.close();
		reader.close();		
		if (mm2Process != null){
			mm2Process.waitFor();
		}

	}
	
	//Take a list of alignments between contigs and a single nanopore read
	// do some interesting things
//	private void processAlignments(ArrayList<AlignmentRecord> list){
//		if(list==null)
//			return;
//		
//		//grouping based on contigs ID
//		HashMap<String,ArrayList<AlignmentRecord>> groups = new HashMap<String,ArrayList<AlignmentRecord>>();
//		for(AlignmentRecord rec:list){
//			String ctgID=rec.contig.getName();
//			ArrayList<AlignmentRecord> tmpList=groups.get(ctgID);
//			if(tmpList==null)
//				tmpList=new ArrayList<AlignmentRecord>();
//			else
//				tmpList.add(rec);
//		}
//		
//		//process each group
//		for (String key:groups.keySet()){
//			ArrayList<AlignmentRecord> tmpList = groups.get(key);
//			if(tmpList.size()>1){
//				for(AlignmentRecord rec:tmpList){
//					int x1=rec.refStart, x2=rec.refEnd,
//						y1=rec.readStart, y2=rec.readEnd;
//					
//					double 	intersect=(x1*y2-x2*y1)/(y2-y1), // ref(x1,x2) <-> read(y1,y2): (x1y2-x2y1)/(y2-y1)
//							angle=(x2-x1)/(y2-y1); //cos
//					
//				}
//			}
//		}
//	}




	/*********************************************************************************/
	protected void addBridge(ReadFilling readSequence, AlignmentRecord a, AlignmentRecord b, double minCov){
		if (a.contig.index > b.contig.index){
			AlignmentRecord t = a;a=b;b=t;
		}
		// Rate of aligned lengths: ref/read (illumina contig/nanopore read)
		int 	alignedReadLen = Math.abs(a.readEnd - a.readStart) + Math.abs(b.readEnd - b.readStart),
				alignedRefLen = Math.abs(a.refEnd - a.refStart) + Math.abs(b.refEnd - b.refStart);
		double rate = 1.0 * alignedRefLen/alignedReadLen;		

		//See if this is reliable
		double score = Math.min(a.score, b.score);
		int alignP = (int) ((b.readStart - a.readStart) * rate);
		int alignD = (a.strand == b.strand)?1:-1;

		//(rough) relative position from ref_b (contig of b) to ref_a (contig of a) in the assembled genome
		int gP = (alignP + (a.strand ? a.refStart:-a.refStart) - (b.strand?b.refStart:-b.refStart));
		if (!a.strand)
			gP = -gP;	
		
		if (	a.contig.getIndex() == b.contig.getIndex() ){
			if(
				alignD > 0
				&& (Math.abs(gP)*1.0 / a.contig.length()) < 1.1 
				&& (Math.abs(gP)*1.0 / a.contig.length()) > 0.9 
//				&& a.readLength < 1.1* a.contig.length()
			)
			{
				if(a.readLength > 1.1* a.contig.length())
					a.contig.cirProb = Double.MIN_VALUE; //tandem not circular
				else if(alignedReadLen*1.0/a.readLength > 0.8 ){ //need more than 80% alignment (error rate of nanopore read)
					a.contig.cirProb += score;	//more likely to be circular
					
					if(verbose) 
						System.out.printf("Potential CIRCULAR contig %s map to read %s(length=%d): (%d,%d) => circular score: %.2f\n"
								, a.contig.getName(), a.readID, a.readLength, gP, alignD, a.contig.cirProb);

				}
				
				
				
			}	
		}else{
//			if(verbose)
//				System.out.println("==> substract cirProb -" + score + " from " + a.contig.getIndex() + " and " + b.contig.getIndex());
			a.contig.cirProb-=score;
			b.contig.cirProb-=score;
		}
		
		// overlap length on aligned read (<0 if not overlap)
		int overlap = Math.min(	a.readAlignmentEnd() - b.readAlignmentStart(), b.readAlignmentEnd() - a.readAlignmentStart());				

		if (	overlap > Math.min(	.5 * Math.min(a.readAlignmentEnd()-a.readAlignmentStart(), b.readAlignmentEnd()-b.readAlignmentStart()),
				minContigLength)      
				|| a.contig.getCoverage() < minCov	// filter out contigs with inappropriate cov
				|| b.contig.getCoverage() < minCov
				)
		{		
//			if(verbose)
//				System.out.println("...ignoring " + a.contig.getName() + "#" + b.contig.getName() + " (overlap=" + overlap + ")");
			return;
		}

		ScaffoldVector trans = new ScaffoldVector(gP, alignD);		

		int count = 0;
		ContigBridge bridge, bridge_rev;
		while (true){
			int brgID = count, revID = count;
			if(a.contig.getIndex()==b.contig.getIndex()){
				brgID = 2*count;
				revID = brgID+1;
			}
			String 	hash = ContigBridge.makeHash(a.contig.index, b.contig.index, brgID),
					hash_rev = ContigBridge.makeHash(b.contig.index, a.contig.index, revID);

			bridge = bridgeMap.get(hash);
			bridge_rev = bridgeMap.get(hash_rev);
			if (bridge == null){
				assert bridge_rev==null:hash_rev + " not null!";
				bridge = new ContigBridge(a.contig, b.contig, brgID);
				bridge_rev = new ContigBridge(b.contig, a.contig, revID);

				bridge.addConnection(readSequence, a, b, trans, score);
				bridge_rev.addConnection(readSequence, b, a, ScaffoldVector.reverse(trans), score);

//				if(verbose)
//					System.out.println("...creating " + bridge.hashKey + " and " + bridge_rev.hashKey);
				bridgesFromContig.get(a.contig.getIndex()).add(bridge);
				bridgesFromContig.get(b.contig.getIndex()).add(bridge_rev);

				bridgeMap.put(hash, bridge);
				
				if(trans.equalsTo(ScaffoldVector.reverse(trans)) && hash.equals(hash_rev))
					break;
				
				bridgeMap.put(hash_rev, bridge_rev);

				break;
			}
			if ((a.contig.getIndex() != b.contig.getIndex()) && bridge.consistentWith(trans)){
				assert bridge_rev!=null:hash_rev + "is null!";
				bridge.addConnection(readSequence, a, b, trans, score);
				bridge_rev.addConnection(readSequence, b, a, ScaffoldVector.reverse(trans), score);
//				if(verbose)
//					System.out.println("...adding connections to " + bridge.hashKey + " and " + bridge_rev.hashKey);
				
				break;
			}
			if(a.contig.getIndex() == b.contig.getIndex()){
				assert bridge_rev!=null:hash_rev + "is null";
				
				if(bridge.consistentWith(trans)){
					bridge.addConnection(readSequence, a, b, trans, score);
					if(!trans.equalsTo(ScaffoldVector.reverse(trans)))
						bridge_rev.addConnection(readSequence, b, a, ScaffoldVector.reverse(trans), score);
					break;
				}
				if(bridge.consistentWith(ScaffoldVector.reverse(trans))){
					bridge.addConnection(readSequence, b, a, ScaffoldVector.reverse(trans), score);
					if(!trans.equalsTo(ScaffoldVector.reverse(trans)))
						bridge_rev.addConnection(readSequence, a, b, trans, score);		
					break;
				}		

//				if(verbose)
//					System.out.println("...adding connections to " + bridge.hashKey + " and " + bridge_rev.hashKey);
			}
			count ++;
		}//while

	}
	public static ArrayList<ContigBridge> getListOfBridgesFromContig(Contig ctg){
		return bridgesFromContig.get(ctg.getIndex());
	}
	/**********************************************************************************************/
	public ContigBridge getReversedBridge(ContigBridge bridge){
		String hash = ContigBridge.makeHash(bridge.secondContig.index, bridge.firstContig.index, bridge.orderIndex);
		return bridgeMap.get(hash);
	}
	/**********************************************************************************************/
	/*
	 * Check if it's possible to extend from *contig* with *bridge* to another extended-already contig (contigF) 
	 * use for markers and unique bridge only. 
	 * This is a pre-step to join 2 scaffolds: scaffoldT going to scaffoldF
	 * @param 	Contig: a contig to start with
	 * 			ContigBridge: a bridge from given contig to a candidate unique contig for the extension
	 * @return	int: direction on targeted scaffold (scaffoldF) that can be traversed
	 */
	protected int extendDirection(Contig contig, ContigBridge bridge){
		Contig contigF = bridge.secondContig;
		ScaffoldVector trans = bridge.getTransVector(); //contig->contigF 
		int pointer = Integer.signum(trans.magnitude * trans.direction); //pointer < 0 => tail of contigF on bridge
		assert scaffolds[contigF.head].size() > 1 : contigF.head;

		int headF = contigF.head;
		int direction = 0; //direction of extension on scaffoldT (we need to return direction on scaffoldF)
		ScaffoldVector headT2contigF = ScaffoldVector.composition(trans, contig.getVector());
		int rEnd = contig.rightMost(), rEndF = contigF.rightMost(headT2contigF),
				lEnd =  contig.leftMost(), lEndF = contigF.leftMost(headT2contigF);
		if(rEndF > rEnd){
			direction = 1;
		}
		else if(lEndF < lEnd){
			direction = -1;
		}
		else 
			return 0;
		if(verbose)
			System.out.println("Examining extending direction from contig " + contig.getIndex() + " to " + bridge.hashKey);
		Scaffold scaffoldF = scaffolds[headF];
		// Get order-based (order on scaffold other than orientation-based of contig) previous and next marker(unique contig)
		Contig 	prevMarker = scaffoldF.nearestMarker(contigF, false), // previous marker of contigF on *corresponding scaffold*
				nextMarker = scaffoldF.nearestMarker(contigF, true); // next marker of contigF on *corresponding scaffold*		

		ScaffoldVector rev = ScaffoldVector.reverse(contigF.getVector()); //rev = contigF->headF	
		if(prevMarker != null){
			ScaffoldVector toPrev = ScaffoldVector.composition(prevMarker.getVector(),rev); //contigF->prevMarker
			if(scaffoldF.indexOf(prevMarker) > scaffoldF.indexOf(contigF) && scaffoldF.closeBridge != null)
				//toPrev = ScaffoldVector.composition(ScaffoldVector.reverse(scaffoldF.circle), toPrev);
				toPrev = scaffoldF.rotate(toPrev, false);
			ScaffoldVector headT2Prev = ScaffoldVector.composition(toPrev, headT2contigF);
			int rEndPrev = prevMarker.rightMost(headT2Prev),
					lEndPrev = prevMarker.leftMost(headT2Prev);
			if(verbose){
				System.out.printf("Extending from contigT %d to targeted contig (contigF) %d with previous contig (prevMarker) %d \n", contig.getIndex(), contigF.getIndex(), prevMarker.getIndex());
				System.out.println("...headT->contig, contigF and prevMarker: " + contig.getVector() + headT2contigF + headT2Prev);
			}
			if ((direction > 0?rEndPrev > rEndF: lEndPrev < lEndF)){
				//check if the candidate ContigBridge is more confident than the current or not 
				if((pointer<0?contigF.nextScore:contigF.prevScore) < bridge.getScore()){
					if(verbose)
						System.out.printf("=> go from %d to %d to %d backward\n", contig.getIndex(), contigF.getIndex(), prevMarker.getIndex());
					return -1;
				}
				else{
					if(verbose)
						System.out.printf("Bridge score not strong enough: %.2f < %.2f (%.2f)\n", 
								bridge.getScore(), pointer<0?contigF.nextScore:contigF.prevScore,
										pointer<0?contigF.prevScore:contigF.nextScore);

					return 0;
				}
			}else{
				if(verbose)
					System.out.printf("Direction conflict: %d, %d %d or %d %d. Checking otherway... \n", direction, rEndPrev, rEndF, lEndPrev, lEndF);
			}
		}
		if(nextMarker != null){
			ScaffoldVector toNext = ScaffoldVector.composition(nextMarker.getVector(),rev); //contigF->nextMarker
			if(scaffoldF.indexOf(nextMarker) < scaffoldF.indexOf(contigF) && scaffoldF.closeBridge != null)
				//toNext = ScaffoldVector.composition(scaffoldF.circle, toNext);
				toNext = scaffoldF.rotate(toNext, true);
			ScaffoldVector headT2Next = ScaffoldVector.composition(toNext, headT2contigF);
			int rEndNext = nextMarker.rightMost(headT2Next),
					lEndNext = nextMarker.leftMost(headT2Next);
			if(verbose){
				System.out.printf("Extending from contigT %d to targeted contig (contigF) %d with next contig (nextMarker) %d \n", contig.getIndex(), contigF.getIndex(), nextMarker.getIndex());
				System.out.println("...headT->contig, contigF and nextMarker: " + contig.getVector() + headT2contigF + headT2Next);
			}

			if ((direction > 0? rEndNext > rEndF : lEndNext < lEndF)){
				//if((rev.direction<0?contigF.nextScore:contigF.prevScore) < bridge.getScore()){
				if((pointer<0?contigF.nextScore:contigF.prevScore) < bridge.getScore()){	
					if(verbose)
						System.out.printf("=> go from %d to %d to %d forward\n", contig.getIndex(), contigF.getIndex(), nextMarker.getIndex());
					return 1;
				}
				else{
					if(verbose)
						System.out.printf("Bridge score not strong enough: %.2f < %.2f (%.2f)\n", 
								bridge.getScore(), pointer<0?contigF.nextScore:contigF.prevScore,
										pointer<0?contigF.prevScore:contigF.nextScore);
					return 0;

				}
			}else{
				if(verbose)
					System.out.printf("Direction conflict: %d, %d %d or %d %d. End searching! \n", direction, rEndNext, rEndF, lEndNext, lEndF);
			}
		}
		return 0;
	}
	/*********************************************************************************/
	public synchronized boolean joinScaffold(Contig contig, ContigBridge bridge, boolean firstDir, int secondDir){		
		if(verbose) {
			System.out.println("PREPARE TO CONNECT " + bridge.hashKey + " with score " + bridge.getScore() + 
					", size " + bridge.getNumOfConnections() + 
					", vector (" + bridge.getTransVector().toString() + 
					"), distance " + bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig));
			bridge.display();
		}


		Contig contigF = bridge.secondContig, contigT = contig;
		ScaffoldVector trans = bridge.getTransVector();

		int headF = contigF.head,
			headT = contigT.head;
		Scaffold 	scaffoldF = scaffolds[headF],
					scaffoldT = scaffolds[headT];
		int	posT = scaffoldT.isEnd(contigT);
		if (posT == 0){
			if(verbose) 
				System.out.println("Impossible to jump from the middle of a scaffold " + headT + ": contig " + contigT.index);
			return false;
		}

		if(verbose) 
			System.out.println("Before joining " + contigF.index + " (" + headF +") to " + contigT.index 
					+ " (" + headT +") " 
					+ (scaffoldT.getLast().rightMost() - scaffoldT.getFirst().leftMost()) 
					+ "+" + (scaffoldF.getLast().rightMost() - scaffoldF.getFirst().leftMost()) 
					+ "=" + (scaffoldT.getLast().rightMost() - scaffoldT.getFirst().leftMost() + scaffoldF.getLast().rightMost() - scaffoldF.getFirst().leftMost()));
		//===================================================================================================
		int index = scaffoldF.indexOf(contigF),
				count = index;

		ScaffoldVector rev = ScaffoldVector.reverse(contigF.getVector()); //rev = contigF->headF	

		int addScf=-1; 

		if(secondDir == -1){
			if(headF==headT){

				if(firstDir){
					//            .........
					//       <----v       :
					// ===========*=======>
					Contig prevMarker = scaffoldF.nearestMarker(contigF, false);
					
					if(prevMarker!=null){
						Contig ctg = scaffoldF.remove(count--);
						ContigBridge brg = scaffoldF.bridges.remove(count+1);
						Scaffold newScf = new Scaffold(ctg);
						
						while(true){
							if(count<0) break;
							ctg= scaffoldF.remove(count--);
							brg = scaffoldF.bridges.remove(count+1);
							newScf.addBackward(ctg,brg); //must be here, after the assignments of ctg, brg???

							
						}
						changeHead(newScf, prevMarker);
					}
					joinScaffold(contig, bridge, true, -1);					

				}
				else{
					// ..............
					// :       <----v
					// <============*=======
					Contig nextMarker = scaffoldF.nearestMarker(contigF, true);
					if(nextMarker!=null){
						Contig ctg = scaffoldF.remove(index+1);
						Scaffold newScf = new Scaffold(ctg);
						ContigBridge brg = scaffoldF.bridges.remove(index);
						while(true){
							if(scaffoldF.size()==index+1) break;
							ctg= scaffoldF.remove(index+1);
							brg = scaffoldF.bridges.remove(index);
							newScf.addForward(ctg,brg); //must be here, after the assignments of ctg, brg???
						}
						newScf.trim();
						changeHead(newScf, nextMarker);
						addScf=nextMarker.getIndex();
					}
					scaffoldF.setCloseBridge(getReversedBridge(bridge));
					changeHead(scaffoldF, contigF);
				}
				
				
			}else{
				Contig 	ctg = scaffoldF.remove(index);
				ContigBridge brg = getReversedBridge(bridge);
				//extend and connect
				while(true){
					ctg.composite(rev); // contigF->headF + headF->ctg = contigF->ctg
					ctg.composite(trans); // contigT->contigF + contigF->ctg = contigT->ctg
					ctg.composite(contigT.getVector()); //headT->contigT + contigT->ctg = headT->ctg : relative position of this ctg w.r.t headT

					ctg.head = headT;
					//if (posT == 1){
					if(!firstDir){
						scaffoldT.addBackward(ctg,brg);
					}else{
						scaffoldT.addForward(ctg,getReversedBridge(brg));
					}	
					if(count<1) break;
					ctg = scaffoldF.remove(--count);
					brg = scaffoldF.bridges.remove(count);

				}
				if(scaffoldF.closeBridge!=null && !scaffoldF.isEmpty()){
					count = scaffoldF.size()-1;
					ctg = scaffoldF.removeLast();
					brg = scaffoldF.closeBridge;

					while(true){
						//ctg.myVector = ScaffoldVector.composition(ScaffoldVector.reverse(scaffoldF.circle),ctg.myVector);
						ctg.myVector = scaffoldF.rotate(ctg.myVector, false);
						ctg.composite(rev); // contigF->headF + headF->ctg = contigF->ctg
						ctg.composite(trans); // contigT->contigF + contigF->ctg = contigT->ctg
						ctg.composite(contigT.getVector()); //headT->contigT + contigT->ctg = headT->ctg : relative position of this ctg w.r.t headT
						//ctg.composite(ScaffoldVector.reverse(scaffoldF.circle)); //composite co tinh giao hoan k ma de day???
						ctg.head = headT;
						//if (posT == 1){ 
						if(!firstDir){
							scaffoldT.addBackward(ctg,brg);
						}else{
							scaffoldT.addForward(ctg,getReversedBridge(brg));
						}	
						if(count<1) break;
						brg = scaffoldF.bridges.remove(count--);
						ctg = scaffoldF.remove(count);		

					}
				}

				//set the remaining.
				scaffoldT.trim();
				scaffoldF.trim();
				if(!scaffoldF.isEmpty()){
					addScf=scaffoldF.getFirst().getIndex();//getFirst: NoSuchElementException
					changeHead(scaffoldF, scaffoldF.getFirst());
				}
			}
			//now since scaffoldF is empty due to changeHead(), re-initialize it!(do we need this??)
			scaffoldF = new Scaffold(contigs.get(headF));
		}
		else if(secondDir == 1){
			if(headF==headT){
				if(!firstDir){
					// ...........
					// :         v--->
					// <=========*=========			
					Contig nextMarker = scaffoldF.nearestMarker(contigF, true);
					
					if(nextMarker!=null){
						Contig ctg = scaffoldF.remove(index);
						ContigBridge brg = scaffoldF.bridges.remove(index-1);
						Scaffold newScf = new Scaffold(ctg);
						
						while(true){
							if(scaffoldF.size()==index) break;
							ctg= scaffoldF.remove(index);
							brg = scaffoldF.bridges.remove(index-1);
							newScf.addForward(ctg,brg); //must be here, after the assignments of ctg, brg???
							
						}
						changeHead(newScf, nextMarker);
					}
					joinScaffold(contig, bridge, false, 1);	
					
				}
				else{
					//       .............
					//       v--->       :
					// ======*===========>
					Contig prevMarker = scaffoldF.nearestMarker(contigF, false);
					if(prevMarker!=null){
						Contig ctg = scaffoldF.remove(--count);
						Scaffold newScf = new Scaffold(ctg);
						ContigBridge brg = scaffoldF.bridges.remove(count);
						while(true){
							if(count<1) break;
							ctg= scaffoldF.remove(--count);
							brg = scaffoldF.bridges.remove(count);
							newScf.addBackward(ctg,brg); //must be here, after the assignments of ctg, brg???
						}
						newScf.trim();
						changeHead(newScf, prevMarker);
						addScf=prevMarker.getIndex();
					}
					scaffoldF.setCloseBridge(bridge);
					changeHead(scaffoldF, contigF);

				}
			}else{
				System.out.println("Before removing index=" + index);
				scaffoldF.view();
				Contig 	ctg = scaffoldF.remove(index);
				ContigBridge brg = bridge;
				//extend and connect
				while(true){
					ctg.composite(rev); // contigF->headF + headF->ctg = contigF->ctg
					ctg.composite(trans); // contigT->contigF + contigF->ctg = contigT->ctg
					ctg.composite(contigT.getVector()); //headT->contigT + contigT->ctg = headT->ctg : relative position of this ctg w.r.t headT

					ctg.head = headT;
					//if (posT == 1){ 
					if(!firstDir){
						scaffoldT.addBackward(ctg,getReversedBridge(brg));
					}else{
						scaffoldT.addForward(ctg,brg);
					}				
					if(scaffoldF.size()<=index) break;
					System.out.println("Before removing index=" + index);
					scaffoldF.view();
					ctg = scaffoldF.remove(index);
					brg = scaffoldF.bridges.remove(index);
				}
				if(scaffoldF.closeBridge!=null && !scaffoldF.isEmpty()){
					ctg = scaffoldF.removeFirst();
					brg = scaffoldF.closeBridge;
					while(true){
						//ctg.myVector = ScaffoldVector.composition(scaffoldF.circle,ctg.myVector);
						ctg.myVector = scaffoldF.rotate(ctg.myVector, true);
						ctg.composite(rev); // contigF->headF + headF->ctg = contigF->ctg
						ctg.composite(trans); // contigT->contigF + contigF->ctg = contigT->ctg
						ctg.composite(contigT.getVector()); //headT->contigT + contigT->ctg = headT->ctg : relative position of this ctg w.r.t headT
						//ctg.composite(scaffoldF.circle);
						ctg.head = headT;
						//if (posT == 1){ 
						if(!firstDir){
							scaffoldT.addBackward(ctg,getReversedBridge(brg));
						}else{
							scaffoldT.addForward(ctg,brg);
						}	
						if(scaffoldF.size()<1) break;
						brg = scaffoldF.bridges.removeFirst();
						ctg = scaffoldF.removeFirst();		
					}	
				}
				//set the remaining
				scaffoldT.trim();
				scaffoldF.trim();
				if(!scaffoldF.isEmpty()){
					addScf=scaffoldF.getLast().getIndex(); //getLast: NoSuchElementException
					changeHead(scaffoldF, scaffoldF.getLast());
				}
			}
			//now since scaffoldF is empty due to changeHead(), re-initialize it!(do we need this??)
			scaffoldF = new Scaffold(contigs.get(headF));
		}	
		else
			return false;

		//===================================================================================================
		if(verbose){ 
			System.out.println("After Joining: " + (addScf<0?1:2) + " scaffolds!");
			scaffolds[contigF.head].view();
			if(addScf >=0)
				scaffolds[addScf].view();
		}
		return true;
	}
	//change head of scaffold scf to newHead. 
	//This should move the content of scf to scaffolds[newHead.idx], leaving scf=null afterward
	public void changeHead(Scaffold scf, Contig newHead){	
		if(!isMarker(newHead)){
			if(verbose)
				System.out.println("Cannot use repeat as a head! " + newHead.getName());
			return;
		}
		//Scaffold scf = scaffolds[scfIndex];
		int scfIndex = scf.scaffoldIndex;
		int headPos = scf.indexOf(newHead);
		if(headPos < 0){
			if(verbose)
				System.out.printf("Cannot find contig %d in scaffold %d\n" , newHead.getIndex(), scfIndex);
			return;
		}
		Scaffold newScf = new Scaffold(newHead.getIndex());
		ScaffoldVector rev = ScaffoldVector.reverse(newHead.getVector()); //rev = newHead->head	

		if(newHead.getRelDir() == 0){
			if(verbose)
				System.out.printf("Contig %d of scaffold %d got direction 0!\n" , newHead.getIndex(), scfIndex);
			return;
		}
		else if(newHead.getRelDir() > 0){
			while(!scf.isEmpty())
				newScf.add(scf.removeFirst());
			while(!scf.bridges.isEmpty())
				newScf.bridges.add(scf.bridges.removeFirst());
			if(scf.closeBridge != null){
				newScf.closeBridge = scf.closeBridge;
				newScf.circle = scf.circle;
				//then reset these factors
				scf.closeBridge = null;
				scf.circle = null;
			}
		}
		else{
			while(!scf.isEmpty())
				newScf.add(scf.removeLast());
			while(!scf.bridges.isEmpty())
				newScf.bridges.add(getReversedBridge(scf.bridges.removeLast()));
			if(scf.closeBridge != null){
				newScf.closeBridge = getReversedBridge(scf.closeBridge);
				//newScf.circle = ScaffoldVector.reverse(scf.circle);
				newScf.circle = scf.circle; //cuz now circle is always positive

				//then reset these factors
				scf.closeBridge = null;
				scf.circle = null;
			}
		}

		for (Contig ctg:newScf){					
			ctg.composite(rev); // leftmost->head + head->ctg = leftmost->ctg
		}
		newScf.setHead(newHead.getIndex());
		scaffolds[newHead.getIndex()] = newScf;

	}
	public synchronized void printSequences(boolean allOut, boolean isBatch) throws IOException{
		//countOccurence=new HashMap<Integer,Integer>();
		int currentNumberOfContigs = 0,
			currentNumberOfCirculars = 0;	

		if(annotation){
			SequenceOutputStream aout = SequenceOutputStream.makeOutputStream(prefix+".anno.japsa");
			for (int i = 0; i < scaffolds.length;i++){
				if(scaffolds[i].isEmpty()) continue;
				
//				if(select 
//						&& !contigs.get(i).isMapped()
//					)
//					continue;
				
				int len = scaffolds[i].getLast().rightMost() - scaffolds[i].getFirst().leftMost();
				
				if(contigs.get(i).head == i ){
					if(	!reportAll && 
							((!isMarker(contigs.get(i)) && scaffolds[i].closeBridge == null) || len < maxRepeatLength)) //repetitive linear or short sequences
							continue;			
					
					if(scaffolds[i].closeBridge != null ){
						currentNumberOfCirculars++;					
					}
					currentNumberOfContigs++;
					
				}else if(reportAll && !isMarker(contigs.get(i)) && needMore(contigs.get(i)) && contigs.get(i).coverage >= .5*estimatedCov){
					currentNumberOfContigs++;
				}else{
					continue;
				}
				
				if(verbose) 
					System.out.println("Scaffold " + i + " estimated length " + len);
				if(isBatch)
					scaffolds[i].printBatchGenes();
				if(allOut)
					scaffolds[i].viewAnnotation(aout);
			}
			aout.close();
		} else{
			SequenceOutputStream 	fout = SequenceOutputStream.makeOutputStream(prefix+".fin.fasta"),
									jout = SequenceOutputStream.makeOutputStream(prefix+".fin.japsa");
//			for (int i = 0; i < scaffolds.length;i++){
//				if(scaffolds[i].isEmpty()) continue;
//				int len = scaffolds[i].getLast().rightMost() - scaffolds[i].getFirst().leftMost();
//				
//				if(contigs.get(i).head == i){
//					if(scaffolds[i].closeBridge != null ){
//						currentNumberOfContigs++;
//						currentNumberOfCirculars++;					
//					}					
//					else if ((!isRepeat(contigs.get(i)) && len > maxRepeatLength) //here are the big ones
//							|| (reportAll && needMore(contigs.get(i)) && contigs.get(i).coverage > .5*estimatedCov)) //short,repetitive sequences here if required
//						currentNumberOfContigs++;
//					else
//						continue;
//					
//					if(verbose) 
//						System.out.println("Scaffold " + i + " estimated length " + len);
//					if(allOut)
//						scaffolds[i].viewSequence(fout, jout);
//				}
//			}
			for (int i = 0; i < scaffolds.length;i++){
				if(scaffolds[i].isEmpty()) continue;
				
//				if(select 
//						&& !contigs.get(i).isMapped()
//					)
//					continue;
				
				int len = scaffolds[i].getLast().rightMost() - scaffolds[i].getFirst().leftMost();

				if(contigs.get(i).head == i ){
					if(	!reportAll && 
						((!isMarker(contigs.get(i)) && scaffolds[i].closeBridge == null) || len < maxRepeatLength)) //repetitive linear or short sequences
						continue;			
					
					if(scaffolds[i].closeBridge != null ){
						currentNumberOfCirculars++;					
					}
					currentNumberOfContigs++;
					
				}else if(reportAll && !isMarker(contigs.get(i)) && needMore(contigs.get(i)) && contigs.get(i).coverage >= .5*estimatedCov){
					currentNumberOfContigs++;
				}else{
					continue;
				}
				
				if(verbose) 
					System.out.println("Scaffold " + i + " estimated length " + len);
				if(allOut)
					scaffolds[i].viewSequence(fout, jout);
			}
			
			fout.close();
			jout.close();
		}
		scfNum=currentNumberOfContigs;
		cirNum=currentNumberOfCirculars;
	}	
	public synchronized static void oneMore(Contig ctg){
		if(countOccurence.get(ctg.getIndex())==null)
			countOccurence.put(ctg.getIndex(), 1);
		else
			countOccurence.put(ctg.getIndex(), countOccurence.get(ctg.getIndex())+1);
	}
	
	synchronized boolean needMore(Contig ctg) {
		Integer count = countOccurence.get(ctg.getIndex());
		if(count==null) return true;
		else return false; //if not occurred (Minh)
		
//		int estimatedOccurence = (int) Math.floor(ctg.coverage/estimatedCov);
//		if(estimatedOccurence <= Math.floor(.75*count))
//			return true;
//		else
//			return false;
	}

	public synchronized void printRT(long tpoint) throws IOException{
		for (Contig contig:contigs){
			if(contig.oriRep.size() > 0){
				String fname = contig.getName() + ".rtout";
				File f = new File(fname);
				if(!f.exists())
					f.createNewFile();

				//BufferedWriter out = new BufferedWriter(new FileWriter(f.getPath(), true));
				FileWriter fw = new FileWriter(f,true);
				BufferedWriter bw = new BufferedWriter(fw);
				PrintWriter pw = new PrintWriter(bw);

				ArrayList<String> 	ctgList = new ArrayList<String>(),
						origList = new ArrayList<String>(), 
						resList = new ArrayList<String>(),
						genesList = new ArrayList<String>();

				for(Contig ctg:scaffolds[contig.head]){
					ctgList.add(ctg.getName());
					if(ctg.oriRep.size()>0)
						for(JapsaFeature ori:ctg.oriRep)
							origList.add(ori.getID());
					for (JapsaFeature feature:ctg.genes)
						genesList.add(feature.toString());
					for (JapsaFeature feature:ctg.resistanceGenes)
						resList.add(feature.toString());
				}
				float streamData=tpoint/1000000;
				pw.print(">");
				for(String ctg:ctgList)
					pw.printf("%s\t", ctg);

				pw.printf("\n>%.2fMpb\t%d genes\t", streamData, genesList.size());

				for(String ori:origList)
					pw.printf("+%s", ori);

				for(String genes:genesList)
					pw.print(" \n\t"+genes);
				pw.println("");

				for(String res:resList)
					pw.print(" \n\t"+res);
				pw.println("");

				pw.close();

			}
		}


	}	
	
	
//	// To check if this contig is likely a repeat or a singleton. If FALSE: able to be used as a marker.
//	public static boolean isRepeat(Contig ctg){
//		//for the case when no coverage information of contigs is found
//		if(estimatedCov == 1.0 && ctg.getCoverage() == 1.0){
//			if(ctg.length() > maxRepeatLength)
//				return false;
//			else
//				return true;
//		}
//
//		if (ctg.length() < minContigLength || ctg.getCoverage() < .3 * estimatedCov) return true;
//		else if (ctg.length() > maxRepeatLength || ctg.getCoverage() < 1.3 * estimatedCov) 
//			return false; 
//		else if (ctg.getCoverage() > 1.5 * estimatedCov)
//			return true;
//		else{
//			for(ContigBridge bridge:getListOfBridgesFromContig(ctg)){
//				Contig other = bridge.firstContig.getIndex()==ctg.getIndex()?bridge.secondContig:bridge.firstContig;
//				if(other.getIndex()==ctg.getIndex()) continue;
//				int dist=bridge.getTransVector().distance(bridge.firstContig, bridge.secondContig);
//				if( dist<0 && dist>-ctg.length()*.25){
//					if(other.length() > maxRepeatLength || other.getCoverage() < 1.3*estimatedCov)
//						return true;
//				}
//			}
//
//		}
//		if(ctg.length() < 2*minContigLength) // further filter: maybe not repeat but insignificant contig 
//			return true;
//		else 
//			return false;
//	}
	
	// To check if this contig is likely a repeat or a singleton. If FALSE: able to be used as a marker.
	// note: isMarker() = !isRepeat() in old version
	public static boolean isMarker(Contig ctg){
		//for the case when no coverage information of contigs is found
		if(estimatedCov == 0){
			if(ctg.length() > maxRepeatLength)
				return true;
			else
				return false;
		}else if(ctg.getCoverage() < .3*estimatedCov) {
			return false;
		}else
			return !ctg.getRepeatFlag();
	}

	abstract public void connectBridges();

	public int getNumberOfContigs(){
		return scfNum;
	}
	public int getNumberOfCirculars(){
		return cirNum;
	}

}