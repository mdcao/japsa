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
 * 07/09/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.np;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.collect.ImmutableMap;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.bio.phylo.NCBITree;
import japsa.bio.phylo.Slug;
import japsa.seq.SequenceOutputStream;
import japsa.tools.bio.np.ReferenceDB;
import japsa.tools.seq.CachedFastqWriter;
import japsa.tools.seq.CachedOutput;
import japsa.tools.seq.CachedSequenceOutputStream;
import japsa.tools.seq.SequenceUtils;
import japsa.tools.seq.SparseVector;
import japsa.tools.seq.SparseVectorCollection;
import japsa.util.DoubleArray;
import japsa.util.HTSUtilities;
import pal.misc.Identifier;
import pal.tree.Node;

/**
 * @author Minh Duc Cao, Son Hoang Nguyen, Lachlan Coin
 *
 */
public class RealtimeSpeciesTyping {
	public static double start_thresh = 0.0;
	//int int total_reads =0;
	public static boolean writeExcludeFile = false;
	public static double targetOverlap =0.95;
	public static double 	 epsilon = 0.001;
	public static double removeLikelihoodThresh=0.0; // 
	public static int EMiterations = 5000;
//	public static boolean exclude = false; //whether to use exclude file to exclude reads
//	public static boolean writeConsensus = false;
	public static boolean reestimate = true;
	public static double overlap_thresh =0.5;
	public static int mincount=2; // this for designing regions for MSA
	public static int minlength =500; // min length for MSA
	//public static boolean useBases = true;
	public static double pseudo = 0.0;
	public static double maxOverlap = .1; // maximum overlap with excl region to remove read
	public static double base_=2; // lower numbers spread out the probability distribution, however it should be base10
	public static String bases_covered = "bases_covered";
	public static String fraction_covered="fraction_covered";
	
	private static final Logger LOG = LoggerFactory.getLogger(RealtimeSpeciesTyping.class);
	
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
	
	static double interpolate(int cov0, int cov, double perc0, double perc, double d) {
		// TODO Auto-generated method stub
		return (double) cov0 + ((d-perc0)/ (perc-perc0)) * ((double)cov-(double) cov0);
	}
	
	Comparator<Interval> compar_start = new Comparator<Interval>(){
		@Override
		public int compare(Interval o1, Interval o2) {
			return Integer.compare(o1.start,o2.start);
		}
		
	};
	Comparator<Interval> compar_len= new Comparator<Interval>(){
		@Override
		public int compare(Interval o1, Interval o2) {
			return -1*Integer.compare(o1.length(),o2.length());
		}
	};
	public static boolean JSON = false;
	public static double ALPHA=0.05;
	public static int MIN_READS_COUNT=0;
	public static Pattern writeSep = null;
//	public static boolean writeUnmapped = false;
	private RealtimeSpeciesTyper typer;
	private OutputStream outputStream_ = null;
	final File[] outputFile;
	//private BufferedReader indexBufferedReader;
	private HashSet<String> filterSet = new HashSet<String>();
	/**
	 * Minimum quality of alignment
	 */
	private double minQual = 0;
	private boolean twoDOnly = false;
	public CachedOutput[] fqw_unmapped = null;
	CachedOutput fqw_filtered = null;
	//final public File[] unmapped_reads;
	File[] exclFile;
	File[] consensusFile;
	ReferenceDB refDB = null;
	final int[] currentReadCount;
//	Integer currentReadCount = 0;
	int[] currentReadAligned;// = 0;
	int[] currentBaseCount;// = 0L;
	 File[] fastqdir, outdir; // one for each sample
	public final 	File[] exclude_file_out;
	public final File[]  cov_file_out;
	public final	File[] consensus_file_out;
	public final PrintWriter[] reads_out;
//	File referenceFile;
		
	//seq ID to species name (from index ref file)
	
	
	public String[] samples() {
		return sampleID;
	}
	
	//to output binned sequences
	public static boolean OUTSEQ=false;
	List<Coverage[]> species2ReadList = new ArrayList<Coverage[]>();

	 public static int getBases(List<SAMRecord> sams) {
			int res =0;
			for(int i=0; i<sams.size(); i++){
				res +=baseCount(sams.get(i));
			}
			return res;
		}
		public static  int baseCount(SAMRecord sam) {
			return sam.getAlignmentEnd() - sam.getAlignmentStart()+1;
		}
	
	public static void filter(List<SAMRecord> sams, List<SAMRecord> out, int besti){
		int len = sams.size();
		SAMRecord first = sams.get(besti);
		out.add(first);
		for(int i=0; i<len; i++){
			if(i==besti) continue;
			SAMRecord sam= sams.get(i);
			int overlap = overlap(sam, first);
			if(overlap>10){
				continue;
			}
			if(sam.getReadNegativeStrandFlag()!=first.getReadNegativeStrandFlag()){
			  continue;
			}
			out.add(sam);
		}
	}
	public static int overlap(SAMRecord sam1, SAMRecord sam2){
		int st1 = sam1.getAlignmentStart(); int end1 = sam1.getAlignmentEnd(); int st2 = sam2.getAlignmentStart(); int end2  = sam2.getAlignmentEnd();
		st1 = sam1.getReadPositionAtReferencePosition(st1);
		 end1 = sam1.getReadPositionAtReferencePosition(end1);
		 st2 = sam2.getReadPositionAtReferencePosition(st2);
		 end2 = sam2.getReadPositionAtReferencePosition(end2);

		int minlen = Math.min(end1-st1, end2 - st2);
		int overlap = Math.min(end1-st2, end2 - st1);
		return Math.min(minlen, overlap);
	}
	
	public void getOutfiles(List<File> []res) {
		Iterator<Coverage[]> it =  species2ReadList.iterator();
		while(it.hasNext()){
			Coverage[] nxt = it.next();
			for(int i=0; i<res.length; i++){
			if(nxt[i].fqw!=null){
				nxt[i].fqw.getOutFile(res[i]);
			}
			}
		}
	
	}
	
	
	//** sorted based on coverage */
	public static class Interval implements Comparable{
		public int start,end, coverage;
		@Override
		public boolean equals(Object o){
			return start== ((Interval)o).start && end== ((Interval)o).end && coverage== ((Interval)o).coverage;  
		}
		public int length() {
			return end - start+1;
		}
		public String toString(){
			return start+"-"+end+","+bases()+","+coverage;
		}
		public Interval(){
			
		}
		public Interval(int start2, int end2, int i) {
		if(end2 <start2) throw new RuntimeException(start2+","+end2);
			this.start = start2; this.end = end2; this.coverage = i;
		}
		
		public Interval(String[] str) {
			this(str,2,3);
		}
		public Interval(String[] str, int start, int end) {

			this.start = Integer.parseInt(str[start]);
			this.end= Integer.parseInt(str[end]);

		}
		public int overlap(Interval interval){
			int start2 = interval.start; int end2 = interval.end;
			int overl = Math.min(end -start2,end2 -start);
			int minlen = Math.min(end2 - start2, end-start);
			return Math.min(overl, minlen);
		}
	
			
		//** smallest to largest coverage and largest span to smallest span
		@Override
		public int compareTo(Object o) {
			int res =  Integer.compare(coverage, ((Interval)o).coverage);
			if(res==0) {
				res = -1* Integer.compare(bases(), ((Interval)o).bases() );
			}
			return res;
		}

		public Integer bases() {
			int res =  end-start + 1;
			if(res<=0) throw new RuntimeException(" segment size of zero");
			return res;
		}
		public void merge(Interval nxt) {
			this.start = Math.min(start, nxt.start);
			this.end = Math.max(end, nxt.end);
			this.coverage = Math.min(nxt.coverage, coverage);
			
		}
		public double overlap(SAMRecord sam) {
			int start2 = sam.getAlignmentStart(); int end2 =sam.getAlignmentEnd();
			int overl = Math.min(end -start2,end2 -start);
			int minlen = Math.min(end2 - start2, end-start);
			return Math.min(overl, minlen);
		}
	}
	
	 
	/** this class represents the coverage of each species
	 * currently only for single chrom species, otherwise each chrom is considered its own species
	 *  */
	class Coverage{
		
		Coverage(String species,  Node node, File fastqdir, boolean writeSep1, boolean hierarchical, boolean fasta, boolean separateIntoContigs,
				ZipFile coverage_file
			){
			
			this.species = species;
			this.node = node;
			fqw = null;
			if(writeSep1){
				if(hierarchical && node!=null){
					Node parent = node.getParent();
					StringBuffer sd = new StringBuffer();
					while(!parent.isRoot()){
						sd.insert(0,"/"+Slug.toSlug(parent.getIdentifier().getName().replace("+-","").trim(),"_"));
						parent = parent.getParent();
					}
					fastqdir = new File(fastqdir.getAbsolutePath()+sd.toString());
				}
				
				fqw = fasta ? new CachedSequenceOutputStream(fastqdir, species, true, RealtimeSpeciesTyping.alignedOnly) : 
				new CachedFastqWriter(fastqdir, species, separateIntoContigs, RealtimeSpeciesTyping.alignedOnly);

			}
			if(coverage_file!=null){
				try{
			//	Enumeration entries = 	coverage_file.entries();
				ZipEntry entry = coverage_file.getEntry(species);
				if(entry!=null){
				add = false;
				
				BufferedReader br = new BufferedReader(new InputStreamReader(coverage_file.getInputStream(entry)));
				String st = "";
				while((st = br.readLine())!=null){
					String[] str = st.split("\t");
					String refname = str[1];
					int index = this.contig_names.indexOf(refname);
					if(index<0){
						index = contig_names.size();
						contig_names.add(refname);
						this.coverages.add(new TreeMap<Integer, Interval>());
					}
					SortedMap<Integer, Interval> coverage = this.coverages.get(index);
					if(!str[5].equals("-")){
					String[] str1 = str[5].split(";");
					for(int i=0; i<str1.length; i++){
						String[] stra = str1[i].split(",");
						Interval interval = new Interval(stra[0].split("-"),0,1);
						interval.coverage = Integer.parseInt(stra[2]);
						coverage.put(interval.start, interval);
					}
					}
				}
				br.close();
				}
				
				}catch(IOException exc){
					exc.printStackTrace();
				}
			}
			
		}
		
	    boolean add = true;
		
		Node node; // this is the node in the tree.  We dont initialise this until there are at least one read
		String species;
		
	//	int len;
		Interval[] newi = new Interval[3];
		SortedMap<Integer, Integer> mapq = new TreeMap<Integer, Integer>();
		SortedMap<Integer, Integer> mapAlign = new TreeMap<Integer, Integer>();
		SortedMap<Integer, Integer> mapLen = new TreeMap<Integer, Integer>();

		int readCount=0; int baseCount=0;
		// we keep the intervals non nested
		List<SortedMap<Integer,Interval>> coverages = new ArrayList<SortedMap<Integer, Interval>>();
		List<String> contig_names = new ArrayList<String>();
		//SortedMap<Integer, Interval> coverage = new TreeMap<Integer, Interval>(); //based on start
		
		
		
		
		//gets all regions with coverage greater than equal to thresh
				public List<Interval> getRegions(Integer thresh, int index, int merge_thresh) {
					Iterator<Interval> it = this.coverages.get(index).values().iterator();
					List<Interval> coverage = new ArrayList<Interval>();

					while(it.hasNext()){
						Interval nxt = it.next();
						if(nxt.coverage>=thresh){
							coverage.add(nxt);
						}
					}
					if(coverage.size()<=1) return coverage;
					Collections.sort(coverage, compar_start);
					List<Interval> cov1 = new ArrayList<Interval>();
					Interval curr = coverage.get(0);
					cov1.add(curr);
					for(int i=1; i<coverage.size(); i++){
						Interval nxt = coverage.get(i);
						if(curr.overlap(nxt)>=merge_thresh){
							curr.merge(nxt);
						}else{
							curr = nxt;
							cov1.add(curr);
						}
					}
					return cov1;
				}
		
		CachedOutput fqw= null;
		boolean lock = false;
		
		public void updateNodeAndParents(int src_index, double readCount){
			if(node!=null){
				Number[] cnts = (Number[]) this.node.getIdentifier().getAttribute(NCBITree.count_tag);
				Number[] cnts1 = (Number[]) this.node.getIdentifier().getAttribute(NCBITree.count_tag1);
			/*	if(cnts==null){
					this.node.getIdentifier().setAttribute(NCBITree.count_tag, cnts = new Number[] {0} );
					this.node.getIdentifier().setAttribute(NCBITree.count_tag1, cnts1 = new Number[] {0} );

				}*/
				cnts[src_index]=readCount;
				cnts1[src_index] = readCount;
				Node parent = node.getParent();
				while(parent!=null ){
				//adds counts to the tree
					Number[] cnts1_p = (Number[]) parent.getIdentifier().getAttribute(NCBITree.count_tag);
					/*if(cnts1_p==null){
						parent.getIdentifier().setAttribute(NCBITree.count_tag, cnts1_p = new Number[] {0} );
						parent.getIdentifier().setAttribute(NCBITree.count_tag1,  new Number[] {0} );

					}*/
					cnts1_p[src_index] = cnts1_p[src_index].doubleValue() + readCount;
					parent = parent.getParent();
				}
			}
		}
		
		//covs= new TreeMap<Integer, Double>()
		public  synchronized  Double[]  medianReadCoverage(double[] percentiles, double[] vals, int index,
				SortedMap<Integer, Integer> covs , SortedMap<Integer, List<Interval>> intervalMap,
				SortedMap<Integer, Integer> segsMap
				){
			covs.clear();
			intervalMap.clear();
			segsMap.clear();
			SortedMap<Integer, Interval> coverage = this.coverages.get(index);
			int non_zero_bases=0;
			Integer len = refDB.species2Len.get(this.contig_names.get(index));// node==null ? null : (Integer) node.getIdentifier().getAttribute("length");
			while(lock){
				try{
					LOG.debug("waiting on lock 1");
				Thread.sleep(100);
				}catch(InterruptedException exc){
					exc.printStackTrace();
				}
			}
			lock=true;
			int zeroBases=0;
			int prev =1;
			for(Iterator<Interval> it = coverage.values().iterator(); it.hasNext();){
				Interval nxt = it.next();
				int gap = nxt.start - prev;
				zeroBases+=gap;
				prev = nxt.end;
			}
			if(len!=null && len> prev) zeroBases+=len -prev;
			else len = prev;
			List<Interval> l = new ArrayList<Interval> ();
			l.addAll(coverage.values());
			Collections.sort(l);
			lock = false;
			
			
			double cumul = zeroBases;
		    
		if(covs!=null){
			covs.put(0,  zeroBases);
		}
			for(int j=0; j<percentiles.length; j++){
				if(percentiles[j] <=(double) zeroBases/(double) len){
					vals[j]  =0; // more zero bases than the percentile
				}
			}
			int cov0 =0;
			double cumul0 = cumul;;
			
			for(int i=0; i<l.size(); i++){
				Interval li = l.get(i);
				int cov = li.coverage;
				non_zero_bases +=li.bases();
				cumul = cumul0 +li.bases(); //percentage of bases
				if(covs!=null) {
					Integer cnt  = covs.get(cov);
					covs.put(cov, (cnt==null ? 0 :  cnt)+li.bases());
					List<Interval> lis = intervalMap.get(cov);
					if(!intervalMap.containsKey(cov)){
						intervalMap.put(cov,lis = new ArrayList<Interval>());
					}
					lis.add(li);
					Integer segs = segsMap.get(cov);
					segsMap.put(cov, segs==null ? 1 :  segs+1);
					
				}
				double perc0 = cumul0 / (double) len;
				double perc = cumul / (double) len;
				//System.err.println(perc0+","+perc);
				for(int j=0; j<percentiles.length; j++){
					if(percentiles[j] >= perc0  && percentiles[j] <=perc){
						vals[j]  = interpolate(cov0, cov, perc0, perc, percentiles[j]);
					}
				}
				cov0 = cov;
				cumul0 = cumul;
			}
			if(node!=null){
				Identifier id = this.node.getIdentifier();
				id.setAttribute(bases_covered, non_zero_bases);
				id.setAttribute(fraction_covered, (double) non_zero_bases/(double)len);
			}
			return new Double[] {(double)non_zero_bases, (double) non_zero_bases/(double)len, (double) len};
		//	System.err.println("h");
		}
		
		
		

		
		//extra are the extra bits still needed to be analysed (they may overlap elsewhere)
		// tomod are the bits which can be added directly.
		public void split(Interval i_new,Interval i_existing, List<Interval> tomod, Stack<Interval> extra_left, Stack<Interval> extra_right){
			boolean isLeft = i_new.start < i_existing.start ; // is new to left of existing based on start position
			Interval left = isLeft ? i_new: i_existing;
			Interval right =isLeft ? i_existing: i_new;
			boolean nested = right.end <=left.end;
			Interval i1,i2,i3;
			if(nested){ //right nested inside left
				// this replaces left.start and right.start
				i1 = left.start == right.start ?  null : (new Interval(left.start, right.start-1, left.coverage)) ; //left
				i2 = (new Interval(right.start, right.end, left.coverage + right.coverage)); //middle
				i3 = left.end ==right.end ? null : (new Interval(right.end+1, left.end,   left.coverage)); //right
				tomod.add(i2);
				if(isLeft){
					extra_left.push(i1); extra_right.push(i2); 
				}else{
					tomod.add(i1); tomod.add(i2);
				}
			}else{
				//try{
					
					i1 = left.start == right.start ?  null :(new Interval(left.start, right.start-1, left.coverage)); //left
					i2 = (new Interval(right.start, left.end, left.coverage + right.coverage)); //middle
					i3=  left.end ==right.end ? null : (new Interval(left.end+1,right.end,   right.coverage)); //right
					
					tomod.add(i2);
					tomod.add(isLeft ? i3 : i1);
					if(isLeft){
						extra_right.push(i3);
					}else{
						extra_left.push(i1);
					}
				//}catch(Exception exc){
				//	exc.printStackTrace();
				//}
			}
			
		}
		
		boolean inclSuppl = true;
		/** includes supplementary reads */
		public void addRead( List<SAMRecord> sams, int src_index, String readnme){
			
		//	System.err.println(sams.size()+" "+getAlignedFrac(sams.get(0)));
			if(sams.size()==0) return;
			int src_index1 = (Integer) sams.get(0).getAttribute(SequenceUtils.src_tag);
			if(src_index1!=src_index) {
				throw new RuntimeException("!!");
			}
			int len = inclSuppl ? sams.size(): 1;
			for(int i=0; i<len; i++){
				SAMRecord sam= sams.get(i);
				//if(node==null) node = tree.getNode(species);
				
				if(consensus_list!=null && i==0){
					String refName = sam.getReferenceName();
					//double len1 = sam.getAlignmentEnd()-sam.getAlignmentStart()+1;
					SortedSet< Interval>lis = consensus_list[src_index].get(refName);
					//SortedSet< Interval> lis = lis_==null ? null : lis_.tailMap(Math.min(minCoverage, lis_.lastKey()));
					if(lis!=null){
						int minCov1 = Math.min(minCoverage, lis.last().coverage);
						
						for(Iterator<Interval> it1 = lis.iterator(); it1.hasNext();){
							
							Interval interval = it1.next();
							if(interval.coverage <minCov1) continue;
							double overlapR = interval.overlap(sam);///(double )len;
							int len2 = interval.length();
							if(//overlapR/len1>targetOverlap && 
									overlapR/len2 > targetOverlap){ // this overlaps the target interval more than 90%
								fqw.write(sam, this.species, interval);
								//System.err.println("overlaps excluded region "+speciesList.get(seq2Species.get(refName))+" "+sam.getAlignmentStart()+" "+overlapR);
								//continue outer;
							}
						}
					}
				}else if(consensus_list==null && fqw!=null){	
					//sam.setReadName(sam.getReadName()+" "+sam.getReferenceName()+":"+sam.getAlignmentStart()+"-"+sam.getAlignmentEnd()+" "+sam.getMappingQuality());
					if(alignedOnly || i==0)  fqw.write(sam, this.species);
					
				}
		/*	if(RealtimeSpeciesTypingCmd.flipName){
				this.addInterval(sam.getReferenceName(),sam.getReadPositionAtReferencePosition(sam.getAlignmentStart()), 
							sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd()));
			}else{*/
			double[] values = new double[2] ;
			double[] values1 = new double[6];
			this.addInterval(sam.getReferenceName(),sam.getAlignmentStart(), sam.getAlignmentEnd(), values);
//			}
			//System.err.println(this.species);
			int q = sam.getMappingQuality();
			
			byte[] b = sam.getBaseQualities();
			double q1 = Double.NaN;
			if(b.length>0){
				double sump = 0;
				for(int ij=0; ij<b.length; ij++){
					sump+= Math.pow(10, -b[ij]/10.0);
				}
				sump = sump/(double) b.length;
				 q1 = -10.0*Math.log10(sump);
			}
			
			int st = sam.getReadPositionAtReferencePosition(sam.getAlignmentStart());
			int end = sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd());
			int alignL = (end-st+1);
			int readLength = (int) Math.round((double) sam.getReadLength()/10.0) * 10;
			int alignF = 10*(int) Math.round(10*(double)alignL/(double)sam.getReadLength());
			Integer cnt = mapq.get(q);
			Integer cnt1 = mapLen.get(readLength);
			Integer cnt2 = mapAlign.get(alignF);
			this.mapq.put(q, cnt==null ? 1 : cnt+1);
			this.mapLen.put(readLength, cnt1==null ? 1: cnt1+1);
			this.mapAlign.put(alignF, cnt2==null ? 1: cnt2+1);
		//	List<AlignmentBlock> blocks  = sam.getAlignmentBlocks();
			int baseCount_ = RealtimeSpeciesTyping.baseCount(sam);
			this.baseCount += baseCount_;
			this.readCount++;
			if(!add){
				HTSUtilities.identity(sam, values1);
				PrintWriter pw = reads_out[src_index];
				//vals1_str = "read_gap_first\tread_gaps\tread_gap_last\tref_gaps\talignment_blocks\talignment_blocks_rel"
				//pw.println("readnme\tspecies_name\tcontig_name\tquality\tphred\treadLength\talignFraction\tbaseCount\tmeanCov\tvarCov+\t"+vals1_str);
				pw.println(readnme+"\t"+this.species+"\t"+sam.getReferenceName()+"\t"+q+"\t"+String.format("%5.3g", q1).trim()
						+"\t"+readLength+"\t"+alignF+"\t"+baseCount_+getStr(values)+getStr(values1));
			}
			
			}
		}
	 
		private String getStr(double[] values) {
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<values.length; i++){
				sb.append("\t");
				sb.append(String.format("%5.3g",values[i]).trim());
			}
			return sb.toString();
		}

		private double getAlignedFrac(SAMRecord sam) {
			double a = sam.getEnd() - sam.getStart()+1;
			double b = sam.getReadLength();
			return a/b;
					
		}
		Interval firstOverlap(Interval i, int start_pos, int index){
			SortedMap<Integer, Interval> coverage = this.coverages.get(index);
			for(Iterator<Interval> it = coverage.tailMap(start_pos).values().iterator(); it.hasNext();){
				Interval nxt = it.next();
				if(i.overlap(nxt)>0) return nxt;
			}
			return null;
		}
		public String toString(){
			return this.species+"_"+coverages;//+":"+coverages..values().toString();
		}
		//we assuming existing segments non overlapping
		public void addInterval(String refname, int start, int end,  double[] value){
			
			int index = this.contig_names.indexOf(refname);
			if(index<0){
				index = contig_names.size();
				contig_names.add(refname);
				this.coverages.add(new TreeMap<Integer, Interval>());
			}
			SortedMap<Integer, Interval> coverage = this.coverages.get(index);
			
			if(add){
				while(lock){
					try{
						LOG.debug("waiting on lock");
						Thread.sleep(100);
					}catch(InterruptedException exc){
						exc.printStackTrace();
					}
				}
				lock=true;
			}
			List<Interval> toadd = new ArrayList<Interval>();
			Stack<Interval> extra_left = new Stack<Interval>(); // extra are like the leftover after merging
			Stack<Interval> extra_right = new Stack<Interval>();
			extra_right.push(new Interval(start, end,1));
			
//			System.err.println("add "+species+" "+i2 + this.coverage.size());
			int start_pos =0;
			while(extra_right.size()>0){ // extra left should all be non overlapping segments as we starting search from left
				Interval i_new = extra_right.pop();
				Interval i_existing = firstOverlap(i_new, start_pos, index);
				if(i_existing!=null){
					split(i_new, i_existing, toadd, extra_left, extra_right);
					start_pos = i_existing.end+1; 
				}else{
					toadd.add(i_new); // no overlaps, can just add segment
				}
			}
			toadd.addAll(extra_left); 
			// we dont need to remove original entries because we will replace them in this step
			if(add){
				for(Iterator<Interval> it = toadd.iterator(); it.hasNext();){
					Interval i1 = it.next();
					if(i1!=null){
						coverage.put(i1.start, i1);
					}
				}
				lock = false;
			}else{
				double sum=0; double sumSq=0;
				double len=0;
				inner: for(int i=0; i<toadd.size(); i++){
					Interval li = toadd.get(i);
					if(li==null) continue inner;
					int cov = li.coverage;
					
					int non_zero_bases =li.bases();
					sum+=non_zero_bases*cov;
					sumSq+=non_zero_bases*Math.pow(cov,2);
					len+=non_zero_bases;
					
				}
				value[0] = sum/len;
				value[1] = sumSq/len - Math.pow(sum/len,2); //E[X²] - (E[X])²

			}
			//System.err.println("species: "+this.toString());
			
		}
		public void close(){
			if(fqw!=null) this.fqw.close(refDB.species2Len);
		}
		public double readCount() {
			
			return readCount;
		}
		public double baseCount(){
			return baseCount;
		}
		
	}
	
	public static  File[] newF(File[] outdir, String nme, boolean mkDir){
		File[] f = new File[outdir.length];
		for(int i=0; i<f.length; i++){
			f[i] = nme==null ? null : new File(outdir[i], nme);
			if(mkDir)f[i].mkdir();
		}
		return f;
	}
	public RealtimeSpeciesTyping(File[] outdir, ReferenceDB refDB, 
			File[] exclFile,File[] consensus_file_in,  File[] cov_file_out, File[] reads_file_out, boolean writeUnmapped, String outputFile) throws IOException{
		this.outdir = outdir;
		this.sampleID =  getNames(outdir);
		this.num_sources = outdir.length;
		this.reads_out = reads_file_out==null ? null : new PrintWriter[reads_file_out.length];
		String vals1_str = "read_gap_first\tread_gaps\tread_gap_last\tref_gaps\talignment_blocks\talignment_blocks_rel";
		String header = "readnme\tspecies_name\tcontig_name\tquality\tphred\treadLength\talignFraction\tbaseCount\tmeanCov\tvarCov\t"+vals1_str;
		if(reads_out!=null){
			for(int i=0; i<reads_out.length; i++){
				reads_out[i] = new PrintWriter(new FileWriter(reads_file_out[i]));
				reads_out[i].println(header);
			}
		}
		this.currentBaseCount = new int[num_sources];
		this.currentReadAligned = new int[num_sources];
		this.currentReadCount = new int[num_sources];
			exclude_file_out =newF(outdir, "exclude.txt", false);
			consensus_file_out =newF(outdir, "consensus_regions.txt", false);
			
			this.outputFile = newF(outdir,outputFile, false);
			this.fastqdir= newF(outdir,"fastqs", writeSep!=null); 
			this.refDB = refDB;
		//	this.unmapped_reads  = newF(outdir, "unmapped", true);
		
	//	this.tree = tree;
		this.refDB = refDB;
		//this.indexFile = indexFile;
		this.exclFile = exclFile;//==null ? null : new File(exclFile); // this lists regions to exclude from count
		this.cov_file_out = cov_file_out;
		this.consensusFile =  consensus_file_in;
		//this.unmapped_reads = (new File(outdir, "unmapped")).getAbsolutePath();
		if(writeUnmapped){
			this.fqw_unmapped = new CachedFastqWriter[outdir.length] ;
			for(int i=0; i<fqw_unmapped.length ;i++){
				this.fqw_unmapped[i] = new CachedFastqWriter(outdir[i], "unmapped", false, false);
			}
		}
		this.fqw_filtered = null;//new CachedFastqWriter(outdir, "filtered");

	}
	//final NCBITree tree;
	
	//* referenceFile is to get the length map */
	public RealtimeSpeciesTyping(ReferenceDB refDB, File[] exclFile, File[] consensusFile, File[] cov_file,File[] reads_out,
			String outputFile, File[] outdir, File referenceFile, 
			boolean unmapped_reads, boolean keepNames, String[] src_names) throws IOException{
		this(outdir, refDB, exclFile,consensusFile, cov_file, reads_out, unmapped_reads, outputFile);
		//this.referenceFile = referenceFile;
	//	boolean useTaxaAsSlug=false;
		
	//	 addExtraNodesFromSpeciesIndex( tree,  indexFile, null);
		
		//this.indexBufferedReader = SequenceReader.openFile(indexFile);
		typer = new RealtimeSpeciesTyper(this,getNames(outdir));
		preTyping();
		if(reestimate) this.all_reads = new SparseVectorCollection(refDB.speciesList, keepNames, src_names);
	
	}

	private String[] getNames(File[] outdir2) {
		String[] res = new String[outdir2.length];
		for(int i=0; i<res.length; i++){
			res[i] = outdir2[i].getName();
		}
		return res;
	}
	final String[] sampleID;
	public RealtimeSpeciesTyping(ReferenceDB refDB, File[] exclFile,File[] consensusFile,  File[] cov_file, File[] reads_out,  File[] outdir, boolean unmapped_reads) throws IOException {
		this(outdir, refDB, exclFile, consensusFile, cov_file,reads_out,   unmapped_reads, "output.data");
		LOG.debug("string outputstream");
		
	//	this.indexBufferedReader = SequenceReader.openFile(indexFile);
	//	this.outputStream = outputStream;
		
		typer = new RealtimeSpeciesTyper(this, sampleID);
		preTyping();
	}

	/*public RealtimeSpeciesTyping(BufferedReader indexBufferedReader, String outputFile, File outdir) throws IOException {
		this(outdir);
		LOG.debug("bufferedreader string");

		//this.indexBufferedReader = indexBufferedReader;
		this.outputStream = SequenceOutputStream.makeOutputStream(outputFile);
		typer = new RealtimeSpeciesTyper(this, outputStream);
		preTyping();
	}

	public RealtimeSpeciesTyping(BufferedReader indexBufferedReader, OutputStream outputStream, File outdir) throws IOException {
		this(outdir);
		LOG.debug("bufferedreader outputstream");

	//	this.indexBufferedReader = indexBufferedReader;
		this.outputStream = outputStream;
		typer = new RealtimeSpeciesTyper(this, outputStream);
		preTyping();
	}*/

//	static class SpeciesCount implements Comparable<SpeciesCount> {
//		String species;
//		int count = 0;
//
//		SpeciesCount (String s){
//			species = s;
//		}
//
//		/* (non-Javadoc)
//		 * @see java.lang.Comparable#compareTo(java.lang.Object)
//		 */
//		@Override
//		public int compareTo(SpeciesCount o) {		
//			return o.count - count;
//		}
//
//	}
//	SequenceUtils.annotateWithGenomeLength(speciesIndex, seq2Species, species2Len);

public static boolean hierarchical = false;
public static boolean separateIntoContigs = true;//false; // whether to do separate file for each contig
public static boolean alignedOnly = true; // whether just to output the aligned component in the output fastq
public static int minCoverage = 2;
public static boolean fastaOutput = true;

Map<String, SortedSet< Interval>>[] exclude_list  = null;
Map<String, SortedSet<Interval>>[] consensus_list  = null;

public static List<String> speciesToIgnore = null;
	private void preTyping() throws IOException{
		
		exclude_list = new Map[exclFile.length];
		consensus_list= new Map[consensusFile.length];
			for(int i=0; i<consensus_list.length; i++){
				exclude_list[i]= getIntervals(exclFile[i], null);
				Map<String, List<String>> m =new HashMap<String, List<String>>();

				consensus_list[i] = 	getIntervals(consensusFile[i],m);
				System.err.println(m);
			}
	
		System.err.println(consensus_list);
		ZipFile[] zf= null;
		if(cov_file_out!=null){
			zf = new ZipFile[this.cov_file_out.length];
			for(int i=0; i<zf.length; i++){
				if(cov_file_out[i]!=null && cov_file_out[i].exists()){
					zf[i] = new ZipFile(cov_file_out[i]);
				}
			}
		}
		for(int i=0; i<refDB.speciesList.size(); i++){
			String sp = refDB.speciesList.get(i);
			boolean writeSep1 = writeSep !=null
					&& (speciesToIgnore==null || ! speciesToIgnore.contains(sp))
					&& writeSep.matcher(sp).find();
				Node n =refDB.getNode(sp);
				species2ReadList.add(getCoverage(sp,n, 	fastqdir, writeSep1, hierarchical, fastaOutput, separateIntoContigs, zf));			
		}
		if(zf!=null){
			for(int i=0; i<zf.length; i++){
			if(zf[i]!=null) 	zf[i].close();
			}
		}
	//	tree.makeTrees();
		
		//indexBufferedReader.close();
	
		LOG.info(refDB.seq2Species.size() + "   " + species2ReadList.size());
	//	speciesList.addAll(species2ReadList.keySet());
		
		//Write header				
	}

	private Coverage[] getCoverage(String sp, Node n, File[] fastqdir2, boolean writeSep1, boolean hierarchical2,
			boolean fastaOutput2, boolean separateIntoContigs2, ZipFile[] zf) {
		Coverage[] cov = new Coverage[fastqdir2.length];
		for(int i=0; i<cov.length; i++){
			cov[i] = new Coverage(sp,n, 	fastqdir[i], writeSep1, hierarchical, fastaOutput, separateIntoContigs, zf==null ? null : zf[i]);
		}
		return cov;
	}
	/**
	 * @param minQual the minQual to set
	 */
	public void setMinQual(double minQual) {
		this.minQual = minQual;
	}

	/**
	 * @param twoOnly the twoOnly to set
	 */
	public void setTwoOnly(boolean twoOnly) {
		this.twoDOnly = twoOnly;
	}

	/*public void typing(String bamFile, int readNumber, int timeNumber, List<String> species, boolean runAnalysis) throws IOException, InterruptedException {
		InputStream bamInputStream;

		if ("-".equals(bamFile))
			bamInputStream = System.in;
		else
			bamInputStream = new FileInputStream(bamFile);

		typing(bamInputStream, readNumber, timeNumber, species, runAnalysis);
	}*/

	/**\
	 * @param filter the species keywords list (separated by comma) to excluded
	 */
	public void setFilter(String filter) {
		if(filter !=null && !filter.isEmpty()){
			String[] toks = filter.split(";");
			for(String tok:toks)
				if(!tok.isEmpty()){
					filterSet.add(tok);
				}
		}
		
	}
	/*public void typing(InputStream bamInputStream, int readNumber, int timeNumber, List<String> species, boolean runAnalysis) throws IOException, InterruptedException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamInputStream));
		SAMRecordIterator samIter = samReader.iterator();
		typing(samIter, readNumber, timeNumber,  runAnalysis);
		samReader.close();
	}*/
	final int num_sources;
	public static boolean realtimeAnalysis = false;
	SparseVectorCollection all_reads = null ;
	public void typing(Iterator<SAMRecord>[] samIters, int readNumber, 
			int timeNumber,  boolean runAnalysis) throws IOException, InterruptedException{
		//if (readNumber <= 0)
		//	readNumber = 1;		
	
		
		typer.setReadPeriod(readNumber);
		typer.setTimePeriod(timeNumber * 1000);
		LOG.info("Species typing ready at " + new Date());
		String readName = "", refName = "";
	//	int src_index = -1;
		
		if(realtimeAnalysis && runAnalysis){
			Thread thread = new Thread(typer);
			LOG.info("starting RealtimeSpeciesTyper thread");
			thread.start();
			LOG.info("started  RealtimeSpeciesTyper thread");
		}
		HashSet<String> skipList = new HashSet<>();
		//String prevReadName = "";
		//String prevRefName= "";
		AllRecords records = new AllRecords(this.outdir); // for supplementary alignemnts to same reference
		for(int jk=0; jk<samIters.length; jk++){
			final int src_index = jk;
		Iterator<SAMRecord> samIter = samIters[jk];
		//Interval interval = new Interval();
		outer: while (samIter.hasNext()){
			try{
			SAMRecord sam = samIter.next();
		//	int src_index0 = (Integer) sam.getAttribute(SequenceUtils.src_tag);
			sam.setAttribute(SequenceUtils.src_tag, src_index);
			if(sam==null) {
				System.err.println("warning sam record is null");
				break;
			}
			if (this.twoDOnly && !sam.getReadName().contains("twodim")){
				continue;
			}
			if (!sam.getReadName().equals(readName) ){
				readName = sam.getReadName();
				//src_index = src_index0;
				synchronized(this){
					if(records.size()>0){
						int src_index1 = records.transferReads(species2ReadList, all_reads,currentReadCount, currentBaseCount);
						if(src_index1!=src_index) throw new RuntimeException("!!");
						//int  records.src_index;
						//if(src_index>=0){
						currentReadCount[src_index] ++;
						currentBaseCount[src_index] += sam.getReadLength();
					}
					//}
				}
			} 
			
			{
				if(start_thresh > 0){
			if(sam.getStart()>start_thresh * sam.getLengthOnReference()){
				continue ;
			}
			}
/*		System.err.println(len +" "+len1+" "+st+","+end+" "+st1+","+end1+  "   "+readLen+ " "+unc_start+" "+unc_end);
		System.err.println("h");*/
			}
		//	Integer src_index = (Integer)sam.getAttribute(SequenceUtils.src_tag);
			 refName = sam.getReferenceName();
			if (sam.getReadUnmappedFlag()){
				LOG.debug("failed unmapped check");
				if(fqw_unmapped!=null) {
					this.fqw_unmapped[src_index].write(sam,"unmapped");
				}
				continue;			
			}
			int len = sam.getAlignmentEnd() - sam.getAlignmentStart()+1;
			SortedSet<Interval>lis = exclude_list[src_index].get(refName);
			if(lis!=null){
				for(Iterator<Interval>it1 = lis.iterator(); it1.hasNext();){
					Interval int1 = it1.next();
					double overlapR = int1.overlap(sam);///(double )len;
					if(overlapR/len>overlap_thresh ){
						System.err.println("overlaps excluded region "+refDB.speciesList.get(refDB.seq2Species.get(refName))+" "+sam.getAlignmentStart()+" "+overlapR);
						continue outer;
					}
				}
			}
			//System.err.println(sam.getReadName()+"\t"+sam.getReferenceName()+"\t"+records.size()+"\t"+sam.isSecondaryOrSupplementary());
			int mq = sam.getMappingQuality();
			if (mq < this.minQual) {
				
				LOG.debug("failed minQual check "+mq);
				if(!sam.isSecondaryOrSupplementary()){
					if(fqw_unmapped!=null) this.fqw_unmapped[src_index].write(sam,"unmapped");
			
				}
				continue;
			}
			

			if(skipList.contains(readName)){
				
				LOG.debug("filter {}", readName);
				continue;
			}
			
			
			if(filterSet.contains(refDB.seq2Species.get(refName))){
			
				if(!sam.isSecondaryOrSupplementary()){
					skipList.add(readName);
					if(fqw_filtered!=null) this.fqw_filtered.write(sam, "filtered");
				}
				continue;
			}
			
		
			/* this probably not necessary  - also unclear what effect this has on coordinates
			if(sam.getReadNegativeStrandFlag()){
				//System.err.println("switching strand");
				SequenceUtils.flip(sam, false); // switch read but keep flag. This corrects minimaps correction when it maps
			}*/

			//if(readList.readList.size()==0 || !readList.readList.contains(readName))
				synchronized(this) {
					currentReadAligned[src_index] ++;	
					/*if(tree!=null){
						tree.addRead(this.seq2Species.get(sam.getReferenceName()));
					}*/
					String refname = sam.getReferenceName();
					//String readname = sam.getReadName()
					Integer specI = refDB.seq2Species.get(refname);
					if(specI==null){
						System.err.println("WARNING no species for "+refname+" "+this.sampleID[src_index]);
					}
					if(specI!=null) records.add(sam, specI);
				//	readList.addRead(sam);
					//readList.add(readName);
	
				}
		}catch(Exception exc){
				exc.printStackTrace();
			}
		}//while
		//
			records.transferReads(species2ReadList, all_reads, currentReadCount, currentBaseCount);
		records.close();
		

		//	for(int src_index1=0; src_index1 < num_sources; src_index1++){
				if(this.reads_out!=null) reads_out[src_index].close();
			
			if(fqw_unmapped!=null && fqw_unmapped[src_index]!=null) {
				this.fqw_unmapped[src_index].close();
			}
	}// samIters
		if(fqw_filtered!=null) this.fqw_filtered.close();
			//boolean hasnext = samIter.hasNext();
		typer.stopWaiting();//Tell typer to stop
		if(runAnalysis){ // This is the E-M step
			double[]v = new double[2];
			if( all_reads!=null){
				//double minv = 0.0;
				this.all_reads.maximisation(v, pseudo,RealtimeSpeciesTyping.EMiterations, null);
				//all_reads.check(all_reads.abundance);
				all_reads.getAbundances();
			}
			if(all_reads.keepNames){
				if(refDB.mostLikely==null){
					refDB.makeMostLikely(all_reads.num_sources);
				}
				this.all_reads.printMostLikely(refDB.speciesList, refDB.mostLikely);	
			}
			Integer[] nonZero = all_reads.nonZero(0.01);
			//String[] species = new String[nonZero.length];
		/*	if(species!=null){
			for(int i=0; i<nonZero.length; i++){
				species.add(refDB.speciesList.get(nonZero[i]));
			}
			}*/
			
			if( !realtimeAnalysis){
				typer.run();
			}
		}//if(runAnalysis)
		//samIter.close();
		
	}	

	private int findMinAbs(Double[] vals) {
		int min_i=0;
		double minv = Math.abs(vals[0]);
		if(vals.length>1){
		for(int j=1; j<vals.length; j++){
			double v = Math.abs(vals[j]);
			if(v<minv){
				minv = v;
				min_i= j;
			}
		}
		}
		return min_i;
	}
	private Map<String,SortedSet<Interval>> getIntervals(File excl, Map<String, List<String>> species2Seq ) {
		//boolean consensus = excl.getName().equals("consensus.txt")""
		Map<String, SortedSet<Interval>> m = new HashMap<String, SortedSet<Interval>>();
		Map<String, String>seq2Species = new HashMap<String, String>();
		boolean consensus = false;
		if(excl==null || !excl.exists()) return m;
		try{
			BufferedReader br = new BufferedReader(new FileReader(excl));
			String st = "";
			while((st = br.readLine())!=null){
				String[] str = st.split("\t");
				if(str.length>5) consensus = true;
			//	if(str.length==1) continue;
				String seq = str[1];
				SortedSet<Interval> li = m.get(seq);
				if(li==null) m.put(seq, li = new TreeSet<Interval>());
				
				seq2Species.put(seq, str[0]);
			//	List<String> seqs = species2Seq.get(str[0]);
				//if(seqs==null) species2Seq.put(str[0], seqs = new ArrayList<String>());
				if(consensus){ // this for consensus.txt
					String[] str1 = str[5].split(";");
					int coverage = Integer.parseInt(str[2]);
					if(coverage >=2){
					inner: for(int k=0; k<str1.length; k++){
						Interval in = new Interval(str1[k].split(","),0,1);
						if(in.length() > minlength){
							in.coverage = coverage;
							for(Iterator<Interval> it1 = li.iterator(); it1.hasNext();){
								if(it1.next().overlap(in) > 0.90 * in.length()) { // we only want to consider new region which does not significantly overlap current region
								//	System.err.println(li.get(j)+ " vs "+in);
									continue inner;
								}
							}
							li.add(in);
						}
					}
					}
				}else{ // this for exclude.txt
					li.add(new Interval(str));
				}
				
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
		Iterator<String> it = (new HashSet<String>(m.keySet())).iterator();
		while(it.hasNext()){
			String st = it.next();
			if(m.get(st).size()==0){
				m.remove(st);seq2Species.remove(st);
			}
		}
		Iterator<String> it1 = seq2Species.keySet().iterator();
		if(species2Seq!=null){
			while(it1.hasNext()){
				String st = it1.next();
				String st1 = seq2Species.get(st);
				List<String> li = species2Seq.get(st1);
				if(li==null) species2Seq.put(st1, li = new ArrayList<String>());
				li.add(st);
			}
		}
		return m;
	}


	public class RealtimeSpeciesTyper extends RealtimeAnalysis {
		MultinomialCI rengine;
		RealtimeSpeciesTyping typing;
		public SequenceOutputStream[] countsOS = null;
		File[] krakenResults; //kraken formatted results
	//	final String[] sampleID;
		//File coverageOutput;
		File[] outdir; 
		
		public int numSources(){
			return num_sources;
		}
		public void initOutput() throws IOException{
			countsOS = new SequenceOutputStream[num_sources];

			for(int i=0; i<num_sources; i++){
				String outf = outputFile[i].getAbsolutePath();
			countsOS[i] = SequenceOutputStream.makeOutputStream(outf);
			if(!JSON)
				countsOS[i].print("sampleID\ttime\tstep\treads\tbases\tspecies\tprob\terr\ttAligned\tsAligned\tbases_covered\tfraction_covered\tlength_best_contig\tcoverage_percentiles\tmapQ\tlength\talignFrac\tprop_to_most_cov_contig\thighest_cov_contig\tE_M_estimate\tlog_EM_estimate\n");
			}
		}
		
		public RealtimeSpeciesTyper(RealtimeSpeciesTyping t,  String[] sampleID) throws IOException {
			typing = t;
			this.outdir = t.outdir;
			krakenResults = newF(t.outdir,"results.krkn", false);
			
			rengine = new MultinomialCI(ALPHA);
//this.sampleID = sampleID;
		}
		//	countsOS = new SequenceOutputStream(outputStream);
		
		double[] perc = new double[] {0.5};// 0.5, 0.75, 0.9, 0.95, 0.99}; // percentiles for printing median
		double[] percQ = new double[] {0.5};
		double[] percL = new double[] {0.5};
		DoubleArray countArray = new DoubleArray();
		List<String> medianArray = new ArrayList<String>();
		List<Double> logEM = new ArrayList<Double>();
		ArrayList<String> speciesArray = new ArrayList<String> ();
		double[][] results = null;
		int[] order = null;
		public boolean lock = false;
		Long step;
		 double  frac = 1.0;
		private void simpleAnalysisCurrent(int src_index)  {
		//	for(int src_index=0; src_index<num_sources; src_index++){
			
				 SparseVector abundance_k = all_reads.abundances[src_index];
			lock = true; // so that the datastructure doesnt change while we doing calculation 
			double[] vals = new double[perc.length];
			double[] valsQ = new double[percQ.length];
			double[] valsL = new double[percL.length];
			step = (lastTime - startTime)/1000;//convert to second

			//int sum = 0;
			int count_len = typing.refDB.speciesList.size();
			SparseVector count = new SparseVector();

			for (int i = 0; i < count_len;i++){
			//	String spec_name = typing.speciesList.get(i);
				Coverage cov = typing.species2ReadList.get(i)[src_index];
				double counti  =  cov.readCount();		
				if(counti>0) count.addToEntry(i, counti);
				
			}
			int sum = (int) count.valsum();
			int minCount = MIN_READS_COUNT>0?MIN_READS_COUNT:Math.max(1,sum/1000);

			
			countArray.clear();medianArray.clear();speciesArray.clear();

			SortedMap<Integer,Integer> covMap = null;
			SortedMap<Integer,Integer> segsMap = null;

			SortedMap<Integer,List<Interval>> intervalMap = null;
			OutputStreamWriter  coverage_out = null; 
			  ZipOutputStream outS = null;
			PrintWriter  regions_to_exclude = null; 
			PrintWriter regions_to_use= null;
			if(final_analysis){
				covMap = new TreeMap<Integer, Integer>();// this can capture the distribution of bases against depth
				intervalMap = new TreeMap<Integer, List<Interval>>();
				segsMap = new TreeMap<Integer, Integer>();// this can capture the distribution of bases against depth
				try{
					File cov_out = cov_file_out[src_index];
						
					File excl_out = (exclude_file_out[src_index]);
					if(!cov_out.exists()){
						
						OutputStream dest =  new FileOutputStream(cov_out);
						CheckedOutputStream checksum = new   CheckedOutputStream(dest, new Adler32());
				       outS = new   ZipOutputStream(new   BufferedOutputStream(checksum));
				        coverage_out  = new OutputStreamWriter(outS);
				         outS.setMethod(ZipOutputStream.DEFLATED);
						
						
				//coverage_out = new PrintWriter(new FileWriter(cov_out));
					}
					if(!excl_out.exists()){
				if(writeExcludeFile) regions_to_exclude = new PrintWriter(new FileWriter(exclude_file_out[src_index]));
					}
					if(!consensus_file_out[src_index].exists()){
				regions_to_use=  new PrintWriter(new FileWriter(consensus_file_out[src_index]));
					}
				}catch(IOException exc){
					exc.printStackTrace();
				}
			}
			for (Iterator<Integer> it = count.keyIt(); it.hasNext();){
				Integer i = it.next();
				double counti = count.get(i).doubleValue();

				double abund_i = abundance_k.get(i).doubleValue();
				try{
				if (counti >= minCount){
					countArray.add(counti);
					String spec_name = typing.refDB.speciesList.get(i);
					ZipEntry ent = new ZipEntry(spec_name);
					 outS.putNextEntry(ent);
					
					speciesArray.add(spec_name);
					LOG.info(step+" : " + spec_name+ " == " + counti);
					//if(count[i]>0){
						Coverage cov = typing.species2ReadList.get(i)[src_index];
						
						 //SortedMap<Integer, Double>  covHist = 
						double max =0; int maxj=0; double tot_bases=0;
						double max1=0; int maxj_1 =0; 
						for(int j =0; j<cov.contig_names.size(); j++){
							
						  Double[] stats = cov.medianReadCoverage(new double[0], new double[0], j,covMap, intervalMap, segsMap);
						  
						  if(covMap!=null && coverage_out!=null){
							  Entry<Integer,Integer> prev = null;
							  Integer thresh = null;
							  for(Iterator<Entry<Integer, Integer>> covs = covMap.entrySet().iterator();covs.hasNext(); ){
								  Entry<Integer, Integer> nxt = covs.next();
								 
								  if(thresh==null && prev!=null){
									  if(//nxt.getKey()>prev.getKey()+1 ||
											 nxt.getValue() > frac* prev.getValue()){
										  thresh =nxt.getKey();
									  }
								  }
								  prev = nxt;
								  
								 List< Interval> iv = intervalMap.get(nxt.getKey());
								  String iv_str = iv==null ?  "-": toString(iv);
								  coverage_out.write(spec_name+"\t"+cov.contig_names.get(j)+"\t"+nxt.getKey()+"\t"+nxt.getValue()
								  +"\t"+segsMap.get(nxt.getKey())+"\t"+iv_str+"\n"
										 );
							  }
							  ;
							  if(thresh!=null && regions_to_exclude!=null){
								 List<Interval> exclusions =  cov.getRegions(thresh,j, -100);
								 for(int jk=0; jk<exclusions.size(); jk++){
									 Interval iv = exclusions.get(jk);
									regions_to_exclude.println(spec_name+"\t"+cov.contig_names.get(j)+"\t"+ iv.start+"\t"+iv.end+"\t"+iv.coverage);
								 }
							  }
							 if(regions_to_use!=null){
								 //print out regions to use
							List<Integer> keys = new ArrayList<Integer>(covMap.keySet());
							Collections.sort(keys);
							for(int jk=keys.size()-1; jk>=0; jk--){
							  int keyj  =keys.get(jk);
							  if(keyj>=mincount){
								  List<Interval> intervals = cov.getRegions(keyj,j, -10);
								  if(intervals.size()>0){
									  Collections.sort(intervals, compar_len);
									  int maxlength = intervals.get(0).length();
									
									  if( maxlength>=minlength){
										 
										  int ik=0;
										  while(ik<intervals.size() && intervals.get(ik).length()>=minlength) ik++;
										  regions_to_use.print(spec_name+"\t"+cov.contig_names.get(j)+"\t"+keyj+"\t"+ik+"\t"+maxlength+"\t");
										  for(int ik1=0; ik1<ik; ik1++){
											  Interval iv = intervals.get(ik1);
											   regions_to_use.print(iv.start+","+iv.end);
											   if(ik1<ik-1) regions_to_use.print(";");
											  
										  }
										  regions_to_use.println();
									  }
								  }
							  }
							  
							}
							 }
						  }
						  tot_bases+= stats[0];
						  if(stats[0] >max){
							  max = stats[0];
							  maxj=j;
						  }
						  if(stats[1] > max1){
							  max1 = stats[1];
							  maxj_1 =j;
						  }
						}
						Double proportion = max/ tot_bases;
						String nme = cov.contig_names.get(maxj);
						String nme1 = cov.contig_names.get(maxj_1);
						Double[] stats =  cov.medianReadCoverage(perc, vals,maxj, covMap, intervalMap,segsMap);
						 String st = combine(perc,vals);
						 medianQ(percQ, valsQ,cov.mapq);
						 medianQ(percL, valsL, cov.mapLen);
						 String st1 = combine(percQ,valsQ);
						 String st2 = combine(percL,valsL);
						 medianQ(percL, valsL, cov.mapAlign);
						 String st3 = combine(percL,valsL);
						 
						if(abund_i>0)						 {
							cov.updateNodeAndParents(src_index, abund_i*currentReadCount[src_index]);
						}

						 String adj_res = all_reads==null ? "":""+String.format("%5.3g",abund_i).trim();
						 String adj_res_log =all_reads==null ? "": ""+String.format("%5.2g",Math.log10(abund_i)).trim();
						 medianArray.add(String.format("%5.3g", stats[0]).trim()+"\t"+String.format("%5.3g",stats[1]).trim()+"\t"+String.format("%5.3g",stats[2]).trim()
						 +"\t"+st+"\t"+st1+"\t"+st2+"\t"+st3
								 +"\t"+nme+":"+String.format("%5.3g",proportion).trim()+"\t"+nme1+":"+max1+"\t"+adj_res+"\t"+adj_res_log);
						 logEM.add(abund_i);
						 if(adj_res.equals("0")){
							 System.err.println("h");
						 }
				/*	if(count[i]>10 * vals[vals.length-1]){
						System.err.println(count[i]+" vs "+st);
						if(covHist!=null) System.err.println(covHist);
					}*/
					LOG.info("medians "+step+" : "+spec_name+" "+ st+ " vs count "+counti);
					//}
					coverage_out.flush();
				    outS.closeEntry();
				} // if
				}catch(IOException exc){
					exc.printStackTrace();
				}
			}		
			//if (countArray.size() > 10) return;
			countArray.add(1); medianArray.add("--");
			speciesArray.add("others");		

			lock =false;
		//	int l1 = countArray.size();
		//	int l2 = logEM.size();
			rengine.assignCount(countArray.toArray());
			rengine.eval();        
			//REXP tab  = rengine.eval("tab",true);  
			results=rengine.tab();
			this.order = rengine.rank();
			//order = rank(medianArray)
			try{
			if(coverage_out!=null) coverage_out.close();
			}catch(IOException exc){
				exc.printStackTrace();
			}
			if(regions_to_exclude!=null) regions_to_exclude.close();
			if(regions_to_use!=null) regions_to_use.close();
			//}
		}
		
		

		private String toString(List<Interval> iv) {
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<iv.size(); i++){
				sb.append((iv.get(i).toString()));
				if(i<iv.size()-1) sb.append(";");
			}
			return sb.toString();
		}
		@Override
		protected void writeFinalResults(int src_index) {
			//for(int src_index=0; src_index < num_sources; src_index++){
			try{
		//	if(typing.refDB.tree!=null) typing.refDB.tree.zeroCounts(0, 1);; // add arrays to nodes for recording counts, or reset to zero
		//	Iterator<Coverage> it = typing.species2ReadList.iterator();
			/*while(it.hasNext()){
				Coverage cov = it.next();
				if(cov.readCount()> 0){		
					cov.updateNodeAndParents();
				}
				
			}*/
			
			}catch(Exception exc){
				exc.printStackTrace();
			}
			//}
			
		}
		
		
		
		
		private void writeResults(double min_thresh , int src_index) throws IOException {
			Gson gson = new GsonBuilder().serializeNulls().create();
			List<JsonObject> data = new ArrayList<JsonObject>();
			if(countsOS==null) this.initOutput();
			for (int i1 = 0; i1 < order.length;i1++){
				int i = order[i1];
				if (results[i][0] <= min_thresh)
					continue;

				Double mid = (results[i][0] + results[i][1])/2;
				Double err = mid - results[i][0];
				if(!JSON) {
					String med_a = medianArray.get(i);
					String line = sampleID[src_index]+"\t"+lastTime + "\t" + step + "\t" + lastReadNumber + "\t" + typing.currentBaseCount 
							+ "\t" + speciesArray.get(i).replaceAll("_", " ") + "\t" + mid + "\t" + err + "\t" + typing.currentReadAligned + 
							"\t" + countArray.get(i)+"\t"+med_a;
					
					System.err.println(med_a);
					countsOS[src_index].print(line);
					countsOS[src_index].println();
				}
				else {
					JsonObject jo = new JsonObject();
					jo.addProperty("species", speciesArray.get(i).replaceAll("_", " "));
					jo.addProperty("step", step.toString());
					jo.addProperty("reads", lastReadNumber.toString());
					jo.addProperty("bases", typing.currentBaseCount.toString());
					jo.addProperty("prob", mid.toString());
					jo.addProperty("err", err.toString());
					jo.addProperty("tAligned", typing.currentReadAligned.toString());
					jo.addProperty("sAligned", Double.valueOf(countArray.get(i)).toString());
					data.add(jo);

				}
			}

			if(JSON) {
				countsOS[src_index].print(gson.toJson(ImmutableMap.of(
						"timestamp", lastTime,
						"data", data
				)));
				countsOS[src_index].println();

			}
			countsOS[src_index].flush();
			LOG.info(step+"  " + countArray.size());
		}

		private String combine(double[] perc2, double[] vals) {
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<perc2.length; i++){
				sb.append(String.format("%5.3g=%5.3g", new Double[] {perc2[i], vals[i]}).replaceAll(" ", "")+","); //perc2[i]+"="+vals[i]+" ");
			}
			return sb.toString().trim();
		}

		protected void close(){
			try{
				//rengine.end();
				for(int src_index=0; src_index<num_sources; src_index++){
				if(countsOS!=null) countsOS[src_index].close();
				}
			}catch (Exception e){
				e.printStackTrace();
			}
			
			Iterator<Coverage[]> it = this.typing.species2ReadList.iterator();
			while(it.hasNext()){
				Coverage[] cov = it.next();
				for(int i=0; i<cov.length; i++)cov[i].close();
			}
			//print out
			
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#analysis()
		 */
		@Override
		protected void analysis(int i){
			simpleAnalysisCurrent(i);
			try{
				
				this.writeResults( 0.000001, i);
			}catch (IOException e){
				LOG.warn(e.getMessage());
			}
		}
		
		public boolean final_analysis = false;// flag if final analysis
		@Override
		protected void lastAnalysis(int i){
			this.final_analysis=true;
			simpleAnalysisCurrent(i);
			
			try{
				
				this.writeResults(-1, i);//to write all results
			}catch (IOException e){
				LOG.warn(e.getMessage());
			}
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#getCurrentRead()
		 */
		@Override
		protected int getCurrentRead() {
			int cnt =0;
			for(int i=0; i<currentReadCount.length; i++){
				cnt+= typing.currentReadCount[i];
			}
			
			return cnt;
//			return typing.currentReadCount.as
		}

		
	}



	

	

}
