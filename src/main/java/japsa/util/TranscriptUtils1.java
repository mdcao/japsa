package japsa.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
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

import org.apache.commons.math3.Field;
import org.apache.commons.math3.FieldElement;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseFieldMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import pal.distance.DistanceMatrix;
import pal.misc.IdGroup;
import pal.misc.SimpleIdGroup;
import pal.tree.NeighborJoiningTree;
import pal.tree.NodeUtils;

public class TranscriptUtils1 {
	
	public static int break_thresh = 10;
	
	public static class IntegerField implements FieldElement<Integer>{
		Integer x;
		IntegerField(int x) {
			this.x = x;
		}
		public Integer getVal(){
			return x;
		}
		@Override
		public Integer add(Integer a) throws NullArgumentException {
			
			return  x+a;
		}

		@Override
		public Integer subtract(Integer a) throws NullArgumentException {
			return x - a;
		}

		@Override
		public Integer negate() {
			return -x;
		}

		@Override
		public Integer multiply(int n) {
			return x*n;
		}

		@Override
		public Integer multiply(Integer a) throws NullArgumentException {
			return x*a;
		}

		@Override
		public Integer divide(Integer a) throws NullArgumentException, MathArithmeticException {
			return null;
		}

		@Override
		public Integer reciprocal() throws MathArithmeticException {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public Field<Integer> getField() {
			// TODO Auto-generated method stub
			return null;
		}

		
	}
	public static class IField implements Field<IntegerField>{
		
			IntegerField zero = new IntegerField(0);
			IntegerField one = new IntegerField(1);
			@Override
			public IntegerField getZero() {
				// TODO Auto-generated method stub
				return zero;
			}

			@Override
			public IntegerField getOne() {
				// TODO Auto-generated method stub
				return one;
			}

			@Override
			public Class getRuntimeClass() {
				// TODO Auto-generated method stub
				return IntegerField.class;
			}

			
			
		
	}
	
	static int round(int pos, int round) {
		int res = (int) Math.floor((double) (pos+5) / (double) round);
		return res*round;
	}
	
	static class SparseVector{
		private SortedMap<Integer, Integer> m = new TreeMap<Integer, Integer>();
		int valsum=0;
		//double valsumByKey=0;
		
		public void addToEntry(Integer position, int i) {
			Integer val = m.get(position);
			valsum+=i;
			//valsumByKey+=position.doubleValue()* (double) i;
			m.put(position, val==null ? i : val + i);
		}
		public String toString(){
			return m.keySet().toString();
		}

		public List<Integer> keys() {
			List<Integer> l= new ArrayList<Integer>(this.m.keySet());
			Collections.sort(l);
			return l;
		}
		public Iterator<Integer> keyIt(){
			return m.keySet().iterator();
		}

		public Integer get(Integer val) {
			Integer res =  this.m.get(val);
			if(res==null) return 0;
			else return res;
		}
		
		
		public double similarity(SparseVector sv_B) {
			double intersection = 0;
			double Asize = m.size();//valsum;// m.size();
			double BnotA = 0;
			for (Iterator<Integer> it = sv_B.keyIt(); it.hasNext();) {
				if (m.containsKey(it.next())) {
					intersection+=1; //already in union
				} else {
					BnotA +=1;
				}
			}
			double union = Asize + BnotA;
			return ((double) intersection) / union;
		}

		public void clear() {
			m.clear();
			
		}
		public int getDepth(Integer i) {
			Integer val = m.get(i);
			return val==null ? 0: val;
		}
		public Iterator<Integer> tailKeys(Integer st) {
			return m.tailMap(st).keySet().iterator();
		}
		void merge(SparseVector source){
			for(Iterator<Integer> it = source.keyIt();it.hasNext();){
				Integer key = it.next();
				this.addToEntry(key, source.get(key));
			}
			
				//return this.valsum;
			}
		
		/*
		public static double similarity(Map<Integer, Integer> map, Map<Integer, Integer> m1) {
			int intersection = 0;
			int union = map.size();
			for (Iterator<Integer> it = m1.keySet().iterator(); it.hasNext();) {
				if (map.containsKey(it.next())) {
					intersection++;

				} else {
					union++;
				}
			}
			return ((double) intersection) / (double) union;
		}

		static int merge(Map<Integer, Integer> target, Map<Integer, Integer> source){
		//	int source_sum = sum(source);
		//	int target_sum = sum(target);
		//	int sum0 = source_sum+target_sum;
			Iterator<Entry<Integer, Integer>> it = source.entrySet().iterator();
		
			while (it.hasNext()) {
				Entry<Integer, Integer> entry = it.next();
				Integer key = entry.getKey();
				int curr = target.containsKey(key) ? target.get(key) : 0;
				int newv = curr + entry.getValue();
				target.put(key, newv);
				
			}
		
			return sum(target);
		}*/
	}
 static int[] printedLines =new int[] {0,0,0,0};
	public static int overlap(int st1, int end1, int st2, int end2){
		int minlen = Math.min(end1-st1, end2 - st2);
		int overlap = Math.min(end1-st2, end2 - st1);
		return Math.max(0, Math.min(minlen, overlap) +1);
	}
	public static int union(int st1, int end1, int st2, int end2, int overlap){
		int len1 = end1 - st1;
		int len2 = end2 - st2;
		if(overlap<=0) {
			return len1+len2;
		}else{
			int union =  Math.max(end1-st2, end2 - st1);
			int maxlen = Math.max(len1, len2);
			return Math.max(union, maxlen)+1 - overlap;
		}
	}
	
	static String getString(int[] c) {
		StringBuffer sb = new StringBuffer(c[0]+"");
		for(int i=1; i<c.length;i++) {
			sb.append(",");sb.append(c[i]);
		}
		return sb.toString();
	}
	static String getString(double[] c) {
		String formatstr = "%5.3g";
		StringBuffer sb = new StringBuffer(String.format(formatstr,c[0]).trim()+"");
		for(int i=1; i<c.length;i++) {
			sb.append(",");sb.append(String.format(formatstr, c[i]).trim());
		}
		return sb.toString();
	}
	 static String getString(String string, int num_sources2, boolean print_index) {
		 StringBuffer sb = new StringBuffer(string);
		 if(print_index) sb.append(0);
			for(int i=1; i<num_sources2;i++) {
				sb.append(",");sb.append(string);
				if(print_index)sb.append(i);
			}
			String res = sb.toString();
			return res;
	}
//	static double round2 = 100;
	
	static String[] nmes = "5_3:5_no3:no5_3:no5_no3".split(":");
	
	
	
	public static class CigarClusters {
		
		
		
		final double thresh;
		
		CigarClusters(double thresh){
			this.thresh = thresh;
		}
		public DistanceMatrix getDistanceMatrix( PrintWriter pw){
			CigarCluster[] l1 = this.l.values().toArray(new CigarCluster[0]);
			int len  = l1.length;
			double[][] res = new double[len][];
			String[] labels = new String[len];
			for(int i=0; i<len; i++) {
				
				res[i] = new double[len];
				CigarCluster cc = l1[i];
				labels[i] = cc.id;
				res[i][i] =0;
				for(int j=0; j<i; j++) {
					CigarCluster cc_j = l1[j];
					double dist = 1-cc.similarity(cc_j);
					res[i][j] = dist;
					res[j][i] = dist;
				}
				
			}
			for(int i=0; i<len; i++) {
				pw.print(labels[i]+","+l1[i].index+",");
				pw.println(getString(res[i]));
			}
		
			IdGroup group = new  SimpleIdGroup(labels);
			DistanceMatrix dm = new DistanceMatrix(res, group);
			return dm;
		}
		Map<CigarHash, CigarCluster> l = new HashMap<CigarHash, CigarCluster>();

	
		
		public String matchCluster(CigarCluster c1, int index, int source_index, int num_sources) {
		
			String clusterID="";
		
			
			if(!l.containsKey(c1.breaks)){
				CigarCluster newc = new CigarCluster(l.keySet().size()+"", index,num_sources);
				newc.setBreaks(c1.breaks);
				newc.addReadCount(source_index);
				newc.merge(c1);
				clusterID = newc.id;
				l.put(newc.breaks, newc);
				//System.err.println("new cluster " +" "+newc.id+" "+index);

			}else{
				CigarCluster clust = l.get(c1.breaks);
				clust.merge(c1);
				clusterID = clust.id;
			}
			return clusterID;
		}

		public void getConsensus(Annotation annot, Sequence refseq,  PrintWriter[] exonP ,
				PrintWriter[] transcriptsP, SequenceOutputStream[] seqFasta, PrintWriter[] clusterW, int[] depth, int num_sources) throws IOException{
			int[] first_last = new int[2];
			for(int i=0; i<exonP.length; i++){
				exonP[i].println("ID,index,start,end");
			//	transcriptsP[i].println("ID,index,start,end,startPos,endPos,totLen,countTotal,"+getString("count", num_sources,true));
				transcriptsP[i].println("ID,index,start,end,breaks,startPos,endPos,totLen,countTotal,"+getString("count", num_sources,true)
				+","+getString("depth", num_sources, true)+","+getString("errors", num_sources, true));
				clusterW[i].println("pos,type,"+getString("depth", num_sources, true)+","+getString("errors", num_sources, true));

			}
			
			//Collections.sort(l);
			int startPos = 0;
			for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();) {
				CigarCluster cc = it.next();
				startPos = printedLines[cc.index];
				int[][] exons = cc.getExons( 0.3,10, depth, clusterW[cc.index]);
				String id = cc.id;
				String read_count = getString(cc.readCount);
				StringBuffer descline = new StringBuffer();//cc.index+","+read_count);
				StringBuffer subseq= new StringBuffer();
				StringBuffer annotline = new StringBuffer();
				int transcript_len =0;
				int endPos = printedLines[cc.index];
				//transcriptsP[cc.index].println(cc.id+","+cc.index+","+cc.start+","+cc.end+","+startPos+","+endPos+","+cc.totLen+","+cc.readCountSum+","+read_count);
				transcriptsP[cc.index].println(cc.id+","+cc.index+","+cc.start+","+cc.end+","+cc.breaks.toString()+","+startPos+","+endPos+","+cc.totLen+","+cc.readCountSum+","+read_count
						+","+cc.getTotDepthSt()+","+cc.getTotErrorSt());
				for(int j=0; j<exons.length; j++) {
					int start = exons[j][0];
					int end = exons[j][1];
					exonP[cc.index].println(id+","+cc.index+","+start+","+end+","+read_count);

					
					annotline.append(annot.calcORFOverlap(start, end, first_last, transcript_len));

					int len = end-start+1;
					descline.append(";");
					descline.append(start); descline.append("-"); descline.append(end); descline.append(","); descline.append(len);
					
					//descline.append("|");descline.append(annot.getInfo(first_last[0]));
					//descline.append("|");descline.append(annot.getInfo(first_last[1]));
					subseq.append(refseq.subSequence(start, end).toString());
					
					//System.err.println(subseq.length());
					//System.err.println(subseq);
					//System.err.println("h");
					transcript_len += len;
					//seqline.append(subseq.toString());
				}
				Sequence subseq1 = new Sequence(refseq.alphabet(),subseq.toString().toCharArray(), id);
			//	subseq1.setName(id);
				descline.append(" "); descline.append(annotline);
				subseq1.setDesc(descline.toString());
				subseq1.writeFasta(seqFasta[cc.index]);
			//	seqFasta.println(idline.toString());
			//	seqFasta.println(seqline.toString());
			}
			
		}

	}

	public static class CigarHash extends ArrayList<Integer>{
		public static int round = 100;
		public String toString(){
			StringBuffer sb = new StringBuffer();
			sb.append(get(0));
			for(int i=1; i<size(); i++){
				sb.append(";");
				sb.append(get(i));
			}
			return sb.toString();
		}
		@Override
		public boolean add(Integer i){
			if(this.contains(i)) throw new RuntimeException("error");
			return (super.add(i));
		}
		
		@Override
		public int hashCode(){
			int hash = this.stream() .reduce(0, Integer::sum);
			return hash;
		}
		
		void roundBreaks(){
	   		for(int i =0; i<size(); i++){
	   			set(i,  round(get(i), round));
	   		}
	   	}
		/*@Override
public boolean equals(Object o){
			CigarHash cc = (CigarHash) o;
		 return this.equals(cc.breaks);	
		}
	}*/
	}
	
	public static class CigarCluster  {
		final private int index;
		
		
		static int round2 = 100;
		
		
		final String id;

		int start=0;
		int end=0;
		int prev_position =0;
		int break_point_cluster = -1;		
		boolean first = true;
		
		/*@Override
		public int compareTo(Object o) {
			CigarCluster ic1 = (CigarCluster)o;
			if(ic1.readCountSum==readCountSum) return 0;
			else return ic1.readCountSum<readCountSum ? -1 : 1;
		}*/
		
		public void addReadCount(int source_index) {
			readCount[source_index]++;
			this.readCountSum++;
			
		}
		
		
   	CigarHash breaks = new CigarHash() ;
 
   	
   	

	
		
		public void setBreaks(List<Integer> breaks){
			this.breaks.addAll(breaks);

		}
		public CigarCluster(String id, int index, int num_sources) {
			this.id = id;
			this.index = index;
			this.readCount = new int[num_sources];
			this.maps = new SparseVector[num_sources];
			this.errors= new SparseVector[num_sources];
			for(int i=0; i<maps.length; i++){
				maps[i] = new SparseVector();
				errors[i] = new SparseVector();
			}
		}

		private SparseVector map = new SparseVector(); //coverate at high res
	
		final private SparseVector[] maps, errors;
		public void clear(int source_index) {
			map.clear();
			for(int i=0; i<maps.length; i++){
				maps[i].clear();
				errors[i].clear();
			}
			Arrays.fill(readCount, 0);
			readCount[source_index]=1;
			readCountSum=1;
			start =0;
			end=0;
			this.prev_position = 0;
			this.break_point_cluster = -1;
			first = false;
			this.breaks.clear();
		}
	

		public void add(int pos, int src_index, boolean match) {
			if(first){
				first = false; 
				start = pos;
				breaks.add(pos);
			} else{
				end = pos; 
				if(match && (pos-prev_position>break_thresh)){
					this.breaks.add(prev_position);
					this.breaks.add(pos);
				}
			}
			prev_position = pos;
			
			map.addToEntry(pos, 1);
		
			//int round1 = (int) Math.floor((double)pos/round2);
		//	map100.addToEntry(round(pos,round2),1);
		
			maps[src_index].addToEntry(pos, 1);
			
			if(!match){
				errors[src_index].addToEntry(pos, 1);
			}
		}

		public String toString() {
			return this.breaks.toString();
		}
/*
		public Iterator<Integer> keys() {
			return map100.keyIt();
		}
*/
		/*public Iterator<Integer> tailKeys(int st) {
			return map100.tailKeys(st);
		}*/

		
		
		int getDepth(Integer i) {
			return this.map.getDepth(i);//this.map.containsKey(i) ?  map.get(i) :  0;
		}
		
		String getDepthSt(Integer i) {
			StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<maps.length; src_index++){
				if(src_index>0)sb.append(",");
				sb.append(this.maps[src_index].getDepth(i));//.containsKey(i) ?  maps[src_index].get(i) :  0);
			}
			return sb.toString();
		}
		
		String getTotDepthSt() {
			StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<maps.length; src_index++){
				if(src_index>0)sb.append(",");
				sb.append(this.maps[src_index].valsum);
			}
			return sb.toString();
		}
		String getErrorSt(Integer i) {
			StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<errors.length; src_index++){
				if(src_index>0)sb.append(",");
				sb.append(this.errors[src_index].get(i));
			}
			return sb.toString();
		}
		String getTotErrorSt() {
			StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<errors.length; src_index++){
				if(src_index>0)sb.append(",");
				sb.append(this.errors[src_index].valsum);
			}
			return sb.toString();
		}
		int[][] exons;
		
		public int[][] getExons( double threshPerc, int numsteps, int[] depth, PrintWriter clusterW) {
			if(exons!=null) return exons;
			List<Integer> start1 = new ArrayList<Integer>();
			List<Integer> end1 = new ArrayList<Integer>();
			double thresh = (double) readCountSum*threshPerc;
			boolean in =false;
			Arrays.fill(depth, 0);
			numPos =0;
			totLen =0;
			boolean prev0 = start>1;
			boolean printPrev = false;
			for(int i=this.start; i<this.end; i++) {
				depth[i] = getDepth(i);
				if(depth[i]>0){
					if(prev0 && !printPrev){
						numPos++;
						//clusterW.println((i-1)+","+this.id+","+getDepthSt(i-1));
						clusterW.println((i-1)+","+this.id+","+getDepthSt(i-1)+","+getErrorSt(i-1));
						printedLines[index]++;
					}
					numPos++;
					clusterW.println(i+","+this.id+","+getDepthSt(i)+","+getErrorSt(i));
//					clusterW.println(i+","+this.id+","+getDepthSt(i));
					printedLines[index]++;
					prev0 = false;
					printPrev = true;
				}else if(!prev0){
					numPos++;
					clusterW.println(i+","+this.id+","+getDepthSt(i)+","+getErrorSt(i));
				//	clusterW.println(i+","+this.id+","+getDepthSt(i));
					printedLines[index]++;
					printPrev = true;
					prev0=true;
				}else{
					printPrev = false;
					prev0 = true;
				}
			}
			outer: for(int i=start; i<=end; i++) {
				double dep = depth[i];
				if(!in && dep>=thresh) {
					for(int j = 1; j<numsteps && i+j < depth.length; j++) {
						if(depth[i+j]<thresh) continue outer; // no longer jumping in
					}
					in = true; 
					start1.add(i);
				}
				if(in && dep<thresh) {
					for(int j = 1; j<numsteps && i+j < depth.length; j++) {
						if(depth[i+j]>=thresh) continue outer; // no longer jumping out
					}
					in  = false;
					end1.add(i-1);
				}
			}
			if(end1.size() < start1.size()) end1.add(end);
			 exons = new int[end1.size()][];
			for(int i=0; i<end1.size(); i++) {
				exons[i] = new int[] {start1.get(i), end1.get(i)};
				totLen += end1.get(i) - start1.get(i)+1;
			}
			return exons;
		}
	
		
		
		int[] readCount; int readCountSum;
		int numPos =-1;
		int totLen = -1;
		
		
		/** if its going to be less than thresh we return zero */
		public double similarity(CigarCluster c1,int index,  double thresh) {
			if(this.index !=index) return 0;
			int overlap = overlap(c1.start,c1.end,start, end);
		    if(overlap<0) return 0;
		    else{
		    	double union = (double) union(c1.start, c1.end, start, end, overlap) ;
		    	if(overlap / union < thresh) return 0;
		    }
		//	double sim = map100.similarity(c1.map100);//this.similarity(map100, c1.map100);
			//if(sim<thresh ) return 0;
			double sim = map.similarity( c1.map);
			//System.err.println(highRes+" "+sim);
			return sim;
		}
		
		
		
		public static double similarity(int[][] exons1, int[][] exons2 ){
			int overlap =0;
			int union =0;
			for(int i=0; i<exons1.length; i++){
				int st1 = exons1[i][0];
				int end1 = exons1[i][1];
				for(int j=0; j<exons2.length; j++){
					int st2  =  exons2[j][0];
					int end2 = exons2[j][1];
					int overl = overlap(st1, end1,st2, end2);
					int unio = union(st1, end1,st2, end2, overl);
					overlap+=overl;
					union+= unio;
				}
			}
			return (double) overlap/(double) union;
		}
		public double exonSimilarity(CigarCluster c1){
			return similarity(c1.exons, this.exons);
		}
		public double similarity(CigarCluster c1) {
			int overlap = overlap(c1.start, c1.end,start, end);
		    if(overlap<=0) return 0;
		    int union = union(c1.start, c1.end,start, end, overlap);
		    double sim = (double) overlap/(double) union;
	    	if(sim < 0.5) return 0;
	    	//double sim1 = map100.similarity(c1.map100);
	    	//if(sim1<0.5) return sim1;
	    	 return  map.similarity( c1.map);  
		}
		
		
	/*	static int sum(Map<Integer, Integer> source){
			return source.values().stream() .reduce(0, Integer::sum);
		}*/
		
		public void merge(CigarCluster c1) {
			if(c1.start < start) start = c1.start;
			if(c1.end > end) end = c1.end;
			for(int i=0; i<this.readCount.length;i++) {
				readCount[i]+=c1.readCount[i];
			}
			this.readCountSum+=c1.readCountSum;
			//System.err.println(this.breaks.toString()+"\t"+c1.index+" "+readCountSum);
					
			 map.merge( c1.map);
		//	int sum2 = map100.merge(c1.map100);
			for(int i=0; i<maps.length; i++){
				maps[i].merge( c1.maps[i]);
				errors[i].merge(c1.errors[i]);
			}
			//if(sum1!=sum2){
				//throw new RuntimeException("maps not concordant");
			//}
		}
	}

	public static class Annotation{
		List<String> genes= new ArrayList<String>();
		List<Integer> start = new ArrayList<Integer>();
		List<Integer> end = new ArrayList<Integer>();
		
		int[][] orfs; 
		//Name,Type,Minimum,Maximum,Length,Direction,gene
	//	5'UTR,5'UTR,1,265,265,forward,none

	public	Annotation(File f) throws IOException{
			
			BufferedReader br = new BufferedReader(new FileReader(f)) ;
			List<String> header = Arrays.asList(br.readLine().split(","));
			int gene_ind = header.indexOf("gene");
			int start_ind = header.indexOf("Minimum");
			int end_ind = header.indexOf("Maximum");
			String str = "";
			while((str=br.readLine())!=null){
				String[] line = str.split(",");
				String gene = line[gene_ind];
				if(!gene.equals("none")) {
					genes.add(gene);
					int st = Integer.parseInt(line[start_ind]);
					int en = Integer.parseInt(line[end_ind]);
					start.add(st);
					end.add(en);
				
				}
			}
			br.close();
			orfs = new int[start.size()][];
			overlap = new double[start.size()];
			orf_len = new int[start.size()];
			
			for(int i=0; i<orf_len.length; i++) {
				orf_len[i] = end.get(i) - start.get(i)+1;
				orfs[i] = new int[] {start.get(i),end.get(i)};
			}
			
		}
	/*public String getInfo(int i) {
			if(i<0 || overlap[i] <0) return "NA";
			else return this.genes.get(i)+","+this.start_offset[i]+","+this.end_offset[i];
		}*/
	public final double[] overlap;
	public final int[] orf_len;
	//public final int[] start_offset; //relative start position of ORF
//	public final int[] end_offset; //relative end position of ORF
		public String calcORFOverlap(int start, int end, int[] first_last, int transcript_len) {
			//could consider correcting here to keep in-fram
			int first = -1;
			int last = -1;
			StringBuffer sb = new StringBuffer();
			for(int i=0 ; i<orfs.length; i++) {
				int[] st_en = orfs[i];
			
				overlap[i] =  (double) overlap(start, end, st_en[0], st_en[1])/ (double)(st_en[1] - st_en[0] +1);
				if(overlap[i]>0) {
					if(first<0) first = i;
					last = i;
					sb.append(";");
					sb.append(this.genes.get(i));
					sb.append(",");
					sb.append(String.format( "%5.3g",overlap[i]).trim());
					sb.append(",");
					sb.append(st_en[0] - start + transcript_len); // how far into read ORF starts;
					sb.append(",");
					sb.append(st_en[1] - start+transcript_len); // how far into read ORF ends
				}
			}
			first_last[0] =first;
			first_last[1] = last;
			return sb.toString();
		}
	}
	
	public static class IdentityProfile1 {
		final File outfile, outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9;
		//OpenMapRealMatrix
		SparseFieldMatrix sm;
		
		
		
		
		public IdentityProfile1(Sequence refSeq, File resDir, List<Integer[]> positions, int num_sources, int genome_index,  boolean calculateCoExpression, double overlapThresh, int startThresh, int endThresh) throws IOException {
			//File readClusterFile = new File(outdir, "readclusters.txt.gz");
			if(positions.size()>0){
				int rowlen = refSeq.length();
				sm = new SparseFieldMatrix(new IField(), rowlen, rowlen);
			for(int i=0; i<positions.size(); i++){
				Integer[] pos = positions.get(i);
				sm.setEntry(pos[0], pos[1],  new IntegerField(pos[2]));
			}
			}
			
			//this.bin = (double) bin;
			this.num_sources = num_sources;
			this.coRefPositions = new CigarCluster("reuseable",0,num_sources);
	//	 TranscriptUtils.round2 = 100.0/round;
			this.genome = refSeq;
			this.startThresh = startThresh; this.endThresh = endThresh;
			this.source_index = 0;		
			this.calculateCoExpression = calculateCoExpression;
			 outfile = new File(resDir,genome_index+ ".txt");
			 outfile1 = new File(resDir, genome_index+ "coref.txt");
			 outfile2 = new File(resDir, genome_index+"clusters.txt");
			 outfile3 = new File(resDir,genome_index+ "readToCluster.txt.gz");
			 outfile4 = new File(resDir,genome_index+ "exons.txt");
			 outfile5 = new File(resDir,genome_index+ "clusters.fa");
			 outfile6 = new File(resDir,genome_index+ "tree.txt.gz");
			 outfile7 = new File(resDir,genome_index+ "dist.txt.gz");
			 outfile8 = new File(resDir,genome_index+ "transcripts.txt");
			 outfile9 = new File(resDir,genome_index+ "breakpoints.txt");

			 readClusters = new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile3))));
				//this.readClusters.println("readID,clusterID,index,source_index");//+clusterID+","+index+","+source_index);
			readClusters.println("readID,clusterId,type,source,break_cluster");
			
			/*
			readClipped = 0;
			numDel = 0;
			numIns = 0;
			match = new int[refSeq.length()];
			mismatch = new int[refSeq.length()];
			refClipped = new int[refSeq.length()];
			baseDel = new int[refSeq.length()];
			baseIns = new int[refSeq.length()];

			Arrays.fill(match, 0);
			Arrays.fill(mismatch, 0);
			Arrays.fill(baseDel, 0);
			Arrays.fill(baseIns, 0);
			Arrays.fill(refClipped, 0);
			 */
			refBase = 0;
			readBase = 0;
			// the number of bases from ref and read
			// following should be made more efficient
			int seqlen = refSeq.length();
			//Set<Integer> roundedPos = new HashSet<Integer>();
			//for (int i = 0; i < refSeq.length(); i++) {
			//	roundedPos.add(i);
			//}
		//	roundedPositions = roundedPos.toArray(new Integer[0]);
			//this.depth = new int[roundedPos.size()];
		
			all_clusters =new CigarClusters(overlapThresh);
			this.breakpoints = new SparseRealMatrix[this.num_sources];
			this.breakSt = new SparseVector[this.num_sources];
			this.breakEnd = new SparseVector[this.num_sources];
			for(int i=0; i<breakpoints.length; i++){
				this.breakSt[i] = new SparseVector();
				this.breakEnd[i] = new SparseVector();
				breakpoints[i] = new OpenMapRealMatrix(seqlen, seqlen);
			}
			
		}

		

		final int startThresh, endThresh;
		
		//final  double bin ;
		private final boolean calculateCoExpression;

		final String[] nmes =  "5_3:5_no3:no5_3:no5_no3".split(":");

		//this after all postions in a read processed
		
		static double break_round = 10.0;
		public void processRefPositions(int startPos, int distToEnd, String id, boolean cluster_reads) {
			List<Integer> breaks  = coRefPositions.breaks;
			breaks.add(coRefPositions.end);
			
			int maxg = 0;
			int maxg_ind = -1;
			for(int i=1; i<breaks.size()-1; i+=2){
				int gap = breaks.get(i+1)-breaks.get(i);
				if(gap>maxg){
					maxg = gap;
					maxg_ind = i;
				}
			}
			
			if(maxg>100){
				Integer prev_position = breaks.get(maxg_ind);
				Integer position = breaks.get(maxg_ind+1);
				if(this.sm!=null){
					 coRefPositions.break_point_cluster = ((IntegerField)this.sm.getEntry(prev_position, position)).getVal();
					
				}
				this.breakpoints[this.source_index].addToEntry(prev_position, position, 1);
				this.breakSt[this.source_index].addToEntry(prev_position, 1);
				this.breakEnd[this.source_index].addToEntry(position,  1);
			}
			
			int index = 0;
			if (startPos < startThresh)
				index = distToEnd < endThresh ? 0 : 1;
			else
				index = distToEnd < endThresh ? 2 : 3;
		//	Iterator<Integer> it = coRefPositions.keys();
			coRefPositions.breaks.roundBreaks();
		//	coRefPositions.index = index;
			String clusterID;
			if(cluster_reads) clusterID = this.all_clusters.matchCluster(coRefPositions,index, this.source_index, this.num_sources); // this also clears current cluster
			else clusterID = "NA";
			this.readClusters.println(id+","+clusterID+","+index+","+source_index+","+coRefPositions.break_point_cluster);
		}
		
		public void addRefPositions(int position, boolean match) {
			coRefPositions.add(position, this.source_index, match);
			
		}

		
		private final CigarCluster coRefPositions;
		private CigarClusters all_clusters;
		final public SparseRealMatrix[] breakpoints;
		final public SparseVector[] breakSt, breakEnd;
		//private Integer[] roundedPositions;// , corefSum;
		private int[] depth; // convenience matrix for depth
		//public int[] match, mismatch, refClipped, baseDel, baseIns;
		//public int numIns, numDel, readClipped,
		public int refBase, readBase;

		private final PrintWriter readClusters;
		private final Sequence genome;
		private final int num_sources;
		public int source_index; //source index
		public void updateSourceIndex(int i) {
			this.source_index = i;
		}
		
		/*public void print(PrintWriter pw, Sequence seq) {
			pw.println("pos,base,match,mismatch,refClipped,baseDel,baseIns");
			for (int i = 0; i < match.length; i++) {
				pw.println((i + 1) + "," + seq.charAt(i) + "," + match[i] + "," + mismatch[i] + "," + refClipped[i]
						+ "," + baseDel[i] + "," + baseIns[i]);
			}
			pw.flush();
		}*/

		/*public void printClusters(File outfile1) throws IOException {
			
				PrintWriter pw = new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1))));
				StringBuffer sb = new StringBuffer("pos,NA,"+getString("NA",num_sources,false));
				for (int i = 0; i < this.roundedPositions.length; i++) {
					sb.append(",");
					sb.append(roundedPositions[i] * round + 1);
				}
				pw.println(sb.toString());
				for (int i = 0; i < this.all_clusters.l.size(); i++) {
					this.all_clusters.l.get(i).print(pw);
				//	pw.println(this.all_clusters.l.get(i).summary(this.roundedPositions));
				}
				pw.close();
		}*/

		private void printBreakPoints(File outfile1)  throws IOException {
			// TODO Auto-generated method stub
			for(int i=0; i<this.breakpoints.length; i++){
			File outfile1_ = new File(outfile1.getParentFile(),
					outfile1.getName() + "." + i + ".gz");
			printMatrix(this.breakpoints[i],this.breakSt[i], this.breakEnd[i],  outfile1_);
			}
		}
		 void printMatrix(SparseRealMatrix cod, SparseVector breakSt2, SparseVector breakEnd2, File outfile1_) throws IOException{
			
			PrintWriter pw = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_))));
			StringBuffer secondLine = new StringBuffer();
			List<Integer> rows  = breakSt2.keys();
			List<Integer> cols =  breakEnd2.keys();
			for(Iterator<Integer> it = cols.iterator(); it.hasNext();){
				pw.print(",");
				Integer val = it.next();
				pw.print(val);
				secondLine.append(",");
				secondLine.append(breakEnd2.get(val));
			}
			pw.println();
			pw.println(secondLine.toString());
			//Set<Integer> cols = breakEnd2
			for (Iterator<Integer> it =rows.iterator(); it.hasNext();) {
				Integer row = it.next(); // nonZeroRows.get(i);
				pw.print(row);
				pw.print(",");
				pw.print(breakSt2.get(row));
				for(Iterator<Integer> it1 = cols.iterator(); it1.hasNext();){
					Integer col = it1.next();// nonZeroRows.get(j);
					int val =  (int) cod.getEntry(row, col);// : 0;
					pw.print(",");
					pw.print(val);
				}
				pw.println();
			}
			pw.close();
		}

		
		
		public void getConsensus(Annotation annot) throws IOException {
				int num_types = this.nmes.length;
				this.depth = new int[this.genome.length()];
				PrintWriter[] exonsP = new PrintWriter[num_types]; 
				PrintWriter[] transcriptsP = new PrintWriter[num_types]; 
				SequenceOutputStream[] seqFasta = new SequenceOutputStream[num_types];
				PrintWriter[] clusterW = new PrintWriter[num_types];
				for(int i=0; i<num_types; i++){
				
						exonsP[i] =new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile4+"."+nmes[i]+".gz"))));
				
				seqFasta [i]=  new SequenceOutputStream(new GZIPOutputStream(new FileOutputStream(outfile5+"."+nmes[i]+".gz")));
				transcriptsP [i]=  new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile8+"."+nmes[i]+".gz"))));
				clusterW[i]= new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile2+"."+nmes[i]+".gz"))));
				}
				this.all_clusters.getConsensus(annot, this.genome, exonsP, transcriptsP, seqFasta,clusterW, this.depth, this.num_sources);
				for(int i=0; i<num_types; i++){
				clusterW[i].close();
				exonsP[i].close();
				seqFasta[i].close();
				transcriptsP[i].close();
				}
			
			
		}
		public void printTree() throws IOException{
			PrintWriter treeP =  new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile6))));
			PrintWriter distP =  new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile7))));
			System.err.println("calculating distance matrix..");
			DistanceMatrix dm = this.all_clusters.getDistanceMatrix( distP);
			System.err.println("..done");
			NeighborJoiningTree tree = new NeighborJoiningTree(dm);
			//treeP.print(tree.toString());
			NodeUtils.printNH(treeP, tree.getRoot(), true, false, 0, true);
			treeP.close();
			
			
			distP.close();
			
		}

		public void finalise()  throws IOException{
			this.readClusters.close();
			IdentityProfile1 pr1 = this;
			PrintWriter pw = new PrintWriter(new FileWriter(outfile));
		//	pr1.print(pw, genome);
			pw.close();
		//	pr1.printCoRef(outfile1);
			pr1.printBreakPoints(outfile9);
		//	pr1.printClusters(outfile2);
			System.out.println("========================= TOTAL ============================");
			
		}

		

		public void newRead(int source_index2) {
			this.coRefPositions.clear(source_index2);
		}

	}

	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * 
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public static void identity1(Sequence refSeq, Sequence readSeq, SAMRecord sam, IdentityProfile1 profile, int source_index, boolean cluster_reads) {


		int readPos = 0;// start from 0
		int refPos = sam.getAlignmentStart() - 1;// convert to 0-based index
		//String id = sam.getHeader().getId();
		String id = sam.getReadName();
		profile.newRead(source_index);
		for (final CigarElement e : sam.getCigar().getCigarElements()) {
			final int length = e.getLength();

			switch (e.getOperator()) {
			case H:
				// nothing todo
			//	profile.readClipped += length;
				break; // ignore hard clips
			case P:
			//	profile.readClipped += length;
				// pad is a kind of hard clipped ??
				break; // ignore pads
			case S:
				// advance on the reference
				//profile.readClipped += length;
				readPos += length;
				break; // soft clip read bases
			case N:
				// System.err.println(length);
				refPos += length;
				/*for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
					profile.refClipped[refPos + i] += 1;
				}*/
				// profile.refClipped += length;
				break; // reference skip

			case D:// deletion
				refPos += length;
				profile.refBase += length;
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
				//	profile.baseDel[refPos + i] += 1;
					profile.addRefPositions(refPos + i, false);
				}
				//profile.numDel++;
				break;

			case I:
				readPos += length;
				profile.readBase += length;

			//	profile.baseIns[refPos] += length;
			//	profile.numIns++;
				break;
			case M:
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
					
					if (refSeq.getBase(refPos + i) == readSeq.getBase(readPos + i)){
				//		profile.match[refPos + i]++;
						profile.addRefPositions(refPos + i, true);
					}
					else{
					//	profile.mismatch[refPos + i]++;
						profile.addRefPositions(refPos + i, false);
					}
				}
				profile.readBase += length;
				profile.refBase += length;

				readPos += length;
				refPos += length;
				break;

			case EQ:
				readPos += length;
				refPos += length;
				for (int i = 0; i < length && (refPos + i) < refSeq.length(); i++) {
					profile.addRefPositions(refPos+i, true);
				}
				profile.readBase += length;
				profile.refBase += length;
				//profile.match[refPos] += length;
				break;

			case X:
				readPos += length;
				refPos += length;

				profile.readBase += length;
				profile.refBase += length;
				for (int i = 0; i < length && (refPos + i) < refSeq.length(); i++) {
					profile.addRefPositions(refPos+i, false);
				}
				//profile.addRefPositions(refPos);
			//	profile.mismatch[refPos] += length;
				break;
			default:
				throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}// case
		} // for
		profile.processRefPositions(sam.getAlignmentStart(), refSeq.length() - sam.getAlignmentEnd(), id, cluster_reads);
		// \return profile;

	}
	
	

}