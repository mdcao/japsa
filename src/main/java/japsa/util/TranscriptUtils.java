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
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import pal.distance.DistanceMatrix;
import pal.misc.IdGroup;
import pal.misc.SimpleIdGroup;
import pal.tree.NeighborJoiningTree;
import pal.tree.NodeUtils;

public class TranscriptUtils {

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
	static double round2 = 100;
	
	static String[] nmes = "5_3:5_no3:no5_3:no5_no3".split(":");
	
	public static class CigarClusters {
		final double thresh;
		
		CigarClusters(double thresh){
			this.thresh = thresh;
		}

		DistanceMatrix getDistanceMatrix( PrintWriter pw){
			int len  = l.size();
			double[][] res = new double[len][];
			String[] labels = new String[len];
			for(int i=0; i<len; i++) {
				
				res[i] = new double[len];
				CigarCluster cc = l.get(i);
				labels[i] = cc.id;
				res[i][i] =0;
				for(int j=0; j<i; j++) {
					CigarCluster cc_j = l.get(j);
					double dist = 1-cc.similarity(cc_j);
					res[i][j] = dist;
					res[j][i] = dist;
				}
				
			}
			for(int i=0; i<len; i++) {
				pw.print(labels[i]+","+l.get(i).index+",");
				pw.println(getString(res[i]));
			}
		
			IdGroup group = new  SimpleIdGroup(labels);
			DistanceMatrix dm = new DistanceMatrix(res, group);
			return dm;
		}
		List<CigarCluster> l = new ArrayList<CigarCluster>();

	
		
		public String matchCluster(CigarCluster c1, int index, int source_index, int num_sources) {
		
			String clusterID="";
			double best_sc0=0;
			int best_index = -1;
			for (int i = 0; i < l.size(); i++) {	
				double sc = l.get(i).similarity(c1,index,   thresh);
				
					if(sc> best_sc0) {
						best_sc0 = sc;
						best_index = i;
					}
				
			}
			
			if (best_sc0 >= thresh) {
				CigarCluster clust = l.get(best_index);
				clust.merge(c1);
				clusterID = clust.id;
			//	System.err.println("merged");
			} else {
				CigarCluster newc = new CigarCluster("ID"+ l.size(), index,num_sources);
				newc.readCount[source_index]++;
				newc.merge(c1);
				clusterID = newc.id;
				l.add(newc);
				System.err.println("new cluster " + best_sc0 + " " + best_index+" "+newc.id+" "+index);
			//	System.err.println(l.size());
			}
			
			return clusterID;
		}

		public void getConsensus(Annotation annot, Sequence refseq,  PrintWriter exonP , SequenceOutputStream seqFasta, double[] depth, int num_sources) throws IOException{
			int[] first_last = new int[2];
			
			exonP.println("ID,index,start,end,"+getString("count", num_sources,true));
			for(int i=0; i<l.size(); i++) {
				CigarCluster cc = l.get(i);
				int[][] exons = cc.getExons( 0.3,10, depth);
				String id = cc.id;
				String read_count = getString(cc.readCount);
				StringBuffer descline = new StringBuffer();//cc.index+","+read_count);
				StringBuffer subseq= new StringBuffer();
				StringBuffer annotline = new StringBuffer();
				int transcript_len =0;
				for(int j=0; j<exons.length; j++) {
					int start = exons[j][0];
					int end = exons[j][1];
					exonP.println(id+","+cc.index+","+start+","+end+","+read_count);

					
					annotline.append(annot.calcORFOverlap(start, end, first_last, transcript_len));

					int len = end-start+1;
					descline.append(";");
					descline.append(start); descline.append("-"); descline.append(end); descline.append(","); descline.append(len);
					
					//descline.append("|");descline.append(annot.getInfo(first_last[0]));
					//descline.append("|");descline.append(annot.getInfo(first_last[1]));
					subseq.append(refseq.subSequence(start, end).toString());
					
					System.err.println(subseq.length());
					System.err.println(subseq);
					System.err.println("h");
					transcript_len += len;
					//seqline.append(subseq.toString());
				}
				Sequence subseq1 = new Sequence(refseq.alphabet(),subseq.toString().toCharArray(), id);
			//	subseq1.setName(id);
				descline.append(" "); descline.append(annotline);
				subseq1.setDesc(descline.toString());
				subseq1.writeFasta(seqFasta);
			//	seqFasta.println(idline.toString());
			//	seqFasta.println(seqline.toString());
			}
			
		}

	}

	public static class CigarCluster {
		final private int index;
		
		final String id;

		int start=0;
		int end=0;
		
		public CigarCluster(String id, int index, int num_sources) {
			this.id = id;
			this.index = index;
			this.readCount = new int[num_sources];
		}

		private SortedMap<Integer, Integer> map = new TreeMap<Integer, Integer>(); //coverate at high res
		
		public void clear(int source_index) {
			map.clear();
			map100.clear();
			Arrays.fill(readCount, 0);
			readCount[source_index]=1;
			start =0;
			end=0;
		}
		
		private SortedMap<Integer, Integer> map100 = new TreeMap<Integer, Integer>(); //coverage at low res (every 100bp)

		public void add(int round) {
			if(round<start) start =round;
			else if(round>end) end = round;
			int round1 = (int) Math.floor((double)round/round2);
			map.put(round, map.containsKey(round) ? map.get(round) + 1 : 1);
			map100.put(round1, map100.containsKey(round1) ? map100.get(round1) + 1 : 1);
		}

		public String toString() {
			return map.keySet().toString();
		}

		public Iterator<Integer> keys() {
			return map.keySet().iterator();
		}

		public Iterator<Integer> tailKeys(int st) {
			return map.tailMap(st).keySet().iterator();
		}

		public String summary(Integer[] positions) {
			StringBuffer sb = new StringBuffer(id+","+index+","+getString(this.readCount));
			for (int i = 0; i < positions.length; i++) {
				int i1 = positions[i];
				int v = map.containsKey(i1) ? map.get(i1) : 0;
				sb.append(",");
				sb.append(v);
			}
			return sb.toString();
		}
		
		double getDepth(Integer i) {
			return this.map.containsKey(i) ? (double) map.get(i) :  0.0;
		}
		int[][] exons;
	
		public int[][] getExons( double threshPerc, int numsteps, double[] depth) {
			if(exons!=null) return exons;
			List<Integer> start1 = new ArrayList<Integer>();
			List<Integer> end1 = new ArrayList<Integer>();
			double thresh = (double) readCountSum()*threshPerc;
			boolean in =false;
			Arrays.fill(depth, 0);
			for(int i=this.start; i<this.end; i++) {
				depth[i] = getDepth(i);
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
			}
			return exons;
		}
		private int readCountSum() {
			int sum=0;
			for(int i=0; i<this.readCount.length; i++) sum = sum+readCount[i];
			return sum;
		}

		int[] readCount;
		
		
		/** if its going to be less than thresh we return zero */
		public double similarity(CigarCluster c1,int index,  double thresh) {
			if(this.index !=index) return 0;
			int overlap = overlap(c1.start,c1.end,start, end);
		    if(overlap<0) return 0;
		    else{
		    	double union = union(c1.start, c1.end, start, end, overlap) ;
		    	if(overlap / union < thresh) return 0;
		    }
			double sim = this.similarity(map100, c1.map100);
			if(sim<thresh ) return 0;
			else sim = this.similarity(map, c1.map);
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
	    	double sim1 = this.similarity(map100, c1.map100);
	    	if(sim1<0.5) return sim1;
	    	else return  this.similarity(map, c1.map);  
		}
		
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
		}
		
		static int sum(Map<Integer, Integer> source){
			return source.values().stream() .reduce(0, Integer::sum);
		}
		
		public void merge(CigarCluster c1) {
			if(c1.start < start) start = c1.start;
			if(c1.end > end) end = c1.end;
			for(int i=0; i<this.readCount.length;i++) {
				readCount[i]+=c1.readCount[i];
			}
			int sum1 = merge(map, c1.map);
			int sum2 = merge(map100,c1.map100);
			if(sum1!=sum2){
				throw new RuntimeException("maps not concordant");
			}
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
		final File outfile, outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7;
		
		
		
		public IdentityProfile1(Sequence refSeq, File resDir, int num_sources, int genome_index, int round, boolean calculateCoExpression, double overlapThresh, int startThresh, int endThresh) throws IOException {
			//File readClusterFile = new File(outdir, "readclusters.txt.gz");
			this.round = (double) round;
			this.num_sources = num_sources;
			this.coRefPositions = new CigarCluster("reuseable",0,num_sources);
		 TranscriptUtils.round2 = 100.0/round;
			this.genome = refSeq;
			this.startThresh = startThresh; this.endThresh = endThresh;
			this.source_index = 0;		
			this.calculateCoExpression = calculateCoExpression;
			 outfile = new File(resDir,genome_index+ ".txt");
			 outfile1 = new File(resDir, genome_index+ "coref.txt");
			 outfile2 = new File(resDir, genome_index+"clusters.txt.gz");
			 outfile3 = new File(resDir,genome_index+ "readToCluster.txt.gz");
			 outfile4 = new File(resDir,genome_index+ "exons.txt.gz");
			 outfile5 = new File(resDir,genome_index+ "clusters.fa.gz");
			 outfile6 = new File(resDir,genome_index+ "tree.txt.gz");
			 outfile7 = new File(resDir,genome_index+ "dist.txt.gz");
			 readClusters = new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile3))));
				this.readClusters.println("readID,clusterID,index,source_index");//+clusterID+","+index+","+source_index);
			readClusters.println("readID,clusterId,type,source");
			
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

			refBase = 0;
			readBase = 0;
			// the number of bases from ref and read
			// following should be made more efficient
			Set<Integer> roundedPos = new HashSet<Integer>();
			for (int i = 0; i < refSeq.length(); i++) {
				roundedPos.add(round(i));
			}
			roundedPositions = roundedPos.toArray(new Integer[0]);
			this.depth = new double[roundedPos.size()];
			/*
			 * for(int i=0; i<roundedPositions.length; i++){ roundedPositions[i] =
			 * 1+i*(int)round; }
			 */
			codepth = new SparseRealMatrix[nmes.length];
			all_clusters =new CigarClusters(overlapThresh);
			if(this.calculateCoExpression) {
				for (int i = 0; i < this.codepth.length; i++) {
					codepth[i] = new OpenMapRealMatrix(roundedPositions.length, roundedPositions.length);
				}
			}
		}

		private int round(int pos) {
			int res = (int) Math.floor((double) pos / round);
			return res;
		}

		final int startThresh, endThresh;
		final  double round ;
		private final boolean calculateCoExpression;

		

		public void processRefPositions(int startPos, int distToEnd, String id) {
			int index = 0;
			if (startPos < startThresh)
				index = distToEnd < endThresh ? 0 : 1;
			else
				index = distToEnd < endThresh ? 2 : 3;
			Iterator<Integer> it = coRefPositions.keys();
			if(calculateCoExpression) {
				while (it.hasNext()) {
					Integer pos1 = it.next();
					Iterator<Integer> it2 = coRefPositions.tailKeys(pos1);
					while (it2.hasNext()) {
						Integer pos2 = it2.next();
						double value = codepth[index].getEntry(pos1, pos2);
						this.codepth[index].setEntry(pos1, pos2, value + 1);
					}
				}
			}
		//	coRefPositions.index = index;
			String clusterID = this.all_clusters.matchCluster(coRefPositions,index, this.source_index, this.num_sources); // this also clears current cluster
			this.readClusters.println(id+","+clusterID+","+index+","+source_index);
		}

		public void addRefPositions(int position) {
			coRefPositions.add(round(position));
		}

		public SparseRealMatrix[] codepth;
		private final CigarCluster coRefPositions;
		private CigarClusters all_clusters;
		private Integer[] roundedPositions;// , corefSum;
		private final double[] depth; // convenience matrix for depth
		public int[] match, mismatch, refClipped, baseDel, baseIns;
		public int numIns, numDel, readClipped, refBase, readBase;

		private final PrintWriter readClusters;
		private final Sequence genome;
		private final int num_sources;
		public int source_index; //source index
		public void updateSourceIndex(int i) {
			this.source_index = i;
		}
		
		public void print(PrintWriter pw, Sequence seq) {
			pw.println("pos,base,match,mismatch,refClipped,baseDel,baseIns");
			for (int i = 0; i < match.length; i++) {
				pw.println((i + 1) + "," + seq.charAt(i) + "," + match[i] + "," + mismatch[i] + "," + refClipped[i]
						+ "," + baseDel[i] + "," + baseIns[i]);
			}
			pw.flush();
		}

		public void printClusters(File outfile1) throws IOException {
			
				PrintWriter pw = new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1))));
				StringBuffer sb = new StringBuffer("pos,NA,"+getString("NA",num_sources,false));

				for (int i = 0; i < this.roundedPositions.length; i++) {
					sb.append(",");
					sb.append(roundedPositions[i] * round + 1);
				}
				pw.println(sb.toString());
				for (int i = 0; i < this.all_clusters.l.size(); i++) {
					pw.println(this.all_clusters.l.get(i).summary(this.roundedPositions));
				}
				pw.close();
		}

		

		public void printCoRef(File outfile1) throws IOException {
			if(calculateCoExpression) {
			for (int index = 0; index < this.codepth.length; index++) {
				SparseRealMatrix cod = this.codepth[index];
				if(cod==null) continue;
				File outfile1_ = new File(outfile1.getParentFile(),
						outfile1.getName() + "." +  nmes[index] + ".gz");
				PrintWriter pw = new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_))));
				int len = this.codepth[index].getRowDimension();

				for (int i = 0; i < len; i++) {
					int row = i; // nonZeroRows.get(i);
					pw.print(roundedPositions[row] * round + 1);
					for (int j = 0; j < len; j++) {
						int col = j;// nonZeroRows.get(j);
						int val = col >= row ? (int) cod.getEntry(row, col) : 0;
						// (int) this.codepth[index].getEntry(col, row);
						pw.print(",");
						pw.print(val);
					}
					pw.println();
				}
				pw.close();
			}
			}

		}
		
		public void getConsensus(Annotation annot) throws IOException {
						
				PrintWriter exonsP =  new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile4))));
				
				
				SequenceOutputStream seqFasta =  new SequenceOutputStream(new GZIPOutputStream(new FileOutputStream(outfile5)));
				this.all_clusters.getConsensus(annot, this.genome, exonsP, seqFasta, this.depth, this.num_sources);
				exonsP.close();
				seqFasta.close();
			
			
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
			pr1.print(pw, genome);
			pw.close();
			pr1.printCoRef(outfile1);
			pr1.printClusters(outfile2);
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
	public static void identity1(Sequence refSeq, Sequence readSeq, SAMRecord sam, IdentityProfile1 profile, int source_index) {


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
				profile.readClipped += length;
				break; // ignore hard clips
			case P:
				profile.readClipped += length;
				// pad is a kind of hard clipped ??
				break; // ignore pads
			case S:
				// advance on the reference
				profile.readClipped += length;
				readPos += length;
				break; // soft clip read bases
			case N:
				// System.err.println(length);
				refPos += length;
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
					profile.refClipped[refPos + i] += 1;
				}
				// profile.refClipped += length;
				break; // reference skip

			case D:// deletion
				refPos += length;
				profile.refBase += length;
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
					profile.baseDel[refPos + i] += 1;
				}
				profile.numDel++;
				break;

			case I:
				readPos += length;
				profile.readBase += length;

				profile.baseIns[refPos] += length;
				profile.numIns++;
				break;
			case M:
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
					profile.addRefPositions(refPos + i);
					if (refSeq.getBase(refPos + i) == readSeq.getBase(readPos + i))
						profile.match[refPos + i]++;
					else
						profile.mismatch[refPos + i]++;
				}
				profile.readBase += length;
				profile.refBase += length;

				readPos += length;
				refPos += length;
				break;

			case EQ:
				readPos += length;
				refPos += length;
				profile.addRefPositions(refPos);
				profile.readBase += length;
				profile.refBase += length;
				profile.match[refPos] += length;
				break;

			case X:
				readPos += length;
				refPos += length;

				profile.readBase += length;
				profile.refBase += length;
				profile.addRefPositions(refPos);
				profile.mismatch[refPos] += length;
				break;
			default:
				throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}// case
		} // for
		profile.processRefPositions(sam.getAlignmentStart(), refSeq.length() - sam.getAlignmentEnd(), id);
		// \return profile;

	}
	
	

}
