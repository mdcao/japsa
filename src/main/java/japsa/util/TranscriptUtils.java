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
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Sequence;

public class TranscriptUtils {

	public static class CigarClusters {
		final double thresh;
		
		CigarClusters(double thresh){
			this.thresh = thresh;
		}

		List<CigarCluster> l = new ArrayList<CigarCluster>();

		public String matchCluster(CigarCluster c1) {
			double[] sc = new double[l.size()];
			int min_index = 0;
			double best_score = 0;
			String clusterID="";
			for (int i = 0; i < l.size(); i++) {
				sc[i] = l.get(i).similarity(c1);
				if (sc[i] > best_score) {
					min_index = i;
					best_score = sc[i];
				}
			}
			if (min_index < sc.length && sc[min_index] >= thresh) {
				CigarCluster clust = l.get(min_index);
				clust.merge(c1);
				clusterID = clust.id;
			} else {

				CigarCluster newc = new CigarCluster("ID" + l.size());
				newc.merge(c1);
				clusterID = newc.id;
				l.add(newc);
				System.err.println("new cluster " + best_score + " " + min_index);
				System.err.println(l.size());

			}
			c1.map.clear();
			return clusterID;
		}

		public void getConsensus(Annotation annot, Sequence refseq, Integer[] positions, PrintWriter exonP , PrintWriter seqFasta) {
			int[] first_last = new int[2];
			for(int i=0; i<l.size(); i++) {
				CigarCluster cc = l.get(i);
				int[][] exons = cc.getExons( positions,0.5);
				String id = cc.id;
				StringBuffer idline = new StringBuffer(">"+id);
				StringBuffer seqline = new StringBuffer();
				for(int j=0; j<exons.length; j++) {
					int start = exons[j][0];
					int end = exons[j][1];
					exonP.println(id+","+start+","+end);

					
					annot.calcORFOverlap(start, end, first_last);

					int len = end-start+1;
					idline.append(";");
					idline.append(start); idline.append("-"); idline.append(end); idline.append(","); idline.append(len);
					
					idline.append("|");idline.append(annot.getInfo(first_last[0]));
					idline.append("|");idline.append(annot.getInfo(first_last[1]));
					
				
					Sequence subseq  = refseq.subSequence(start, end);
				
					seqline.append(subseq.toString());
				}
				seqFasta.println(idline.toString());
				seqFasta.println(seqline.toString());
			}
			
		}

	}

	public static class CigarCluster {

		
		
		

		
		
		final String id;

		public CigarCluster(String id) {
			this.id = id;
		}

		private SortedMap<Integer, Integer> map = new TreeMap<Integer, Integer>();

		public void add(int round) {
			map.put(round, map.containsKey(round) ? map.get(round) + 1 : 1);
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
			StringBuffer sb = new StringBuffer(id+"_"+readCount);
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
		public int[][] getExons( Integer[] positions, double threshPerc) {
			List<Integer> start = new ArrayList<Integer>();
			List<Integer> end = new ArrayList<Integer>();
			double thresh = (double) readCount*threshPerc;
			boolean in =false;
			for(int i=0; i<positions.length; i++) {
				double dep = getDepth(positions[i]);
				if(!in && dep>=thresh) {
					in = true; 
					start.add(positions[i]);
				}
				if(in && dep<thresh) {
					in  = false;
					end.add(positions[i-1]);
				}
			}
			if(end.size() < start.size()) end.add(positions[positions.length-1]);
			int[][] exons = new int[end.size()][];
			for(int i=0; i<end.size(); i++) {
				exons[i] = new int[] {start.get(i), end.get(i)};
			}
			return exons;
		}
		int readCount =1;
		public double similarity(CigarCluster c1) {
			// if(true) return 1.0;
			int intersection = 0;
			int union = map.size();
			for (Iterator<Integer> it = c1.map.keySet().iterator(); it.hasNext();) {
				if (map.containsKey(it.next())) {
					intersection++;

				} else {
					union++;
				}
			}
			return ((double) intersection) / (double) union;
			// Set<Integer> union = Stream.concat(setA.stream(),
			// setB.stream()).collect(Collectors.toSet());
			// Set<Integer> intersect =
			// setA.stream().filter(setB::contains).collect(Collectors.toSet());
			// return ((double)intersect.size())/((double) union.size());
		}

		public void merge(CigarCluster c1) {
			this.readCount++;
			Iterator<Entry<Integer, Integer>> it = c1.map.entrySet().iterator();
			while (it.hasNext()) {
				Entry<Integer, Integer> entry = it.next();
				Integer key = entry.getKey();
				int curr = map.containsKey(key) ? map.get(key) : 0;
				map.put(key, curr + entry.getValue());
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
			overlap = new int[start.size()];
			orf_len = new int[start.size()];
			for(int i=0; i<orf_len.length; i++) {
				orf_len[i] = end.get(i) - start.get(i)+1;
				orfs[i] = new int[] {start.get(i),end.get(i)};
			}
			
		}
	public String getInfo(int i) {
			if(i<0 || overlap[i] <0) return "NA";
			else return this.genes.get(i)+","+this.start_offset[i]+","+this.end_offset[i];
		}
	public int[] overlap;
	public int[] orf_len;
	public int[] start_offset; //relative start position of ORF
	public int[] end_offset; //relative end position of ORF
		public void calcORFOverlap(int start, int end, int[] first_last) {
			//could consider correcting here to keep in-fram
			int first = -1;
			int last = -1;
			for(int i=0 ; i<orfs.length; i++) {
				int[] st_en = orfs[i];
				start_offset[i] = start - st_en[0];
				end_offset[i] = end - st_en[1];
				overlap[i] =  Math.min( end - st_en[0], st_en[1] - start);
				if(overlap[i]>0) {
					if(first<0) first = i;
					last = i;
				}
				//overlap[i] = overl/(st_en[1] - st_en[0]);
			}
			first_last[0] =first;
			first_last[1] = last;
		}
	}
	
	public static class IdentityProfile1 {
		final File outfile, outfile1, outfile2, outfile3, outfile4, outfile5;
		
		
		
		public IdentityProfile1(Sequence refSeq, File bam, int genome_index, int round, boolean calculateCoExpression, double overlapThresh, int startThresh, int endThresh) throws IOException {
			//File readClusterFile = new File(outdir, "readclusters.txt.gz");
			this.round = (double) round;
			this.genome = refSeq;
			this.startThresh = startThresh; this.endThresh = endThresh;
					
			this.calculateCoExpression = calculateCoExpression;
			 outfile = new File(bam.getParentFile(), bam.getName() + genome_index+ ".txt");
			 outfile1 = new File(bam.getParentFile(), bam.getName()  + genome_index+ "coref.txt");
			 outfile2 = new File(bam.getParentFile(), bam.getName()  + genome_index+"clusters.txt");
			 outfile3 = new File(bam.getParentFile(), bam.getName()  + genome_index+ "readToCluster.txt.gz");
			 outfile4 = new File(bam.getParentFile(), bam.getName()  + genome_index+ "exons.txt.gz");
			 outfile5 = new File(bam.getParentFile(), bam.getName()  + genome_index+ "clusters.fa.gz");
			 readClusters = new PrintWriter[nmes.length];
			 for(int index1 =0; index1<nmes.length; index1++) {
				 File outfile3_ = new File(outfile3.getParentFile(),
						outfile3.getName() + "." + round + nmes[index1] + ".gz");
				 this.readClusters[index1] = new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile3_))));
					readClusters[index1].println("readID,clusterId");
			 }
		
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
			readBase = 0;// the number of bases from ref and read
			// following should be made more efficient
			Set<Integer> roundedPos = new HashSet<Integer>();
			for (int i = 0; i < refSeq.length(); i++) {
				roundedPos.add(round(i));
			}
			roundedPositions = roundedPos.toArray(new Integer[0]);
			/*
			 * for(int i=0; i<roundedPositions.length; i++){ roundedPositions[i] =
			 * 1+i*(int)round; }
			 */
			codepth = new SparseRealMatrix[nmes.length];
			all_clusters = new CigarClusters[nmes.length];

			for (int i = 0; i < this.codepth.length; i++) {
				if(this.calculateCoExpression) {
					codepth[i] = new OpenMapRealMatrix(roundedPositions.length, roundedPositions.length);
				}
				all_clusters[i] = new CigarClusters(overlapThresh);
			}
		}

		private int round(int pos) {
			int res = (int) Math.floor((double) pos / round);
			return res;
		}

		final int startThresh, endThresh;
		final  double round ;
		private final boolean calculateCoExpression;

		static String[] nmes = "5_3:5_no3:no5_3:no5_no3".split(":");

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
			String clusterID = this.all_clusters[index].matchCluster(coRefPositions); // this also clears current cluster
			this.readClusters[index].println(id+","+clusterID);
		}

		public void addRefPositions(int position) {
			coRefPositions.add(round(position));
		}

		public SparseRealMatrix[] codepth;
		private CigarCluster coRefPositions = new CigarCluster("reuseable");
		private CigarClusters[] all_clusters;
		private Integer[] roundedPositions;// , corefSum;
		public int[] match, mismatch, refClipped, baseDel, baseIns;
		public int numIns, numDel, readClipped, refBase, readBase;

		private final PrintWriter[]  readClusters;
		private final Sequence genome;
		
		public void print(PrintWriter pw, Sequence seq) {
			pw.println("pos,base,match,mismatch,refClipped,baseDel,baseIns");
			for (int i = 0; i < match.length; i++) {
				pw.println((i + 1) + "," + seq.charAt(i) + "," + match[i] + "," + mismatch[i] + "," + refClipped[i]
						+ "," + baseDel[i] + "," + baseIns[i]);
			}
			pw.flush();
		}

		public void printClusters(File outfile1) throws IOException {
			for (int index = 0; index < this.codepth.length; index++) {
				File outfile1_ = new File(outfile1.getParentFile(),
						outfile1.getName() + nmes[index] + "_" + round + ".gz");
				PrintWriter pw = new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_))));
				StringBuffer sb = new StringBuffer("pos");
				for (int i = 0; i < this.roundedPositions.length; i++) {
					sb.append(",");
					sb.append(roundedPositions[i] * round + 1);
				}
				pw.println(sb.toString());
				for (int i = 0; i < this.all_clusters[index].l.size(); i++) {
					pw.println(this.all_clusters[index].l.get(i).summary(this.roundedPositions));
				}
				pw.close();
			}
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
			
			for(int i=0; i<all_clusters.length; i++) {
				File outfile3_ = new File(outfile4.getParentFile(),
						outfile4.getName() + "." + nmes[i] + ".gz");
				
				PrintWriter exonsP =  new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile3_))));
				
				File outfile5_ = new File(outfile5.getParentFile(),
						outfile5.getName() + "." + nmes[i] + ".gz");
				
				PrintWriter seqFasta =  new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile5_))));
				this.all_clusters[i].getConsensus(annot, this.genome, this.roundedPositions, exonsP, seqFasta);
				exonsP.close();
				seqFasta.close();
			}
			
		}

		public void finalise()  throws IOException{
			for(int i=0; i<readClusters.length; i++) {
			this.readClusters[i].close();
			}
			IdentityProfile1 pr1 = this;
			PrintWriter pw = new PrintWriter(new FileWriter(outfile));
			pr1.print(pw, genome);
			pw.close();
			pr1.printCoRef(outfile1);
			pr1.printClusters(outfile2);
			System.out.println("========================= TOTAL ============================");
			
		}

	}

	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * 
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public static void identity1(Sequence refSeq, Sequence readSeq, SAMRecord sam, IdentityProfile1 profile) {


		int readPos = 0;// start from 0
		int refPos = sam.getAlignmentStart() - 1;// convert to 0-based index
		//String id = sam.getHeader().getId();
		String id = sam.getReadName();
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
