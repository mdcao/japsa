package japsadev.util;

import java.io.File;
import java.io.FileOutputStream;
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


public class HTSUtilities {

public static class CigarClusters{
	static double thresh = 0.9;
	
	List<CigarCluster> l = new ArrayList<CigarCluster>();
	
	public void matchCluster(CigarCluster c1){
		double[] sc = new double[l.size()];
		int min_index = 0;
		double best_score = 0;
		for(int i=0; i<l.size(); i++){
			sc[i] = l.get(i).similarity(c1);
			if(sc[i] > best_score){
				min_index = i;
				best_score = sc[i];
			}
		}
		if(min_index<sc.length && sc[min_index] >= thresh){
			l.get(min_index).merge(c1);
		}else{
		
			CigarCluster newc = new CigarCluster("ID"+l.size());
			newc.merge(c1);
			l.add(newc);
			System.err.println("new cluster "+best_score+" "+min_index);
			System.err.println(l.size());
			
		}
		c1.map.clear();
	}
	
	
}
	
public static class CigarCluster{
	
	final String id;
	
	public CigarCluster(String id){
		this.id = id;
	}
	
	private SortedMap<Integer, Integer> map = new TreeMap<Integer, Integer>();
	public void add(int round) {
		map.put(round, map.containsKey(round) ? map.get(round)+1 : 1);
	}
	public String toString(){
		return map.keySet().toString();
	}
	public Iterator<Integer> keys() {
		return map.keySet().iterator();
	}
	public Iterator<Integer> tailKeys(int st) {
		return map.tailMap(st).keySet().iterator();
	}

	public String summary(Integer[] positions){
		StringBuffer sb = new StringBuffer(id);
		for(int i=0; i<positions.length; i++){
			int i1 = positions[i];
			int v = map.containsKey(i1) ? map.get(i1) : 0;
			sb.append(",");
			sb.append(v);
		}
		return sb.toString();
	}
	
	public double similarity(CigarCluster c1){
		//if(true) return 1.0;
		int intersection =0; 
		int union =map.size();
		for(Iterator<Integer> it = c1.map.keySet().iterator(); it.hasNext(); ){
			if(map.containsKey(it.next())){
				intersection ++;
				
			}else{
				union++;
			}
		}
		return ((double)intersection)/(double) union;
		//Set<Integer> union = Stream.concat(setA.stream(), setB.stream()).collect(Collectors.toSet());
		//Set<Integer> intersect = setA.stream().filter(setB::contains).collect(Collectors.toSet());
		//return ((double)intersect.size())/((double) union.size());
	}
	
	
	
	public void merge(CigarCluster c1){
		Iterator<Entry<Integer, Integer>> it = c1.map.entrySet().iterator();
		while(it.hasNext()){
			Entry<Integer, Integer> entry = it.next();
			Integer key  =entry.getKey();
			int curr =  map.containsKey(key) ? map.get(key): 0;
			map.put(key, curr + entry.getValue());
		}
	}
}
	
public static class IdentityProfile1{
		
		public IdentityProfile1(Sequence refSeq) {
			readClipped = 0;
			numDel = 0;
			numIns = 0;
			match  = new int[refSeq.length()];
			mismatch  = new int[refSeq.length()];
			refClipped = new int[refSeq.length()];
			baseDel = new int[refSeq.length()];
			baseIns = new int[refSeq.length()];
			
			Arrays.fill(match , 0);
			Arrays.fill(mismatch , 0);
			Arrays.fill(baseDel , 0);
			Arrays.fill(baseIns , 0);
			Arrays.fill(refClipped , 0);

			refBase = 0;
			readBase = 0;//the number of bases from ref and read
			//following should be made more efficient
			Set<Integer>roundedPos = new HashSet<Integer>();
			for(int i=0; i<refSeq.length(); i++){
				roundedPos.add(round(i));
			}
			roundedPositions = roundedPos.toArray(new Integer[0]);
			/*for(int i=0; i<roundedPositions.length; i++){
				roundedPositions[i] = 1+i*(int)round;
			}*/
			codepth = new SparseRealMatrix[nmes.length];
			all_clusters = new CigarClusters[nmes.length];
			
			for(int i=0; i<this.codepth.length; i++){
				codepth[i] = new OpenMapRealMatrix(roundedPositions.length, roundedPositions.length);
				all_clusters[i] =  new CigarClusters();
			}
		}
		private int round(int pos){
			int res =  (int) Math.floor((double)pos/round);
			return res;
		}
	
		
		static int refThresh = 80;  //
		static double round = 10.0;
		
		static String[] nmes = "5_3:5_no3:no5_3:no5_no3".split(":");
		
		public void processRefPositions(int startPos, int distToEnd) {
			int index =0;
			if(startPos<refThresh) index = distToEnd < refThresh ? 0 : 1;
			else index = distToEnd < refThresh ? 2 : 3;
			Iterator<Integer> it = coRefPositions.keys();
			while(it.hasNext()){
				Integer pos1 = it.next();
				Iterator<Integer> it2 = coRefPositions.tailKeys(pos1);
				while(it2.hasNext()){
					Integer pos2 = it2.next();
					double value = codepth[index].getEntry(pos1, pos2);
					this.codepth[index].setEntry(pos1, pos2, value+1);
				}
			}
			this.all_clusters[index].matchCluster(coRefPositions);  // this also clears current cluster
		}

	
		public void addRefPositions(int position) {
			coRefPositions.add(round(position));
		}

		public SparseRealMatrix[] codepth;
		private CigarCluster coRefPositions = new CigarCluster("reuseable");
		private CigarClusters[] all_clusters;
		private Integer[] roundedPositions;//, corefSum;
		public int[] match, mismatch, refClipped, baseDel, baseIns;
		public int  numIns, numDel,  readClipped, refBase, readBase;
		
		public void print(PrintWriter pw, Sequence seq) {
			pw.println("pos,base,match,mismatch,refClipped,baseDel,baseIns");
			for(int i=0;i<match.length; i++){
				pw.println((i+1)+","+seq.charAt(i)+","+match[i]+","+mismatch[i]+","+refClipped[i]+","+baseDel[i]+","+baseIns[i]);
			}
			pw.flush();
		}
		
		public void printClusters(File outfile1) throws IOException{
			for(int index=0; index < this.codepth.length; index++){
				File outfile1_ = new File(outfile1.getParentFile(), outfile1.getName()+nmes[index]+"_"+round+".gz");
				PrintWriter pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_))));
				StringBuffer sb = new StringBuffer("pos");
				for(int i=0 ; i<this.roundedPositions.length; i++){
					sb.append(",");
					sb.append(roundedPositions[i]*round+1);
				}
				pw.println(sb.toString());
				for(int i=0; i<this.all_clusters[index].l.size(); i++){
					pw.println(this.all_clusters[index].l.get(i).summary(this.roundedPositions));
				}
				pw.close();
			}
		}
		
		public void printCoRef(File outfile1) throws  IOException {
			for(int index=0; index < this.codepth.length; index++){
				File outfile1_ = new File(outfile1.getParentFile(), outfile1.getName()+"."+round+nmes[index]+".gz");
			PrintWriter pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_))));
			int len = this.codepth[index].getRowDimension();
			
			for(int i=0;i<len; i++){
				int row = i; //nonZeroRows.get(i);
				pw.print(roundedPositions[row]*round+1);
				for(int j=0; j<len; j++){
					int col = j;//nonZeroRows.get(j);
					int val = col>=row ? (int) this.codepth[index].getEntry(row, col) : 0;
						//(int) this.codepth[index].getEntry(col, row);
					pw.print(",");pw.print(val); 
				}
				pw.println();
			}
			pw.close();
			}
			
		}

	}
	
/**
 * Get the identity between a read sequence from a sam and a reference sequence
 * @param refSeq
 * @param sam
 * @return
 */
public static void identity1(Sequence refSeq, Sequence readSeq,  SAMRecord sam, IdentityProfile1 profile){
//	IdentityProfile1 profile = new IdentityProfile1(refSeq);
	//profile.refClipped += sam.getAlignmentStart() + refSeq.length() - sam.getAlignmentEnd();
	
	int readPos = 0;//start from 0					
	int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index				

	

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
			//System.err.println(length);
			refPos += length; 
			for (int i = 0; i < length && refPos + i < refSeq.length(); i++){
				profile.refClipped[refPos+i] += 1;
			}
			//profile.refClipped += length;
			break;  // reference skip

		case D ://deletion      	
			refPos += length;
			profile.refBase += length;
			for (int i = 0; i < length && refPos + i < refSeq.length(); i++){
				profile.baseDel[refPos+i] +=1;
			}
			profile.numDel ++;
			break; 	

		case I :	                	
			readPos += length;
			profile.readBase += length;

			profile.baseIns[refPos] += length;
			profile.numIns ++;
			break;
		case M :
			for (int i = 0; i < length && refPos + i < refSeq.length(); i++){
				profile.addRefPositions(refPos+i);
				if (refSeq.getBase(refPos + i) == readSeq.getBase(readPos + i))
					profile.match[refPos+i]++;
				else
					profile.mismatch[refPos+i]++;
			}
			profile.readBase += length;
			profile.refBase += length;

			readPos += length;
			refPos  += length;
			break;

		case EQ :
			readPos += length;
			refPos  += length;
			profile.addRefPositions(refPos);
			profile.readBase += length;
			profile.refBase += length;
			profile.match[refPos] += length;
			break;

		case X :
			readPos += length;
			refPos  += length;

			profile.readBase += length;
			profile.refBase += length;
			profile.addRefPositions(refPos);
			profile.mismatch[refPos] += length;
			break;
		default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
		}//case
	}//for			
	profile.processRefPositions(sam.getAlignmentStart(), refSeq.length()-sam.getAlignmentEnd());
	//\return profile;

}


}
