package japsadev.util;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Sequence;

public class HTSUtilities {

public static class IdentityProfile1{
		
		public IdentityProfile1(Sequence refSeq) {
			readClipped = 0;
		//	refClipped = sam.getAlignmentStart() + refSeq.length() - sam.getAlignmentEnd();
			
			
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
			roundedPositions = new int[1+round(refSeq.length())];
			for(int i=0; i<roundedPositions.length; i++){
				roundedPositions[i] = 1+i*10;
			}
			codepth = new SparseRealMatrix[nmes.length];
			for(int i=0; i<this.codepth.length; i++){
				codepth[i] = new OpenMapRealMatrix(roundedPositions.length, roundedPositions.length);
			}
			//corefSum = new int[roundedPositions.length];
			//Arrays.fill(corefSum, 0);
		}
		private int round(int pos){
			//return pos;
			int res =  (int) Math.floor((double)pos/10.0);
			//System.err.println(pos+"->"+res);
			return res;
		}
		//static int round = 10;
		
	/*	private void checkDiaganol(int i){
			System.err.println("checking diag" +i);
			int dim = this.codepth.getColumnDimension();
				double sc = codepth.getEntry(i, i);
				for(int j=0; j<dim; j++ ){
					if(codepth.getEntry(i, j)>sc) throw new RuntimeException(sc+" "+i+" "+j+" "+codepth.getEntry(i, j));//+" "+corefSum[i]);
				}
		}*/
		
		static int refThresh = 80;  //
		
		static String[] nmes = "5_3:5_no3:no5_3:no5_no3".split(":");
		
		public void processRefPositions(int startPos, int distToEnd) {
			int index =0;
			if(startPos<refThresh) index = distToEnd < refThresh ? 0 : 1;
			else index = distToEnd < refThresh ? 2 : 3;
			//if(index==1) System.err.println("found non leader read");
			Iterator<Integer> it = coRefPositions.iterator();
			while(it.hasNext()){
				Integer pos1 = it.next();
				Iterator<Integer> it2 = coRefPositions.tailSet(pos1).iterator();
				//corefSum[pos1]+=1;
				while(it2.hasNext()){
					Integer pos2 = it2.next();
					double value = codepth[index].getEntry(pos1, pos2);
					this.codepth[index].setEntry(pos1, pos2, value+1);
				}
			}
			//this.checkDiaganol();
			coRefPositions.clear();
		}

	
		public void addRefPositions(int position) {
			coRefPositions.add(round(position));
			
		}

		public SparseRealMatrix[] codepth;
		private SortedSet<Integer> coRefPositions = new TreeSet<Integer>();
		private int[] roundedPositions;//, corefSum;
		public int[] match, mismatch, refClipped, baseDel, baseIns;
		public int  numIns, numDel,  readClipped, refBase, readBase;
		
		public void print(PrintWriter pw, Sequence seq) {
			pw.println("pos,base,match,mismatch,refClipped,baseDel,baseIns");
			for(int i=0;i<match.length; i++){
				pw.println((i+1)+","+seq.charAt(i)+","+match[i]+","+mismatch[i]+","+refClipped[i]+","+baseDel[i]+","+baseIns[i]);
			}
			pw.flush();
		}
		
		
		
		public void printCoRef(File outfile1) throws  IOException {
			for(int index=0; index < this.codepth.length; index++){
				File outfile1_ = new File(outfile1.getParentFile(), outfile1.getName()+"."+nmes[index]+".gz");
			PrintWriter pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_))));
			
			
			/*List<Integer> nonZeroRows = new ArrayList<Integer>();
			for(int i=0; i<corefSum.length; i++){
				if(corefSum[i]>=thresh){
					nonZeroRows.add(i);
				}
			}*/
			int len = this.codepth[index].getRowDimension();
			
			for(int i=0;i<len; i++){
				int row = i; //nonZeroRows.get(i);
				pw.print(roundedPositions[row]);
			//	this.checkDiaganol(row);
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
