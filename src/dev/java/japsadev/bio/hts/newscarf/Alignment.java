package japsadev.bio.hts.newscarf;

import java.util.ArrayList;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Sequence;

public class Alignment implements Comparable<Alignment> {
	public final static int OVERHANG_THRES=1000; 
	
	int score;

	public String readID;
	BidirectedNode node;

	public int refStart, refEnd;  //1-based position on ref of the start and end of the alignment
	
	//Position on read of the start and end of the alignment (using the direction of read) 
	//readStart map to refStart, readEnd map to refEnd. 
	//readStart < readEnd if strand = true, else readStart > readEnd
	public int readStart = 0, readEnd = 0;	
	
	//read length
	public int readLength = 0;		

	public boolean strand = true;//positive
	public boolean useful = false;
	//SAMRecord mySam;
	
	ArrayList<CigarElement> alignmentCigars = new ArrayList<CigarElement>();
	

	//public int readLeft, readRight, readAlign, refLeft, refRight, refAlign;
	//left and right are in the direction of the reference sequence
	
	public Alignment(String readID, int refStart, int refEnd, int readLength, 
			int readStart, int readEnd, boolean strand, boolean useful, BidirectedNode node, int score){
		this.readID = readID;
		this.node = node;
		this.refStart = refStart;
		this.refEnd = refEnd;
		
		this.readLength = readLength;
		this.readStart = readStart;//1-index
		this.readEnd = readEnd;//1-index
		this.strand = strand;
		this.useful = useful;			
		this.score = score;
	}
	public Alignment(SAMRecord sam, BidirectedNode node) {
//		readID = Integer.parseInt(sam.getReadName().split("_")[0]);
		readID = sam.getReadName();

		this.node = node;

		refStart = sam.getAlignmentStart();
		refEnd = sam.getAlignmentEnd();
		
		Cigar cigar = sam.getCigar();			
		boolean enterAlignment = false;						
		//////////////////////////////////////////////////////////////////////////////////

		for (final CigarElement e : cigar.getCigarElements()) {
			alignmentCigars.add(e);
			final int  length = e.getLength();
			switch (e.getOperator()) {
			case H :
			case S :					
			case P : //pad is a kind of clipped
				if (enterAlignment)
					readEnd = readLength;
				readLength += length;
				break; // soft clip read bases
			case I :	                	
			case M :					
			case EQ :
			case X :
				if (!enterAlignment){
					readStart = readLength + 1;
					enterAlignment = true;
				}
				readLength += length;
				break;
			case D :
			case N :
				if (!enterAlignment){
					readStart = readLength + 1;
					enterAlignment = true;
				}
				break;				
			default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}//case
		}//for
		if (readEnd == 0)
			readEnd = readLength;
		//these temporary variable to determine usefulness
		int readLeft = readStart -1;
		int readRight = readLength - readEnd;

		int refLeft = refStart - 1;
		int refRight = ((Sequence) node.getAttribute("seq")).length() - refEnd;
		
		score = refEnd + 1 - refStart;
		if (sam.getReadNegativeStrandFlag()){			
			strand = false;
			//need to convert the alignment position on read the correct direction 
			readStart = 1 + readLength - readStart;
			readEnd = 1 + readLength - readEnd;
		}

		if (
				(readLeft < OVERHANG_THRES || refLeft < OVERHANG_THRES) &&
				(readRight  < OVERHANG_THRES || refRight < OVERHANG_THRES) &&
				score > BidirectedGraph.getKmerSize()
			)
			useful = true;

	}
	
	
	public int readAlignmentStart(){
		return Math.min(readStart,readEnd);
	
	}
	
	public int readAlignmentEnd(){
		return Math.max(readStart,readEnd);
	}

	public String toString() {
		return node.getId()    
				+ " " + refStart 
				+ " " + refEnd
				+ " " + ((Sequence) node.getAttribute("seq")).length()
				+ " " + readStart 
				+ " " + readEnd
				+ " " + readLength				
				+ " " + strand;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(Alignment o) {			
		return readAlignmentStart() - o.readAlignmentStart();
	}
}
