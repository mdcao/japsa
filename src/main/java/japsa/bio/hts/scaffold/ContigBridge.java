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
 * 19/12/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.hts.scaffold;



import htsjdk.samtools.CigarElement;
import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import japsa.bio.np.ErrorCorrection;

/**
 * Create a bridge that connects two contigs. The bridged can be ranked based
 * on the confidence so that more confident bridge is used before.
 * Note that two contigs can have more than one bridges from circular 
 * sequence or false positives. 
 * @author minhduc
 *
 */

public class ContigBridge implements Comparable<ContigBridge>{

	final Contig firstContig, secondContig;
	final String hashKey;
	final int orderIndex;

	private double score = 0;//more is better
	private ScaffoldVector transVector = null;
	private ArrayList<Connection> connections;//a list of connections that make up this

	public ContigBridge(Contig c1, Contig c2, int ind){
		firstContig = c1;
		secondContig = c2;
		orderIndex = ind;		
		hashKey = makeHash(c1.index,c2.index, orderIndex);

		connections = new  ArrayList<Connection>();
	}

	public static String makeHash(int aIndex, int bIndex, int order){		
		return aIndex+"#"+bIndex + "#" + order;
	}

	public boolean consistentWith(ScaffoldVector aVector){
		int tolerance = firstContig.getIndex()==secondContig.getIndex()?100:250;
		return (aVector.direction == transVector.direction)
				&& (((aVector.magnitude * 1.0 / transVector.magnitude > 0.75)
				&& (aVector.magnitude * 1.0 / transVector.magnitude < 1.25)) 
				|| (Math.abs(aVector.magnitude-transVector.magnitude) < tolerance)
				)
				;
	}

	public double addConnection(ReadFilling readSequence, 
			AlignmentRecord firstAlignment, 
			AlignmentRecord secondAlignment, 
			ScaffoldVector trans, 
			double sc){
	
		//NB: firstAlignment for firstContig, secondAlignment for secondContig

		if (transVector == null){
			transVector = trans;
			score = sc;			
			connections.add(new Connection(readSequence, firstAlignment,secondAlignment,trans));			
		}else{		
			Connection newConnect = new Connection(readSequence, firstAlignment,secondAlignment,trans);
			connections.add(newConnect);			
			transVector.magnitude = (transVector.magnitude * connections.size() + trans.magnitude) / (connections.size() + 1);

			//the metric for bridge score is important!
			score += sc;
			//score = score>sc?score:sc;
		}
		return score;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	public void setScore(double s) {
		score = s;
	}
	
	//TODO: magnitude usually don't help for bridge with repeat.
	// E.g. <--===---------------> prev not next for the both
	public void setContigScores(){
		int 	firstPointer = 0,
				secondPointer = 0;
		
		if(transVector.magnitude < 0){
			firstPointer=-1;
			if(transVector.direction < 0)
				secondPointer=-1;
			else
				secondPointer=1;
		}
		// special case: magnitude < firstContig.length() && transVector.direction < 0;
		else if(transVector.magnitude < firstContig.length() && transVector.direction < 0){
			firstPointer = secondPointer = -1;
		}
		else{
			
			firstPointer=1;
			if(transVector.direction > 0)
				secondPointer=-1;
			else
				secondPointer=1;
		}
		//reset based on the pointers
		if(firstPointer > 0){
			firstContig.nextScore = score;
			if(ScaffoldGraph.verbose)
				System.out.printf("...set nextScore of %s to %.2f\n", firstContig.getName(), score);
		}else{
			firstContig.prevScore = score;
			if(ScaffoldGraph.verbose)
				System.out.printf("...set prevScore of %s to %.2f\n", firstContig.getName(), score);	
		}
		
		if(secondPointer > 0){
			secondContig.nextScore = score;
			if(ScaffoldGraph.verbose)
				System.out.printf("...set nextScore of %s to %.2f\n", secondContig.getName(), score);
		}else{
			secondContig.prevScore = score;
			if(ScaffoldGraph.verbose)
				System.out.printf("...set prevScore of %s to %.2f\n", secondContig.getName(), score);	
		}
		
	}
	// when contig bridge is removed, reset the scores
	public void resetContigScores(){
		int 	firstPointer = 0,
				secondPointer = 0;
		if(ScaffoldGraph.verbose)
			System.out.print("Trans vector " + transVector + " :" );
		if(transVector.magnitude < 0){
			firstPointer=-1;
			if(transVector.direction < 0)
				secondPointer=-1;
			else
				secondPointer=1;
		}
		// special case: magnitude < firstContig.length() && transVector.direction < 0;
		else if(transVector.magnitude < firstContig.length() && transVector.direction < 0){
			firstPointer = secondPointer = -1;
		}
		else{
			firstPointer=1;
			if(transVector.direction > 0)
				secondPointer=-1;
			else
				secondPointer=1;
		}
		//reset based on the pointers
		if(firstPointer > 0){
			firstContig.nextScore = .0;
			if(ScaffoldGraph.verbose)
				System.out.printf("...reset nextScore of %s to 0, ", firstContig.getName());
		}else{
			firstContig.prevScore = .0;
			if(ScaffoldGraph.verbose)
				System.out.printf("...reset prevScore of %s to 0, ", firstContig.getName());	
		}
		
		if(secondPointer > 0){
			secondContig.nextScore = .0;
			if(ScaffoldGraph.verbose)
				System.out.printf("reset nextScore of %s to 0\n", secondContig.getName());
		}else{
			secondContig.prevScore = .0;
			if(ScaffoldGraph.verbose)
				System.out.printf("reset prevScore of %s to 0\n", secondContig.getName());	
		}
			
	}
	/**
	 * @return the transVector
	 */
	public ScaffoldVector getTransVector() {
		return transVector;
	}
	/**
	 * @return the connections
	 */
	public ArrayList<Connection> getConnections() {
		return connections;
	}	

	public Connection fewestGapConnection() throws IOException{
		Collections.sort(connections);
		
		//Find the best connections (has fewest gaps)
		Connection gapsBestConnection = null;
		int gapsBest = Integer.MAX_VALUE;
		for (Connection connection:connections){
			int gapsBt = connection.gapsBetween();
			if (gapsBt < gapsBest){
				gapsBest = gapsBt;
				gapsBestConnection = connection;
			}
		}	
		
		return gapsBestConnection;

	}
	/**
	 * Try to connect contigs with consensus sequence from involved reads
	 * TODO optimized the code
	 * @return 
	 * @throws IOException
	 */
	public Connection consensusConnection(SequenceOutputStream consensusOut) throws IOException{
		int offset = 100; //1-based
		Collections.sort(connections);
		ArrayList<Sequence> readList = new ArrayList<Sequence>(connections.size());
		// locate the offset points on two contigs
		int cutOnFirstContig, cutOnSecondContig;
		int tS = 1, tE = firstContig.length(),
				fS, fE;			
		if (transVector.direction > 0){
			fS = transVector.magnitude;
			fE = transVector.magnitude + secondContig.length();
		}else{
			fE = transVector.magnitude;
			fS = transVector.magnitude - secondContig.length();
		}		
		// tS---|->tE fS<-|--->fE
		if (fS-tE > tS-fE){
			cutOnFirstContig = firstContig.length()>offset?firstContig.length()-offset:firstContig.length();			
			cutOnSecondContig = secondContig.length()>offset?offset:secondContig.length();
			
		}
		// fS<---|->fE tS-|--->tE
		else{
			cutOnFirstContig = firstContig.length()>offset?offset:firstContig.length();
			cutOnSecondContig = secondContig.length()>offset?secondContig.length()-offset:secondContig.length();
			
		}
		// should we check other case (overlapped, contained..)??
		cutOnSecondContig = transVector.direction>0?cutOnSecondContig:secondContig.length()-cutOnSecondContig;
		
		Connection gapsBestConnection = null;
		int gapsBest = Integer.MAX_VALUE;
		int rplStart=0, rplEnd=0;
		for (Connection connection:connections){
			int 	firstCutOnRead=mapToRead(cutOnFirstContig, connection.firstAlignment),
					secondCutOnRead=mapToRead(cutOnSecondContig, connection.secondAlignment);
			Sequence tmp = null;
			try{
//				if(firstCutOnRead < secondCutOnRead)
//					tmp = connection.read.readSequence.subSequence(firstCutOnRead, secondCutOnRead);
//				else
//					tmp = connection.read.readSequence.subSequence(secondCutOnRead, firstCutOnRead);
				ReadFilling tmpRead = connection.read;
				if (firstCutOnRead > secondCutOnRead){
					connection.read = connection.read.reverse();
					connection.firstAlignment=connection.firstAlignment.reverseRead();
					connection.secondAlignment=connection.secondAlignment.reverseRead();
					firstCutOnRead = tmpRead.readSequence.length()-firstCutOnRead+1;
					secondCutOnRead = tmpRead.readSequence.length()-secondCutOnRead+1;
				}
				
				tmp = connection.read.readSequence.subSequence(firstCutOnRead-1, secondCutOnRead-1);
				tmp.setName(tmpRead.readSequence.getName());
				tmp.setDesc(tmpRead.readSequence.getDesc());
				readList.add(tmp);
			}
			catch(Exception e){
				e.printStackTrace();
				System.err.println("Failed attempt to extract (" + firstCutOnRead + ", " + secondCutOnRead 
									+ ") from sequence with length " + connection.read.readSequence.length());
			}
			int gapsBt = connection.gapsBetween();
			if (gapsBt < gapsBest){
				gapsBest = gapsBt;
				gapsBestConnection = connection;
				rplStart = firstCutOnRead;
				rplEnd = secondCutOnRead;
			}
		}

		Sequence consensus = null;
		Connection consensusConnection = gapsBestConnection;
		try {
			consensus = ErrorCorrection.consensusSequence(readList, hashKey, "poa");
			consensus.setName(hashKey);
			consensus.setDesc("Consensus sequence");
			consensus.writeFasta(consensusOut);
			Sequence gapsBestSequence = gapsBestConnection.read.readSequence;
			int len = gapsBestSequence.length()-(rplEnd-rplStart+1)+consensus.length();
			Sequence rpl = new Sequence(Alphabet.DNA16(), len);
			for (int idx=0;idx < len; idx++){
				if (idx < rplStart)
					rpl.setBase(idx, gapsBestSequence.getBase(idx));
				else if (idx >= rplStart+consensus.length())
					rpl.setBase(idx, gapsBestSequence.getBase(idx+rplEnd+1-rplStart-consensus.length()));
				else
					rpl.setBase(idx, consensus.getBase(idx-rplStart));
			}
			System.out.println("---->Changed length (loss): " + (len-gapsBestSequence.length()));
			
			AlignmentRecord newFirst=gapsBestConnection.firstAlignment,
							newSecond=gapsBestConnection.secondAlignment;
			newSecond.readStart+=len-gapsBestSequence.length();
			newSecond.readEnd+=len-gapsBestSequence.length();
			ArrayList<AlignmentRecord> ends = new ArrayList<AlignmentRecord>();
			ends.add(newFirst);
			ends.add(newSecond);
			ReadFilling simple = new ReadFilling(rpl, ends);
			consensusConnection = new Connection(simple,gapsBestConnection.firstAlignment,gapsBestConnection.secondAlignment,gapsBestConnection.trans);
			
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("Can not generate consensus sequence!");
		}
		return consensusConnection;		
		
	}

	/**
	 * Return the position on the reference that corresponds to a given position
	 * on read.
	 *  
	 * @param posInRead
	 * @param record
	 * @return
	 */
	static int positionOnRef(int posOnRead, AlignmentRecord record){
		if (posOnRead < record.readAlignmentStart() || posOnRead > record.readAlignmentEnd())
			return 0;

		if (!record.strand)
			posOnRead = record.readLength - posOnRead + 1; // use direction of ref (forward)


		int pos = record.strand?record.readStart:(record.readLength + 1 - record.readStart);
		int posOnRef = record.refStart;

		//assert pos <= posOnRead
		for (final CigarElement e : record.alignmentCigars) {
			final int  length = e.getLength();
			switch (e.getOperator()) {
			case H :
			case S :					
			case P :
				break; // ignore pads and clips
			case I :				
				//insert
				if (pos + length < posOnRead){
					pos += length;				
				}else{
					return posOnRef;
				}
				break;
			case M ://match or mismatch				
			case EQ://match
			case X ://mismatch
				if (pos + length < posOnRead){
					pos += length;
					posOnRef += length;
				}else{
					return posOnRef + posOnRead - pos;
				}
				break;
			case D :
				posOnRef += length;
				break;
			case N :	
				posOnRef += length;
				break;								
			default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}//casse
		}//for		

		return 0;
	}
	/**
	 * Return the position on the read that corresponds to a given position
	 * on reference.
	 *  
	 * @param posInRef
	 * @param record
	 * @return
	 */
	static int mapToRead(int posOnRef, AlignmentRecord record){
		// read htsjdk.samtools.* API
		int location = -1;
		
		if ((posOnRef - record.refStart)*(posOnRef - record.refEnd) >= 0){
			if (Math.abs(posOnRef-record.refStart) > Math.abs(posOnRef-record.refEnd))
				location = record.strand?record.readEnd+posOnRef-record.refEnd:record.readEnd-posOnRef+record.refEnd;			
			else
				location = record.strand?record.readStart+posOnRef-record.refStart:record.readStart-posOnRef+record.refStart;
		}
		else{
			// current coordinate on read, followed the reference contig's direction
			int posOnRead = record.strand?record.readStart:record.readLength-record.readStart+1;
			 // current position on ref 
			int pos = record.refStart;
			
			for (final CigarElement e : record.alignmentCigars) {
				final int  length = e.getLength();
				switch (e.getOperator()) {
				case H :
				case S :					
				case P :
					break; // ignore pads and clips
				case I :			
					posOnRead += length;
					break;	
				case M ://match or mismatch				
				case EQ://match
				case X ://mismatch
					if (pos + length < posOnRef){
						pos += length;
						posOnRead += length;
					}else{
						location = posOnRef + posOnRead - pos;
					}
					break;
				case D :
				case N :	
					//delete
					if (pos + length < posOnRef){
						pos += length;				
					}else{
						location = posOnRead;
					}
					break;	
				default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
				}//casse
			}//for		
			//convert back to coordinate based on read direction
			location = record.strand?location:record.readLength-location+1;
		}
		
		System.out.println( "Contig (ref): " + record.contig.getName() + " Read: " + record.readID + " Strand: " + record.strand);   
		System.out.println( "\tOn contig: " + record.refStart	+ " -> " + record.refEnd
				+ " Len: " + record.contig.length() + " Cut point: " + posOnRef);
		System.out.println( "\tOn read: " + record.readStart	+ " -> " + record.readEnd
				+ " Len: " + record.readLength + " Alleged cut point: " + location);
		
		location=location>0?location:0;
		location=location<record.readLength?location:record.readLength-1;
		
		System.out.println( "\tOn read: " + record.readStart	+ " -> " + record.readEnd
				+ " Len: " + record.readLength + " Final cut point: " + location);
		return location;
	}

	public boolean isContaining(Contig ctg){
		if(firstContig.getIndex() == ctg.getIndex() || secondContig.getIndex() == ctg.getIndex())
			return true;
		else 
			return false;
	}
	public void display(){
		System.out.printf("##################START########################\n"
				+ "Contig %3d (%d) -> Contig %3d (%d) Vector (%s) score = %f distance = %d\n",
				this.firstContig.index,
				this.firstContig.length(),
				this.secondContig.index,
				this.secondContig.length(),
				transVector.toString(),
				this.score,
				transVector.distance(firstContig, secondContig)				
				);

		Collections.sort(connections);
		for (Connection connect:connections)
			connect.display();
		System.out.println("##################END########################");
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(ContigBridge o) {
		return (int) (o.score - score);

	}	
	public ReadFilling consensusRead() throws IOException{
		int offset = 300; //1-based
		Collections.sort(connections);
		ArrayList<Sequence> readList = new ArrayList<Sequence>(connections.size());
		// locate the offset points on two contigs
		int cutOnFirstContig, cutOnSecondContig;
		int tS = 1, tE = firstContig.length(),
				fS, fE;			
		if (transVector.direction > 0){
			fS = transVector.magnitude;
			fE = transVector.magnitude + secondContig.length();
		}else{
			fE = transVector.magnitude;
			fS = transVector.magnitude - secondContig.length();
		}		
		// tS---|->tE fS<-|--->fE
		if (fS-tE > tS-fE){
			cutOnFirstContig = firstContig.length()>offset?firstContig.length()-offset:firstContig.length();			
			cutOnSecondContig = secondContig.length()>offset?offset:secondContig.length();
			
		}
		// fS<---|->fE tS-|--->tE
		else{
			cutOnFirstContig = firstContig.length()>offset?offset:firstContig.length();
			cutOnSecondContig = secondContig.length()>offset?secondContig.length()-offset:secondContig.length();
			
		}
		// cuz first contig direction was used as base -> adjust coordinate on the second
		cutOnSecondContig = transVector.direction>0?cutOnSecondContig:secondContig.length()-cutOnSecondContig;
		
		for (Connection connection:connections){
			int 	firstCutOnRead=mapToRead(cutOnFirstContig, connection.firstAlignment),
					secondCutOnRead=mapToRead(cutOnSecondContig, connection.secondAlignment);
			Sequence tmp = null;
			try{
				ReadFilling tmpRead = connection.read;
				if (firstCutOnRead > secondCutOnRead){
					tmpRead = connection.read.reverse();
					firstCutOnRead = tmpRead.readSequence.length() -  firstCutOnRead;
					secondCutOnRead = tmpRead.readSequence.length() - secondCutOnRead;
				}
				
				tmp = tmpRead.readSequence.subSequence(firstCutOnRead-1, secondCutOnRead-1);
				tmp.setName(tmpRead.readSequence.getName());
				tmp.setDesc(tmpRead.readSequence.getDesc());
				readList.add(tmp);
			}
			catch(Exception e){
				e.printStackTrace();
				System.err.println("Failed attempt to extract (" + firstCutOnRead + ", " + secondCutOnRead 
									+ ") from sequence with length " + connection.read.readSequence.length());
			}

		}

		Sequence consensus = readList.get(0);
		ReadFilling consensusRead = null;
		try {
			consensus = ErrorCorrection.consensusSequence(readList, hashKey, "poa");
			consensusRead = new ReadFilling(consensus, new ArrayList<AlignmentRecord>());
		} catch (InterruptedException e) {
			e.printStackTrace();
			System.err.println("Can not generate consensus sequence!");
		}
		
		return consensusRead;		
		
	}
	/* 
	 * Fill the scaffold considering all connections (get the consensus)
	 * 
	 */
	public Sequence fillConsensus(AlignmentRecord ttAlign, AlignmentRecord ffAlign){
		int tS = 1, tE = firstContig.length(),
				fS, fE, tC, fC;
		AlignmentRecord tAlign = new AlignmentRecord(), 	
						fAlign = new AlignmentRecord();
		if (transVector.direction > 0){
			fS = transVector.magnitude;
			fE = transVector.magnitude + secondContig.length();
		}else{
			fE = transVector.magnitude;
			fS = transVector.magnitude - secondContig.length();
		}		
		// tS---|->tE fS<-|--->fE
		if (fS-tE > tS-fE){
			int tEnd = firstContig.length()-1, fEnd = transVector.direction>0?0:secondContig.length()-1; //furthest pair
			tC=tEnd;
			fC=fEnd;
			
			for (Connection connection:connections){
				if(Math.min(Math.abs(connection.firstAlignment.refStart-tEnd),
							Math.abs(connection.firstAlignment.refEnd-tEnd))
					> Math.abs(tC-tEnd)){
					tC=	Math.abs(connection.firstAlignment.refStart-tEnd) <
						Math.abs(connection.firstAlignment.refEnd-tEnd)?
								connection.firstAlignment.refStart
								:connection.firstAlignment.refEnd;
					tAlign=connection.firstAlignment;
				
				}
				if(Math.min(Math.abs(connection.secondAlignment.refStart-fEnd),
						Math.abs(connection.secondAlignment.refEnd-fEnd))
					> Math.abs(fC-fEnd)){
					fC=	Math.abs(connection.secondAlignment.refStart-fEnd) <
						Math.abs(connection.secondAlignment.refEnd-fEnd)?
						connection.secondAlignment.refStart
						:connection.secondAlignment.refEnd;
					fAlign=connection.secondAlignment;	
				}
			}
			
		}
		// fS<---|->fE tS-|--->tE
		else{
			int tEnd = 0, fEnd = transVector.direction>0?secondContig.length()-1:0; //furthest pair
			tC=tEnd; 
			fC=fEnd;
			
			for (Connection connection:connections){
				if(Math.min(Math.abs(connection.firstAlignment.refStart-tEnd),
							Math.abs(connection.firstAlignment.refEnd-tEnd))
					> Math.abs(tC-tEnd)){
					tC=	Math.abs(connection.firstAlignment.refStart-tEnd) <
						Math.abs(connection.firstAlignment.refEnd-tEnd)?
								connection.firstAlignment.refStart
								:connection.firstAlignment.refEnd;
					tAlign=connection.firstAlignment;

				}
				if(Math.min(Math.abs(connection.secondAlignment.refStart-fEnd),
						Math.abs(connection.secondAlignment.refEnd-fEnd))
					> Math.abs(fC-fEnd)){
					fC=	Math.abs(connection.secondAlignment.refStart-fEnd) <
						Math.abs(connection.secondAlignment.refEnd-fEnd)?
						connection.secondAlignment.refStart
						:connection.secondAlignment.refEnd;
					fAlign=connection.secondAlignment;
				}
			}
		}
		ttAlign.copy(tAlign);
		ffAlign.copy(fAlign);

		//----------------------------------------------------------------------------------
		ArrayList<Sequence> seqList = new ArrayList<Sequence>();
		
		Contig fromContig = firstContig,
				toContig = secondContig;

		//loop over all connections
		for(Connection connection:connections){
	
			AlignmentRecord fromAlignment = 
					(fromContig == firstContig)?connection.firstAlignment:connection.secondAlignment;
	
			AlignmentRecord toAlignment = 
					(fromContig == firstContig)?connection.secondAlignment:connection.firstAlignment;
	
			ReadFilling readFilling = connection.read;
			if ( fromContig.getRelDir()>0 != fromAlignment.strand){
				//swap
				readFilling = connection.read.reverse();
				readFilling.sortAlignment();
				for (AlignmentRecord record:readFilling.alignments){
					if (record.contig == fromContig)
						fromAlignment = record;
	
					if (record.contig == toContig)
						toAlignment = record;	
				}
			}
			//now readFilling is good to go
			int posReadEnd   = fromAlignment.readAlignmentEnd();
			int posReadFinal = toAlignment.readAlignmentStart();// I need as far as posReadFinal
			
			Sequence seqRead = readFilling.readSequence.subSequence(posReadEnd, posReadFinal);
			seqRead.setName(new String("R_" + connection.readID + "_" + score));
			
			SequenceBuilder seqContig = new SequenceBuilder(Alphabet.DNA16(),1024*1024,"C_" + connection.readID + "_" + score);
			//int curPos=0; //current 1-based position pointer of seqContig
			for (AlignmentRecord record:readFilling.alignments){
				Contig contig = record.contig;
				if (contig == fromContig)
					continue;

				contig.addRange(record.refStart,record.refEnd,record.score);
				if (posReadEnd >= posReadFinal -1)
					//continue;//I can break here, but want to get portionUsed of other contigs
					break;
	
				if (record.readAlignmentEnd() < posReadEnd)
					continue;				
				
				//assert:  posReadEnd < readEnd				
				if (record.readAlignmentStart() > posReadEnd){
					//Really need to fill in using read information
					int newPosReadEnd = Math.min(posReadFinal - 1, record.readAlignmentStart() -1);
					if (newPosReadEnd > posReadEnd){
						seqContig.append(readFilling.readSequence.subSequence(posReadEnd, newPosReadEnd));
//						if (connection.readID%3==2)
//							seqContig.append(readFilling.readSequence.subSequence(posReadEnd, newPosReadEnd));
//						else{
//							char[] n = new char[newPosReadEnd-posReadEnd+1];
//							java.util.Arrays.fill(n,'-');
//							seqContig.append(new Sequence(Alphabet.DNA16(),n,"Filling"));
//						}
						posReadEnd = newPosReadEnd;
				
					}
					if (posReadEnd + 1 >= posReadFinal)
						//continue;//Done
						break;
					//Now get information on the contig from start
					if (contig == toContig)
						//continue;//could break
						break;
					if (record.strand){
						int refLeft = record.refStart;
						int refRight = record.refEnd;
	
						if (posReadFinal <= record.readAlignmentEnd()){
							refRight = positionOnRef(posReadFinal, record) -1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
														
						seqContig.append(contig.contigSequence.subSequence(refLeft - 1, refRight));
					}else{//neg strain
						int refRight = record.refStart;
						int refLeft = record.refEnd;
	
						if (posReadFinal <= record.readAlignmentEnd()){
							refLeft = positionOnRef(posReadFinal, record) + 1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
							
						seqContig.append(Alphabet.DNA.complement(contig.contigSequence.subSequence(refRight - 1, refLeft)));
					}
				}//if record.readAlignmentStart() > posReadEnd
				else{//Now get information on the contig from start
					if (contig == toContig)
						//continue;//could break
						break;
					if (record.strand){
						int refLeft = positionOnRef(posReadEnd, record) + 1;						
						int refRight = record.refEnd;
	
						if (posReadFinal <= record.readAlignmentEnd()){
							refRight = positionOnRef(posReadFinal, record) -1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						
						seqContig.append(contig.contigSequence.subSequence(refLeft - 1, refRight));
					}else{//neg strain						
						int refLeft = positionOnRef(posReadEnd, record) + 1;		
						int refRight = record.refStart;
	
						if (posReadFinal <= record.readAlignmentEnd()){
							refLeft = positionOnRef(posReadFinal, record) + 1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						seqContig.append(Alphabet.DNA.complement(contig.contigSequence.subSequence(refRight - 1, refLeft)));
					}
				}
			}
			seqList.add(seqContig.toSequence());
			if(connection.readID.contains("twodimentional")) //only add 2D reads to calculate the consensus
				seqList.add(seqRead);
		}
		Sequence consensus = null;
		try {
			consensus = ErrorCorrection.consensusSequence(seqList, hashKey, "poa");
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("Can not generate consensus sequence!");
		}
		return consensus;
//		AlignmentRecord 	first=fromContig==firstContig?tAlign.clones():fAlign.clones(),
//							second=fromContig==firstContig?fAlign.clones():tAlign.clones();
//
//
//		SequenceBuilder builder = new SequenceBuilder(Alphabet.DNA16(), 1024*1024, "consensus");
//		
//		Sequence 	ligateToStart = first.contig.contigSequence.subSequence(first.refStart, first.refEnd),
//					ligateToEnd = second.contig.contigSequence.subSequence(second.refStart, second.refEnd);
//		first.readID = second.readID = -1;
//		first.strand=fromContig.getRelDir()>0;
//		second.strand = toContig.getRelDir()>0;
//		//        |--|------------|---->
//		//   -----|--|--        --|----|-----
//		first.readLength = second.readLength = ligateToStart.length() + consensus.length() + ligateToEnd.length();
//		
//		first.readStart = first.strand?0:ligateToStart.length()-1; 
//		first.readEnd = first.strand?ligateToStart.length()-1:0;
//		second.readStart = second.strand?ligateToStart.length()+consensus.length()-1:second.readLength-1; 
//		second.readEnd = second.strand?second.readLength-1:ligateToStart.length()+consensus.length()-1;
//		
//		builder.append(ligateToStart);		
//		builder.append(consensus);
//		builder.append(ligateToEnd);
//		ArrayList<AlignmentRecord> alignList = new ArrayList<AlignmentRecord>();
//		alignList.add(first); alignList.add(second);
//		ReadFilling rf=new ReadFilling(builder.toSequence(), alignList);
//		return new Connection(rf, first, second, ScaffoldVector.composition(second.contig.getVector(),ScaffoldVector.reverse(first.contig.getVector())));
//		
	}

	public class Connection implements Comparable<Connection>{
		ReadFilling read;
		String readID;
		int score;
		ScaffoldVector trans;
		AlignmentRecord firstAlignment, secondAlignment;

		int distanceOnRead = 0;

		Connection(ReadFilling mRead, AlignmentRecord a, AlignmentRecord b, ScaffoldVector trans){
			this.read = mRead;
			this.readID = a.readID;
			this.firstAlignment = a;
			this.secondAlignment = b;

			int aAlign = Math.abs(a.refStart - a.refEnd);
			int bAlign = Math.abs(b.refStart - b.refEnd);

			score = aAlign * bAlign / (aAlign  +bAlign);
			this.trans = trans;			

			distanceOnRead = Math.max(Math.min(b.readStart, b.readEnd) - Math.max(a.readEnd,a.readStart),
					Math.min(a.readStart, a.readEnd) - Math.max(b.readEnd,b.readStart));											
		}

		void display (){
			System.out.printf("[%6d %6d] -> [%6d %6d] : [%6d %6d] -> [%6d %6d] (%s) score=%d Read %s ==> %d [%d]\n", 
					firstAlignment.refStart, firstAlignment.refEnd, secondAlignment.refStart, secondAlignment.refEnd,
					firstAlignment.readStart, firstAlignment.readEnd, secondAlignment.readStart, secondAlignment.readEnd,
					trans.toString(),					
					score, read.readSequence.getName(),
					trans.distance(firstContig, secondContig),
					distanceOnRead);
		}

		/**
		 * Count the number of gaps (that is the number of bases that are not
		 * aligned to a contig) between 2 main contigs 
		 * @return
		 */
		public int gapsBetween(){
			int start = 0, end = read.readSequence.length();

			if (firstAlignment.readAlignmentStart() < secondAlignment.readAlignmentStart()){
				start = firstAlignment.readAlignmentEnd();
				end = secondAlignment.readAlignmentStart();
			}else{
				start = secondAlignment.readAlignmentEnd();
				end = firstAlignment.readAlignmentStart();
			}	

			if (start >= end)
				return 0;

			BitSet bitSet = new BitSet(end);
			bitSet.set(start, end);
			for (AlignmentRecord record:read.alignments){
				bitSet.clear(record.readAlignmentStart(), record.readAlignmentEnd());
			}

			return bitSet.cardinality();
		}

		public int filling(SequenceBuilder seqBuilder, JapsaAnnotation anno){
			Contig fromContig = firstContig;
			Contig toContig = secondContig;

			AlignmentRecord fromAlignment = firstAlignment,
							toAlignment = secondAlignment;

			ReadFilling readFilling = read;
			if ( fromContig.getRelDir()>0 != fromAlignment.strand){
				//swap
				readFilling = read.reverse();
				readFilling.sortAlignment();

				fromAlignment = fromAlignment.reverseRead();
				toAlignment = toAlignment.reverseRead();
			}
			//now readFilling is good to go
			int posReadEnd   = fromAlignment.readAlignmentEnd();
			int posReadFinal = toAlignment.readAlignmentStart();// I need as far as posReadFinal
			// locate the last position being extended...
			int lastExtendedPosition = posReadFinal;
			if(posReadEnd > posReadFinal -1 ){
				lastExtendedPosition = Math.min(posReadEnd,toAlignment.readAlignmentEnd());
				return positionOnRef(lastExtendedPosition, toAlignment);
			}
			if(seqBuilder == null)
				return positionOnRef(lastExtendedPosition, toAlignment);
			
			for (AlignmentRecord record:readFilling.alignments){
				Contig contig = record.contig;
				if (contig.getIndex() == fromContig.getIndex())
					continue;

				if (posReadEnd >= posReadFinal -1)
					//continue;//I can break here, but want to get portionUsed of other contigs
					break; //so do it!


				if (record.readAlignmentEnd() < posReadEnd)
					continue;				
			
				if (record.readAlignmentStart() > posReadEnd){
					//Really need to fill in using read information
					int newPosReadEnd = Math.min(posReadFinal - 1, record.readAlignmentStart() -1);
					if (newPosReadEnd > posReadEnd){
//						JapsaFeature feature = 
//								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() + newPosReadEnd - posReadEnd,
//										"CONTIG",readFilling.readSequence.getName(),'+',"");
//
//						//P=0 get the orignial read name and position
//						feature.addDesc(readFilling.readSequence.getName() + "+("+(posReadEnd + 1) +"," + newPosReadEnd+")");
//						anno.add(feature);
						seqBuilder.append(readFilling.readSequence.subSequence(posReadEnd, newPosReadEnd));
						posReadEnd = newPosReadEnd;
					}
					if (posReadEnd + 1 >= posReadFinal)
						continue;//Done

					//Now get information on the contig from start
					if (contig.getIndex() == toContig.getIndex())
						continue;//could break
					if (record.strand){
						int refLeft = record.refStart;
						int refRight = record.refEnd;

						if (posReadFinal <= record.readAlignmentEnd()){
							refRight = positionOnRef(posReadFinal, record) -1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						
//						JapsaFeature feature = 
//								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() + refRight - refLeft +1,
//										"CONTIG",contig.getName(),'+',"");
//						feature.addDesc(contig.getName() + "+("+(refLeft ) +"," + refRight+")");
//						anno.add(feature);
						
						seqBuilder.append(contig.contigSequence.subSequence(refLeft - 1, refRight));

					}else{//neg strain
						int refRight = record.refStart;
						int refLeft = record.refEnd;

						if (posReadFinal <= record.readAlignmentEnd()){
							refLeft = positionOnRef(posReadFinal, record) + 1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						
//						JapsaFeature feature = 
//								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() - refRight + refLeft +1,
//										"CONTIG",contig.getName(),'+',"");
//						feature.addDesc(contig.getName() + "-("+(refRight ) +"," + refLeft+")");
//						anno.add(feature);
						
						seqBuilder.append(Alphabet.DNA.complement(contig.contigSequence.subSequence(refRight - 1, refLeft)));

					}
				}//if record.readAlignmentStart() > posReadEnd
				else{//Now get information on the contig from start
					if (contig.getIndex() == toContig.getIndex())
						continue;//could break
					if (record.strand){
						int refLeft = positionOnRef(posReadEnd, record) + 1;						
						int refRight = record.refEnd;

						if (posReadFinal <= record.readAlignmentEnd()){
							refRight = positionOnRef(posReadFinal, record) -1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						
//						JapsaFeature feature = 
//								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() + refRight - refLeft +1,
//										"CONTIG",contig.getName(),'+',"");
//						feature.addDesc(contig.getName() + "+("+(refLeft ) +"," + refRight+")");
//						anno.add(feature);
						
						seqBuilder.append(contig.contigSequence.subSequence(refLeft - 1, refRight));
	
					}else{//neg strand						
						int refLeft = positionOnRef(posReadEnd, record) + 1;		
						int refRight = record.refStart;

						if (posReadFinal <= record.readAlignmentEnd()){
							refRight = positionOnRef(posReadFinal, record) + 1; 
							posReadEnd = posReadFinal -1;
						}else{
							posReadEnd = record.readAlignmentEnd();
						}
						
//						JapsaFeature feature = 
//								new JapsaFeature(seqBuilder.length() + 1, seqBuilder.length() - refRight + refLeft +1,
//										"CONTIG",contig.getName(),'+',"");
//						feature.addDesc(contig.getName() + "-("+(refRight ) +"," + refLeft+")");
//						anno.add(feature);
						
						seqBuilder.append(Alphabet.DNA.complement(contig.contigSequence.subSequence(refRight - 1, refLeft)));
				
					}
				}
			}	
			return positionOnRef(lastExtendedPosition, toAlignment);
		}
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(Connection o) {
			return o.score - score;
		}
	}
}
