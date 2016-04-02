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

package japsa.bio.hts.scaffold;

import java.util.ArrayList;
import java.util.Collections;

import japsa.seq.Sequence;
import japsa.seq.JapsaFeature;

public class Contig{
	int index;
	ScaffoldVector myVector;//relative position to the head contig of my scaffold	
	Sequence contigSequence;//the sequence of the contig	
	double   coverage = 1.0;
	int head = -1; //point to the index of its head contig in the scaffold 
	double prevScore=0, nextScore=0;
	boolean isCircular = false;
	//for depth first search
	ArrayList<ContigBridge> bridges;	
	//for annotation
	ArrayList<JapsaFeature> oriRep,				//origin of replication: indicator of plasmid for bacteria
							insertSeq,			// IS
							resistanceGenes;	//list of antibiotic resistance genes found in this contig
		
	public Contig(int index, Sequence seq){
		this.index = index;
		contigSequence = seq;
		myVector = new ScaffoldVector(0,1);
		bridges = new ArrayList<ContigBridge>();
		usedRanges = new ArrayList<Range>();
		
		oriRep = new ArrayList<JapsaFeature>();
		insertSeq = new ArrayList<JapsaFeature>();
		resistanceGenes = new ArrayList<JapsaFeature>();
	}

	public Contig clone(){
		Contig ctg = new Contig(this.index, this.contigSequence);
		ctg.coverage = coverage;
		ctg.bridges = this.bridges;
		ctg.head = this.head; //update later
		ctg.isCircular = this.isCircular;
		ctg.usedRanges = this.usedRanges;

		ctg.oriRep = this.oriRep;
		ctg.insertSeq = this.insertSeq;
		ctg.resistanceGenes = this.resistanceGenes;
		return ctg;
	}
	// Get features in an interval of contig
	public ArrayList<JapsaFeature> getFeatures(ArrayList<JapsaFeature> features, int start, int end){
		
		ArrayList<JapsaFeature> remainFeatures = new ArrayList<JapsaFeature>();
		boolean isReverse= (start>end)?true:false;
		for(JapsaFeature feature:features){
			JapsaFeature cutFeature=feature.cloneFeature();
			int fstart = feature.getStart(),
				fend = feature.getEnd();
			
			//find overlap
			if(Integer.signum(fstart-start)*Integer.signum(fstart-end) <= 0){
				if(Integer.signum(fend-start)*Integer.signum(fend-end) > 0){
					fend = (Math.abs(fend-start) < Math.abs(fend-end))?start:end;
				}
			}else{
				fstart = (Math.abs(fstart-start) < Math.abs(fstart-end))?start:end;
				if(Integer.signum(start-fend)*Integer.signum(start-fstart) <= 0 && Integer.signum(end-fend)*Integer.signum(end-fstart) <= 0)
					fend = (Math.abs(fend-start) < Math.abs(fend-end))?start:end;
				else if(Integer.signum(start-fend)*Integer.signum(start-fstart) > 0 && Integer.signum(end-fend)*Integer.signum(end-fstart) > 0)
					continue;
			}
			//if the contig is reversed complement
			if(isReverse){
				int ostart = fstart;
				fstart= this.length() - fend;
				fend = this.length() - ostart;
				if(cutFeature.getStrand() == '+')
					cutFeature.setStrand('-');
				else
					cutFeature.setStrand('+');
			}
				
			cutFeature.setStart(fstart);
			cutFeature.setEnd(fend);
			double cutRate=(float) Math.abs(cutFeature.getLength())/Math.abs(feature.getLength());
			if(cutRate > .9){
				cutFeature.setScore(feature.getScore()*cutRate);
				remainFeatures.add(cutFeature);
				
			}
		}
		
		return remainFeatures;
	}

	
	public String getName(){
		return contigSequence.getName();
	}
	
	public int getIndex(){
		return index;
	}
	//actually a backward composite
	public void composite(ScaffoldVector aVector){
		myVector = ScaffoldVector.composition(myVector, aVector);		
	}	
	/**
	 * Relative position to the head of the scaffold
	 * @return
	 */
	public int getRelPos(){
		return myVector.magnitude;
	}	
	
	public int getRelDir(){
		return myVector.direction;
	}
		
	/**
	 * Get the left most position if transpose by vector trans
	 * @return
	 */
	public int leftMost(ScaffoldVector trans){
		return trans.magnitude - ((trans.direction > 0)?0:length()); 
	}
	
	/**
	 * Get the right most position if transpose by vector trans
	 * @return
	 */
	public int rightMost(ScaffoldVector trans){
		return trans.magnitude + ((trans.direction > 0)?length():0); 
	}
	
	/**
	 * Get the left most position
	 * @return
	 */
	public int leftMost(){
		return leftMost(myVector); 
	}
	
	/**
	 * Get the right most position
	 * @return
	 */
	public int rightMost(){
		return rightMost(myVector); 
	}
	
	
	public ScaffoldVector getVector(){
		return myVector;
	}
	
	public int length(){
		return contigSequence.length();
	}		
	
	public double getCoverage(){
		return coverage;
	}
	public void setCoverage(double cov){
		coverage =  cov;
	}
	public String toString(){
		return new String(" contig" + getIndex());
	}
	////////////////for tracing the used part////////////////////
	ArrayList<Range> usedRanges;
	class Range implements Comparable<Range> {
		int start, end, score;
		Range(){
			start = end = score = 0;
		}
		Range(int start, int end, int score){
			this.start = start<end?start:end;
			this.end = start+end-this.start;
			this.score = score;
		}
		public int getLen(){
			return Math.abs(end-start)+1;
		}
		public String toString(){
			return new String(start + " --> " + end + ": " + score);
		}
		@Override
		public int compareTo(Range rg) {
			// TODO Auto-generated method stub
			if(this.start!=rg.start)
				return (this.start-rg.start);
			else if(this.end != rg.end)
				return (this.end-rg.end);
			else
				return (this.score-rg.score);
		}
	}
	public void addRange(int start, int end, int score){
		usedRanges.add(new Range(start,end,score));
	}
	public void display(){
		Collections.sort(usedRanges);
		System.out.println("Contig " + this.getName());
		for(Range rg:usedRanges)
			System.out.println("used " + rg);
		
		int prevEnd = 0, minLen = 100;
		for (Range rg:usedRanges){
			if(rg.start > prevEnd + minLen){
				System.out.println("\tuncovered: " + prevEnd + " --> " + rg.start + " ( " + (rg.start-prevEnd+1) +" )");
			}
			if(prevEnd < rg.end)
				prevEnd = rg.end;
		}
		if (prevEnd + minLen < length()-1)
			System.out.println("\tuncovered: " + prevEnd + " --> " + (length()-1) + " ( " + (length()-prevEnd) +" )");
	}
}