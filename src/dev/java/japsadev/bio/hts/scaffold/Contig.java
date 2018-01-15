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

import java.util.ArrayList;
import java.util.stream.IntStream;

import japsa.seq.Sequence;
import japsa.util.Logging;
import japsa.seq.JapsaFeature;
import japsa.util.Range;

public class Contig{
	int index;
	ScaffoldVector myVector;//relative position to the head contig of my scaffold	
	Sequence contigSequence;//the sequence of the contig	
	double coverage = 1.0;
	private double aStats=0;
	private boolean estRepeatFlag;
	int head = -1; //point to the index of its head contig in the scaffold 
	double prevScore=0, nextScore=0;
	ArrayList<Range> lowConfidentRegions;
	double cirProb = -1.0; //measure how likely the contig itself is circular
	//int[] isMapped; //which bases is mapped by any long reads
	//for annotation
	ArrayList<JapsaFeature> genes,				//genes list
							oriRep,				//origin of replication: indicator of plasmid for bacteria
							insertSeq,			//Insertion Sequence
							resistanceGenes;	//list of antibiotic resistance genes found in this contig
	//a contig is composed of edges from assembly graph

	static Graph asGraph=null;
	public static void setGraph(Graph g){
		asGraph=g;
	}
	public static boolean hasGraph(){
		return asGraph!=null;
	}

	ArrayList<Path> paths;
	
	public Contig(int index, Sequence seq){
		this.index = index;
		contigSequence = seq;
//		isMapped = new int[seq.length()];
		
		myVector = new ScaffoldVector(0,1);
		
		genes = new ArrayList<JapsaFeature>();
		oriRep = new ArrayList<JapsaFeature>();
		insertSeq = new ArrayList<JapsaFeature>();
		resistanceGenes = new ArrayList<JapsaFeature>();
		lowConfidentRegions = new ArrayList<Range>();
		paths = new ArrayList<Path>();
	}

	
	public Contig clone(){
		Contig ctg = new Contig(this.index, this.contigSequence);
		ctg.coverage = coverage;

		ctg.head = this.head; //update later
		ctg.cirProb = this.cirProb;
//		ctg.isMapped = this.isMapped;
		
		ctg.genes = this.genes;
		ctg.oriRep = this.oriRep;
		ctg.insertSeq = this.insertSeq;
		ctg.resistanceGenes = this.resistanceGenes;
		ctg.lowConfidentRegions = this.lowConfidentRegions;
		ctg.paths = new ArrayList<Path>();
		for(Path p:paths)
			ctg.paths.add(p);
		
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

	//get the SPAdes name (out of MicroManage name maybe)
	public String getName(){
		return contigSequence.getName();
	}
	public String getDesc(){
		return contigSequence.getDesc();
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
	
	public boolean isCircular(){
		return !ScaffoldGraph.eukaryotic && (cirProb > 1);
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

	public double getAstatistics(){
		return aStats;
	}
	public void setAstatistics(double stats){
		aStats =  stats;
	}
	
	public boolean getRepeatFlag(){
		return estRepeatFlag;
	}
	public void setRepeatFlag(boolean flag){
		estRepeatFlag = flag;
	}
//	public boolean isMapped(){
//		int sum = IntStream.of(isMapped).sum();
//		return ((double)sum/length()) > .8;
//	}
	/*
	 * Operators related to Path
	 */
	public ArrayList<Path> getPaths(){
		return paths;
	}
	public void setPath(Path path){
		this.paths.add(path);
	}
	

	public String toString(){
		return new String(" contig" + getIndex());
	}
	
	public void addLowConfidentRegion(Range r){
		if(lowConfidentRegions.size()==0){
			lowConfidentRegions.add(r);
			return;
		}
		int index=0;
		
		while(lowConfidentRegions.size()>index){
			if(lowConfidentRegions.get(index).getRight() >= r.getLeft())
				break;
			index++;
		}
		if(index==lowConfidentRegions.size())
			lowConfidentRegions.add(r);
		else {
			Range cur = lowConfidentRegions.get(index++);
			int start=index;
			if(cur.merge(r)){
				while(lowConfidentRegions.size()>index){
					if(!cur.merge(lowConfidentRegions.get(index)))
						break;
					index++;
				}
				for(int i=start;i<index;i++)
					lowConfidentRegions.remove(start);
			}
			else
				lowConfidentRegions.add(index-1, r);

		}
			
		
	}
	//calculate the intersection between a range and the list of contig's low confident regions
	public int countLowBases(Range r){
		int retval=0;
		for(Range range:lowConfidentRegions){
			if(range.getRight()<r.getLeft())
				continue;
			if(range.getLeft()>r.getRight())
				break;
			retval+=Math.min(r.getRight(), range.getRight())-Math.max(r.getLeft(), range.getLeft());
			
		}
		return retval;
	}
	public String displayLowConfidentRegions(){
		String retval="";
		for(Range r:lowConfidentRegions){
			retval+=r.toString()+" ; ";
		}
		return retval;
	}
	
}