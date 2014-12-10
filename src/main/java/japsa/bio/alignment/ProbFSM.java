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
 * 08/12/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.alignment;


import java.util.HashMap;
import java.util.Random;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.util.ByteArray;
import japsa.util.JapsaMath;

/**
 * Implementation of alignment using a probabilistic finite state machine
 * @author minhduc
 *
 */
public abstract class ProbFSM {	
	/**
	 * List of states: 
	 */
	MachineState [] states;

	Sequence mSeq;//sequence belong to the model
	double insEmissionCost = 2, changeEmissionCost = JapsaMath.log2(3);
	
	Alphabet alphabet = Alphabet.DNA();
	
	/**
	 * Reset all counts before each learning step
	 */
	public void resetCount(){
		for (int i = 0; i < states.length;i++){
			states[i].countCopy	= states[i].countMutate = states[i].countIns = states[i].countDel = 0;  
		}
	}

	/**
	 * Re-estimate parameters based on counts
	 */
	public void reEstimate(){
		int countM = 0, countC = 0;
		for (int i = 0; i < states.length;i++){		
			countM += states[i].countMutate;
			countC += states[i].countCopy ;			
			states[i].setTransitionProb(states[i].countMutate + states[i].countCopy + 1.0, states[i].countIns + 1.0, states[i].countDel + 1.0);
		}
		double probC = (countC + 1.0) / (countM + countC + 2.0);
		
		for (int i = 0; i < states.length;i++){
			states[i].setCopyProb(probC);
			System.out.printf("State %s: %3d %3d %3d %3d %8.4f %8.4f %8.4f %8.4f %8.4f\n", states[i].name, states[i].countCopy, states[i].countMutate, states[i].countIns, states[i].countDel, states[i].matchProb,  states[i].insProb,states[i].delProb, states[i].copyProb, states[i].changeProb);
		}
			

			
		
	}
	
	/**
	 * Show all the parameters of the machine
	 */
	
	public void showProb(){
		for (int i = 0; i < states.length;i++){			
			System.out.printf("Prob state %s : [%8.4f %8.4f %8.4f] [%8.4f %8.4f]\n", states[i].name, states[i].matchProb,  states[i].insProb,states[i].delProb, states[i].copyProb, states[i].changeProb);
			//System.out.printf("Cost state %s : %8.4f %8.4f %8.4f %8.4f %8.4f\n", states[i].name, states[i].matchCost,  states[i].insCost,states[i].delCost, states[i].copyCost, states[i].changeCost);
		}

		
	}
	
	

	/**
	 * Generate a sequence by the machine
	 * @param rnd
	 * @return
	 */
	public Sequence generate(Random rnd){		
		int indexSrc = 0;

		ByteArray byteArray = new ByteArray(mSeq.length() * 2);
		MachineState currentState = states[0];

		double toss;
		double cost = 0;
		while (indexSrc < mSeq.length()){
			toss = rnd.nextDouble();//chosing the transition
			//System.out.print(toss + " ==> ");

			if (toss < currentState.delProb && currentState.delState != null){
				//System.out.println("Gen D " + indexSrc);
				cost += currentState.delCost;
				indexSrc ++;
			}else if (toss < currentState.insProb + currentState.delProb && currentState.insState != null){
				//System.out.println("Gen I " + indexSrc);
				//Insertion
				cost += currentState.insCost;

				//cost of emitting the inserted base
				byteArray.add((byte)rnd.nextInt(4));
				cost += insEmissionCost;				
			}else{//match
				cost += currentState.matchCost;
				toss = rnd.nextDouble();//copy or change
				byte base = mSeq.getBase(indexSrc);
				if (toss < currentState.copyProb){
					//copy
					//System.out.println("Gen C " + indexSrc);
					cost += currentState.copyCost;
					byteArray.add(base);		
				}else{
					//System.out.println("Gen M " + indexSrc);
					//change
					cost += currentState.changeCost;

					//need to toss again
					byteArray.add((byte)((base + 1 + rnd.nextInt(3)) % 4));
					cost += changeEmissionCost;								
				}
				indexSrc ++;
			}
		}

		System.out.println("Generating cost " + cost);
		return new Sequence(alphabet, byteArray, "gene");
	}

	/**
	 * Update counts based on the path of best alignment (backward pass)
	 * @param emiss
	 * @return
	 */
	public int updateCount(Emission emiss){
		int countEmis = 0;
		while (true){

			Emission bwdEmission = emiss.bwdEmission;
			if (bwdEmission == null)
				break;

			switch (emiss.type){			
			case INSERTION:bwdEmission.toState.countIns ++;break;				
			case DELETION:bwdEmission.toState.countDel ++;break;			
			case COPY:bwdEmission.toState.countCopy ++;break;			
			case MUTATE:bwdEmission.toState.countMutate ++;break;				
			}
			countEmis ++;
			emiss = bwdEmission;			
		}
		return countEmis;
	}

	/**
	 * Backward pass: only print out the path in reverse. Consider updateCount
	 * for learning
	 * @param emiss
	 * @return
	 */
	public int backward(Emission emiss){
		int countEmis = 0;
		while (true){

			Emission bwdEmission = emiss.bwdEmission;
			if (bwdEmission == null)
				break;

			switch (emiss.type){			
			case INSERTION:System.out.println(bwdEmission.toState.name + " I");break;				
			case DELETION:System.out.println(bwdEmission.toState.name + " D");break;
			case COPY:System.out.println(bwdEmission.toState.name + " C");break;
			case MUTATE:System.out.println(bwdEmission.toState.name + " M");break;		
			}
			countEmis ++;
			emiss = bwdEmission;			
		}
		return countEmis;
	}

	/************************************************************************/
	/**
	 * Forward pass: find the best path to align a sequence
	 * @param genSeq
	 * @return
	 */
	public Emission align(Sequence genSeq){
		//return state
		Emission retEmission = new Emission(states[0], mSeq.length()-1, genSeq.length() -1);
		retEmission.myCost = genSeq.length() * (insEmissionCost + 4);


		Emission currentEmission, finalEmission;//current pointer and last pointer on the linked-list		
		currentEmission = finalEmission = new Emission(states[0],-1,-1);
		currentEmission.myCost = 0;		

		HashMap<String, Emission> hash = new HashMap<String, Emission>();

		while (currentEmission != null){//linked list not exhausted
			if (currentEmission.gPos >= genSeq.length() - 1){
				//done generating genSeq
				if (currentEmission.myCost < retEmission.myCost){
					retEmission = currentEmission;					
				}
			}else if (currentEmission.myCost < retEmission.myCost){
				String hashKey;
				Emission nextEmission;
				double cost;

				//1. consider deletion if profile has something to offer			
				if (currentEmission.mPos + 1 < mSeq.length() && currentEmission.toState.delState != null){
					cost = currentEmission.myCost + currentEmission.toState.delCost;

					hashKey= Emission.hashKey(currentEmission.toState.delState.name, currentEmission.mPos + 1, currentEmission.gPos);
					nextEmission = hash.get(hashKey);
					if (nextEmission == null){
						nextEmission = new Emission(currentEmission.toState.delState, currentEmission.mPos + 1, currentEmission.gPos);
						nextEmission.type = EmissionType.DELETION;

						nextEmission.myCost = cost;
						hash.put(hashKey, nextEmission);

						finalEmission.next = nextEmission;
						finalEmission = nextEmission;
						nextEmission.bwdEmission = currentEmission;
					}else{
						if (nextEmission.myCost > cost){
							nextEmission.myCost = cost;
							nextEmission.bwdEmission = currentEmission;
							nextEmission.type = EmissionType.DELETION;
						}
					}//else - if nextstate != null
				}//if mPos

				//2. insertion			
				if (currentEmission.gPos + 1 < genSeq.length() && currentEmission.toState.insState != null){
					cost = currentEmission.myCost + currentEmission.toState.insCost + this.insEmissionCost;

					hashKey= Emission.hashKey(currentEmission.toState.insState.name, currentEmission.mPos, currentEmission.gPos + 1);
					nextEmission = hash.get(hashKey);
					if (nextEmission == null){
						nextEmission = new Emission(currentEmission.toState.insState, currentEmission.mPos, currentEmission.gPos + 1);
						nextEmission.myCost = cost;
						nextEmission.type = EmissionType.INSERTION;
						hash.put(hashKey, nextEmission);

						finalEmission.next = nextEmission;
						finalEmission = nextEmission;
						nextEmission.bwdEmission = currentEmission;

					}else{
						if (nextEmission.myCost > cost){
							nextEmission.myCost = cost;
							nextEmission.bwdEmission = currentEmission;
							nextEmission.type = EmissionType.INSERTION;
						}
					}
				}

				//3.Match
				if (currentEmission.gPos + 1 < genSeq.length() && currentEmission.mPos + 1 < mSeq.length()){
					EmissionType type = EmissionType.COPY;

					if (mSeq.getBase(currentEmission.mPos + 1) == genSeq.getBase(currentEmission.gPos + 1)){
						cost = currentEmission.myCost + currentEmission.toState.matchCost + currentEmission.toState.copyCost; 					
					}else{
						cost = currentEmission.myCost + currentEmission.toState.matchCost + currentEmission.toState.changeCost + this.changeEmissionCost;
						type = EmissionType.MUTATE;
					}

					if (cost < retEmission.myCost){
						hashKey= Emission.hashKey(currentEmission.toState.matchState.name, currentEmission.mPos + 1, currentEmission.gPos +1);
						nextEmission = hash.get(hashKey);
						if (nextEmission == null){
							nextEmission = new Emission(currentEmission.toState.matchState, currentEmission.mPos + 1, currentEmission.gPos +1);
							nextEmission.type = type;

							nextEmission.myCost = cost;
							hash.put(hashKey, nextEmission);
							finalEmission.next = nextEmission;
							finalEmission = nextEmission;
							nextEmission.bwdEmission = currentEmission;
						}else{
							if (nextEmission.myCost > cost){
								nextEmission.myCost = cost;
								nextEmission.bwdEmission = currentEmission;
								nextEmission.type = type;
							}
						}
					}
				}//match
			}

			//helping GC to gabbabe collect current state
			Emission tmp = currentEmission.next;
			currentEmission.next = null;
			currentEmission = tmp;
		}
		return retEmission;
	}
	/********************************************************************/

	public static enum EmissionType {COPY, MUTATE,  INSERTION, DELETION};
	public static class Emission{
		//This pointer is for implementing dynamic programming using a linked-list
		Emission next = null;		

		EmissionType type;

		public int gPos, mPos;

		MachineState toState;

		//These are for backward/forward
		//public Emission fwdEmission = null;
		public Emission bwdEmission = null;

		public double myCost;

		String hashKey;

		public Emission(MachineState state, int mP, int gP){
			toState = state;
			mPos = mP;
			gPos = gP;
			myCost = Double.MAX_VALUE;
		}

		static public String hashKey(String name, int mPos, int gPos){
			return name + "_" + mPos + "_" + gPos;

		}
		public double getScore(){
			return myCost;
		}			
	}

	public static class MachineState{	
		double copyProb = 0.9, changeProb = 0.1;
		double matchProb = 0.8, insProb = 0.14, delProb = 0.06;

		double  copyCost, changeCost, matchCost,  insCost,	delCost;

		//For training
		int countIns = 0, countDel = 0, countCopy = 0, countMutate = 0;

		MachineState insState, delState, matchState;
		String name;

		public MachineState(String name){
			this.name = name;

		}

		public void setTransitionProb(double mP, double iP, double dP){
			double sum = mP + iP + dP;
			matchProb = mP/ sum;
			insProb   = iP/sum;
			delProb   = dP /sum;

			matchCost = (matchProb > 0)? -JapsaMath.log2(matchProb):Double.MAX_VALUE;
			insCost = (insProb > 0)? -JapsaMath.log2(insProb):Double.MAX_VALUE;
			delCost = (delProb > 0)? -JapsaMath.log2(delProb):Double.MAX_VALUE;
		}
		public void setCopyProb(double cP){
			copyProb = cP;
			changeProb = 1 - copyProb;

			copyCost =  (copyProb > 0)? -JapsaMath.log2(copyProb):Double.MAX_VALUE;
			changeCost =  (changeProb > 0)? -JapsaMath.log2(changeProb):Double.MAX_VALUE;
		}
	}


	public static class ProbThreeSM extends ProbFSM{
		public ProbThreeSM(Sequence seq){
			states = new MachineState[3];

			states[0] = new MachineState("S");
			states[1] = new MachineState("I");
			states[2] = new MachineState("D");

			states[0].matchState = states[0];
			states[0].insState = states[1];			
			states[0].delState = states[2];
			states[0].setCopyProb(0.9);
			states[0].setTransitionProb(.8, 0.1, 0.1);


			states[1].matchState = states[0];
			states[1].insState = states[1];
			states[1].delState = null;
			states[1].setCopyProb(0.9);
			states[1].setTransitionProb(.8, 0.2, 0);
			

			states[2].matchState = states[0];
			states[2].insState = null;
			states[2].delState = states[2];
			
			states[2].setCopyProb(0.9);
			states[2].setTransitionProb(.8, 0, 0.2);
			this.mSeq = seq;
		}
		
		public MachineState getMatState(){
			return states[0];
		}
		
		public MachineState getInsState(){
			return states[1];
		}
		
		public MachineState getDelState(){
			return states[2];
		}
	}

	public static class ProbOneSM extends ProbFSM{
		public ProbOneSM(Sequence seq){
			states = new MachineState[1];
			states[0] = new MachineState("S");
			states[0].insState = states[0].delState = states[0].matchState = states[0];			
			
			states[0].setCopyProb(0.9);
			states[0].setTransitionProb(.85, 0.07, 0.08);
			
			this.mSeq = seq;
		}
		
		public MachineState getState(){
			return states[0];
		}
	}



}
