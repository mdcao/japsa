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
 * 10/08/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.alignment;

import java.io.IOException;
import java.util.HashMap;
import java.util.Random;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.util.ByteArray;
import japsa.util.JapsaMath;


/**
 * This is a one-state finite automata
 * @author minhduc
 *
 */
public class ProfileDP {	
	// cost from a match state
	double  matProb,
	insProb,
	delProb;

	double  matchProb,misMatchProb;

	private double 	matCost, delCost, insCost;

	public double getMatCost() {
		return matCost;
	}

	public double getDelCost() {
		return delCost;
	}

	public double getInsCost() {
		return insCost;
	}

	public double getMatchCost() {
		return matchCost;
	}

	public double getMisMatchCost() {
		return misMatchCost;
	}

	private double  matchCost, misMatchCost;

	Alphabet alphabet = Alphabet.DNA6();
	Sequence profileSeq;
	int repStart, repEnd;	

	/**
	 * @param args
	 * retStart: inclusive
	 * repEnd:incluive
	 *  
	 */
	public ProfileDP(Sequence seq, int repStart, int repEnd){
		profileSeq    = seq;
		this.repEnd   = repEnd;
		this.repStart = repStart;
		setTransitionProbability(0.95,0.025,0.025);
		setMatchProbability(0.85/0.95);
	}

	/**************************************************/
	public void setTransitionProbability(double matP, double insP, double delP){		

		double sum = matP + insP + delP;
		matProb = matP / sum;
		insProb = insP / sum;
		delProb = delP / sum;

		matCost = -JapsaMath.log2(matProb);		
		insCost = -JapsaMath.log2(insProb) + 2;
		delCost = -JapsaMath.log2(delProb);
	}

	public void setMatchProbability(double matP){
		matchProb = matP;
		misMatchProb = (1.0 - matchProb);

		matchCost  = - JapsaMath.log2(matchProb);
		misMatchCost = -JapsaMath.log2(misMatchProb)- JapsaMath.log2(1.0/3.0);;
	}

	public int getProfileLength(){
		return profileSeq.length();
	}
	/**************************************************/

	public EmissionState align(Sequence seq){
		//return state
		EmissionState retState = new EmissionState(profileSeq.length() -1 , seq.length()-1, 0);
		retState.score = seq.length() * (insCost + 2);

		EmissionState currentState, lastState;
		currentState = lastState = new EmissionState(-1,-1,0);
		currentState.score = 0;


		HashMap<String, EmissionState> hash = new HashMap<String, EmissionState>();		
		while (currentState != null){
			if (currentState.seqPos >= seq.length() - 1){
				if (currentState.score < retState.score){
					retState = currentState;
					currentState = currentState.next;					
					continue;
				}
			}

			String hashKey;
			EmissionState nextState;
			double cost;

			int iterAdvance = 0;				
			//if it is about the enter the repeat
			if (currentState.profilePos + 1 == repStart){
				iterAdvance = 1;
			}

			//1. consider deletion if profile has something to offer			
			if (currentState.profilePos + 1 < profileSeq.length()){
				cost = currentState.score + delCost;

				if (cost < retState.score){ 
					hashKey= EmissionState.hashKey(currentState.seqPos, currentState.profilePos + 1, currentState.iter + iterAdvance);
					nextState = hash.get(hashKey);
					if (nextState == null){
						nextState = new EmissionState(currentState.seqPos, currentState.profilePos + 1, currentState.iter + iterAdvance);
						nextState.score = cost;
						hash.put(hashKey, nextState);
						lastState.next = nextState;
						lastState = nextState;
						nextState.bwdState = currentState;

						nextState.countDel = currentState.countDel + 1;
						nextState.countIns = currentState.countIns;
						nextState.countMG =  currentState.countMG;
						nextState.countMB =  currentState.countMB;

					}else{
						if (nextState.score > cost){
							nextState.score = cost;
							nextState.bwdState = currentState;

							nextState.countDel = currentState.countDel + 1;
							nextState.countIns = currentState.countIns;
							nextState.countMG =  currentState.countMG;
							nextState.countMB =  currentState.countMB;
						}
					}//else - if nextstate != null
				}//if cost
			}//if profile

			//2. insertion			
			if (currentState.seqPos + 1 < seq.length()){
				cost = currentState.score + insCost;
				if (cost < retState.score){
					//note: this does not advance on the profile, thus no need to add iterAdvance
					hashKey= EmissionState.hashKey(currentState.seqPos + 1, currentState.profilePos, currentState.iter);
					nextState = hash.get(hashKey);
					if (nextState == null){
						nextState = new EmissionState(currentState.seqPos + 1, currentState.profilePos, currentState.iter);
						nextState.score = cost;
						hash.put(hashKey, nextState);
						lastState.next = nextState;
						lastState = nextState;
						nextState.bwdState = currentState;

						nextState.countDel = currentState.countDel;
						nextState.countIns = currentState.countIns + 1;
						nextState.countMG =  currentState.countMG;
						nextState.countMB =  currentState.countMB;
					}else{
						if (nextState.score > cost){
							nextState.score = cost;
							nextState.bwdState = currentState;

							nextState.countDel = currentState.countDel;
							nextState.countIns = currentState.countIns + 1;
							nextState.countMG =  currentState.countMG;
							nextState.countMB =  currentState.countMB;
						}
					}
				}
			}

			//3.Match
			if (currentState.seqPos + 1 < seq.length() && currentState.profilePos + 1 < profileSeq.length()){
				cost = currentState.score + matCost + (seq.getBase(currentState.seqPos + 1) == profileSeq.getBase(currentState.profilePos + 1)? matchCost:misMatchCost);
				if (cost < retState.score){
					hashKey= EmissionState.hashKey(currentState.seqPos + 1, currentState.profilePos +1, currentState.iter + iterAdvance);
					nextState = hash.get(hashKey);
					if (nextState == null){
						nextState = new EmissionState(currentState.seqPos + 1, currentState.profilePos+1, currentState.iter + iterAdvance);
						nextState.score = cost;
						hash.put(hashKey, nextState);
						lastState.next = nextState;
						lastState = nextState;
						nextState.bwdState = currentState;

						nextState.countDel = currentState.countDel;
						nextState.countIns = currentState.countIns;
						nextState.countMG =  currentState.countMG;
						nextState.countMB =  currentState.countMB;

						if (seq.getBase(currentState.seqPos + 1) == profileSeq.getBase(currentState.profilePos + 1)){
							nextState.countMG ++;
						}else
							nextState.countMB ++;							


					}else{
						if (nextState.score > cost){
							nextState.score = cost;
							nextState.bwdState = currentState;

							nextState.countDel = currentState.countDel;
							nextState.countIns = currentState.countIns;
							nextState.countMG =  currentState.countMG;
							nextState.countMB =  currentState.countMB;

							if (seq.getBase(currentState.seqPos + 1) == profileSeq.getBase(currentState.profilePos + 1)){
								nextState.countMG ++;
							}else
								nextState.countMB ++;
						}
					}
				}
			}//match


			//Consider jumping to the beginning of the rep
			//Note I dont need iterAdvance here (it is 0 anyway)
			if (currentState.profilePos == repEnd){
				cost = currentState.score + delCost;
				if (cost < retState.score){ 
					hashKey= EmissionState.hashKey(currentState.seqPos, repStart, currentState.iter + 1);
					nextState = hash.get(hashKey);
					if (nextState == null){
						nextState = new EmissionState(currentState.seqPos, repStart, currentState.iter + 1);
						nextState.score = cost;
						hash.put(hashKey, nextState);
						lastState.next = nextState;
						lastState = nextState;
						nextState.bwdState = currentState;

						nextState.countDel = currentState.countDel + 1;
						nextState.countIns = currentState.countIns;
						nextState.countMG =  currentState.countMG;
						nextState.countMB =  currentState.countMB;						
					}else{
						if (nextState.score > cost){
							nextState.score = cost;
							nextState.bwdState = currentState;

							nextState.countDel = currentState.countDel + 1;
							nextState.countIns = currentState.countIns;
							nextState.countMG =  currentState.countMG;
							nextState.countMB =  currentState.countMB;
						}
					}
				}//if (del)		

				//3.Match
				if (currentState.seqPos + 1 < seq.length() && currentState.profilePos + 1 < profileSeq.length()){
					cost = currentState.score + matCost + (seq.getBase(currentState.seqPos + 1) == profileSeq.getBase(repStart)? matchCost:misMatchCost);
					if (cost < retState.score){
						hashKey= EmissionState.hashKey(currentState.seqPos + 1, repStart, currentState.iter + 1);
						nextState = hash.get(hashKey);
						if (nextState == null){
							nextState = new EmissionState(currentState.seqPos + 1, repStart, currentState.iter + 1);
							nextState.score = cost;
							hash.put(hashKey, nextState);
							lastState.next = nextState;
							lastState = nextState;
							nextState.bwdState = currentState;


							nextState.countDel = currentState.countDel;
							nextState.countIns = currentState.countIns;
							nextState.countMG =  currentState.countMG;
							nextState.countMB =  currentState.countMB;

							if (seq.getBase(currentState.seqPos + 1) == profileSeq.getBase(repStart)){
								nextState.countMG ++;
							}else
								nextState.countMB ++;				

						}else{
							if (nextState.score > cost){
								nextState.score = cost;
								nextState.bwdState = currentState;

								nextState.countDel = currentState.countDel;
								nextState.countIns = currentState.countIns;
								nextState.countMG =  currentState.countMG;
								nextState.countMB =  currentState.countMB;

								if (seq.getBase(currentState.seqPos + 1) == profileSeq.getBase(repStart)){
									nextState.countMG ++;
								}else
									nextState.countMB ++;
							}
						}
					}
				}
			}
			EmissionState tmp = currentState.next;

			//helping GC to gabbabe collect current state			
			hash.remove(EmissionState.hashKey(currentState.seqPos, currentState.profilePos, currentState.iter));			
			currentState.next = null;
			currentState = tmp;
		}
		//System.out.printf("Consider %d states\n",complexity);

		System.out.printf("Estimate: %2d %3d %3d %3d %3d %8.4f %8.4f\n",retState.iter, retState.countMG, retState.countMB, 
				retState.countIns, retState.countDel, retState.score,
				retState.countMG * (matCost + matchCost) + retState.countMB *(matCost + misMatchCost) + 
				retState.countIns * insCost + retState.countDel * delCost);	return retState;//bestScore;
	}

	public static class EmissionState{
		int countDel = 0, countIns = 0, countMG = 0, countMB = 0;

		EmissionState next = null;
		public int seqPos, profilePos, iter;

		public EmissionState bwdState = null;
		public double score;

		String hashKey;
		//		int count = 0;//how many in

		public EmissionState(int sPos, int pPos, int i){
			seqPos = sPos;
			profilePos = pPos;
			iter = i;
			score = Double.MAX_VALUE;			
			//hashKey = hashKey(gPos, mPos, iter);
		}

		static public String hashKey(int sPos, int pPos, int iter){
			return String.format("%5d  %4d  %2d",sPos, pPos, iter);			
		}
		public double getScore(){
			return score;
		}
		public int getIter(){
			return iter;
		}

		/**
		 * @return the countDel
		 */
		public int getCountDel() {
			return countDel;
		}

		/**
		 * @return the countIns
		 */
		public int getCountIns() {
			return countIns;
		}

		/**
		 * @return the countMG
		 */
		public int getCountMG() {
			return countMG;
		}

		/**
		 * @return the countMB
		 */
		public int getCountMB() {
			return countMB;
		}		
	}

	public Sequence generate(int iter, Random rnd) throws IOException{
		ByteArray bArray = new ByteArray();
		int profilePos = 0;		

		double cost = 0;
		byte nuc;
		double prob;
		int myIter = 0;
		int countMG = 0, countMB = 0, countDel = 0, countIns = 0;
		while (profilePos < this.profileSeq.length()){			
			prob = rnd.nextDouble();
			if (prob < delProb){
				cost += delCost;
				nuc = this.profileSeq.getBase(profilePos);
				//datOutGen.print("D " + mPos + " " + bArray.size() + " " + alphabet.int2char(nuc)+"\n");
				profilePos ++;
				countDel ++;				
			}else if (prob < insProb + delProb){
				cost += insCost;
				nuc = (byte)rnd.nextInt(4);
				//datOutGen.print("I " + mPos + " " + bArray.size() + " " + alphabet.int2char(nuc)+"\n");
				bArray.add(nuc);
				countIns ++;
			}else{
				//match
				prob = rnd.nextDouble();
				byte nextByte = this.profileSeq.getBase(profilePos);
				//System.out.println(mPos + " " + bArray.size());

				if (prob < matchProb){
					cost += matCost + matchCost;
					//	System.out.println(mPos + " " + bArray.size());
					nuc = this.profileSeq.getBase(profilePos);
					//datOutGen.print("= " + mPos + " " + bArray.size() + " " + alphabet.int2char(nuc)+'\n');
					bArray.add(nextByte);
					countMG ++;					
				}else{
					cost += matCost + misMatchCost;
					nuc  = (byte) ((nextByte + rnd.nextInt(3) + 1) % 4);
					//datOutGen.print("X " + mPos + " " + bArray.size() + " " + alphabet.int2char(nextByte) + " " + alphabet.int2char(nuc)+'\n');
					//System.out.println("B " + bArray.size() + " " + mPos + " " + this.profileSeq.getBase(mPos) + " " + nextByte);
					bArray.add(nuc);
					countMB ++;					
				}
				profilePos ++;
			}

			if (profilePos == this.repEnd && myIter < iter){
				profilePos = this.repStart;
				myIter ++;
			}
			//System.out.println();
		}		
		//datOutGen.print("GEN: " + countMG + " " +countMB + " " + countIns + " " + countDel + " " + cost + "\n");


		System.out.printf("Generate: %2d %3d %3d %3d %3d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n", myIter, countMG, countMB, countIns, countDel, cost,
				countMG * (matCost + matchCost) + countMB *(matCost + misMatchCost) + 
				countIns * insCost + countDel * delCost,
				matProb, insProb, delProb, matchProb, misMatchProb);
		return new Sequence(Alphabet.DNA4(), bArray, "gen");		
	}
}
