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
 * 31/01/2017 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsa.bio.alignment.ppfsm;

import java.util.ArrayList;
import java.util.HashMap;

import japsa.bio.alignment.ppfsm.state.MachineState;
import japsa.bio.alignment.ppfsm.transition.Transition;
import japsa.seq.Sequence;
import japsa.util.JapsaMath;

/**
 * Implementing a flexible profile ProbFSM. This machine would have each 3 
 * state corresponding to a base of the underlying sequence
 * @author minhduc
 *
 */
public class ProfilePFSM{
	MachineState startState;
	MachineState endState;	

	static double costCopyCopy = 1, costCopyMismatch = 2, costCopyDelete = 4, costCopyInsert = 4,
			costDeleteCopy = 1, costDeleteMismatch = 2, costDeleteDelete = 4,
			costInsertCopy = 1, costInsertMismatch = 2, costInsertInsert = 4;

	static{
		double  matProb = 0.95,
				insProb = 0.025,
				delProb = 0.025;

		costInsertInsert = costCopyInsert = -JapsaMath.log2(insProb);
		costDeleteDelete = costCopyDelete = -JapsaMath.log2(delProb);

		double matCost = -JapsaMath.log2(matProb);

		costDeleteCopy = costInsertCopy = costCopyCopy  = matCost - JapsaMath.log2(0.9);
		costDeleteMismatch = costInsertMismatch = costCopyMismatch = matCost - JapsaMath.log2(0.1) - JapsaMath.log2(1.0/3.0);

	}

	public ProfilePFSM(MachineState start, MachineState end){		
		this.startState = start;
		this.endState = end;
	}
	
	public Emission align(Sequence genSeq) {
		Emission retEmission = new Emission(endState, genSeq.length()-1, 0);
		retEmission.myCost = genSeq.length() * (costCopyInsert + 2);//plus 2 because of insert

		Emission currentEmission, lastEmission;
		currentEmission = lastEmission = new Emission(startState, -1, 0);
		currentEmission.myCost = 0;

		HashMap<String, Emission> hash = new HashMap<String, Emission>();
		int count = 0;
		while (currentEmission != null){
			if (currentEmission.gPos > 30 && currentEmission.myCost > 2.5 * currentEmission.gPos){
				Emission tmp = currentEmission.next;
				//helping GC to gabbabe collect current state			
				hash.remove(Emission.hashKey(currentEmission.toState, currentEmission.gPos, currentEmission.iteration));			
				currentEmission.next = null;
				currentEmission = tmp;

				continue;
			}

			count ++;
			if (count % 10000 == 0){				
				System.out.println(count + " " + hash.size() + "  " + 
						currentEmission.toState.getName() + " " + 
						currentEmission.iteration + " " +
						currentEmission.gPos + " " + 
						currentEmission.myCost + "###" +
						retEmission.myCost + " " + retEmission.iteration);

			}
			if (currentEmission.toState == endState){//If I am at the endState, check with the current best
				if (currentEmission.myCost < retEmission.myCost){
					//TODO: May need to remove the emission links of the existing ret to make thing easy for GC
					retEmission = currentEmission;					
				}
				Emission tmp = currentEmission.next;
				//helping GC to gabbabe collect current state			
				hash.remove(Emission.hashKey(currentEmission.toState, currentEmission.gPos, currentEmission.iteration));			
				currentEmission.next = null;
				currentEmission = tmp;

				continue;
			}
			//assert: the currentState is NOT the end
			String hashKey;
			Emission nextEmission;
			double cost;

			ArrayList<Transition> transitions = currentEmission.toState.getTransitions();			
			for (Transition tran: transitions){
				int nextPos = -2;
				byte nextBase = -1;
				int nextInt = currentEmission.iteration + tran.getIterIncrease();
				if(tran.getState() == endState){					
					//only consider this if the 
					if (currentEmission.gPos == genSeq.length() - 1){//done
						nextPos = currentEmission.gPos;
						//nextByte = -1;//dont care
					}
				}else if (tran instanceof Transition.InsertTransition){
					if (currentEmission.gPos < genSeq.length() - 1){
						nextPos = currentEmission.gPos + 1;
						nextBase = genSeq.getBase(nextPos);
					}					
				}else if (tran instanceof Transition.DeleteTransition){
					nextPos = currentEmission.gPos;
					//nextBase = genSeq.getBase(nextPos);					
				}else if (tran instanceof Transition.CopyTransition){
					if (currentEmission.gPos < genSeq.length() - 1){
						nextPos = currentEmission.gPos + 1;
						nextBase = genSeq.getBase(nextPos);
					}
				}else{
					//have to be free transation
					nextPos = currentEmission.gPos;
					//nextBase = genSeq.getBase(nextPos);					
				}
				if (nextPos >= -1){
					//only considered
					cost = currentEmission.myCost + tran.emissionCost(nextBase);

					if (cost < retEmission.myCost){ 
						hashKey= Emission.hashKey(tran.getState(), nextPos, nextInt);
						nextEmission = hash.get(hashKey);
						if (nextEmission == null){
							nextEmission = new Emission(tran.getState(), nextPos, nextInt);
							nextEmission.myCost = cost;
							hash.put(hashKey, nextEmission);
							lastEmission.next = nextEmission;
							lastEmission = nextEmission;
							nextEmission.bwdEmission = currentEmission;

							//nextEmission.countDel = currentEmission.countDel + 1;
							//nextEmission.countIns = currentEmission.countIns;
							//nextEmission.countMG =  currentEmission.countMG;
							//nextEmission.countMB =  currentEmission.countMB;

						}else{
							if (nextEmission.myCost > cost){
								nextEmission.myCost = cost;
								nextEmission.bwdEmission = currentEmission;

								//nextEmission.countDel = currentEmission.countDel + 1;
								//nextEmission.countIns = currentEmission.countIns;
								//nextEmission.countMG =  currentEmission.countMG;
								//nextEmission.countMB =  currentEmission.countMB;
							}
						}//else - if nextstate != null
					}//if cost					
				}//nextPos
			}//for

			Emission tmp = currentEmission.next;

			//helping GC to gabbabe collect current state			
			hash.remove(Emission.hashKey(currentEmission.toState, currentEmission.gPos, currentEmission.iteration));			
			currentEmission.next = null;
			currentEmission = tmp;
		}
		//System.out.printf("Consider %d states\n",complexity);

		//System.out.printf("Estimate: %2d %3d %3d %3d %3d %8.4f %8.4f\n",retState.iter, retState.countMG, retState.countMB, 
		//		retState.countIns, retState.countDel, retState.score,
		//		retState.countMG * (matCost + matchCost) + retState.countMB *(matCost + misMatchCost) + 
		//		retState.countIns * insCost + retState.countDel * delCost);	return retState;//bestScore;

		//timer.systemInfo();
		//Runtime.getRuntime().gc();
		//timer.systemInfo();

		//timer.mark("toc");
		return retEmission;
		/*******************************************************************/
	}
}