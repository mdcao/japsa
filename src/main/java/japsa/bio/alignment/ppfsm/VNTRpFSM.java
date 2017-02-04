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

/*                           Revision History                                
 * 31/01/2017 - Minh Duc Cao: Created                                        
 ****************************************************************************/
package japsa.bio.alignment.ppfsm;

import java.util.ArrayList;

import japsa.bio.alignment.ppfsm.state.MachineState;
import japsa.bio.alignment.ppfsm.transition.Transition;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;

/**
 * @author minhduc
 */
public class VNTRpFSM {
	MachineState.StartState startState;
	MachineState.EndState endState;

	MachineState.CopyState[] copyState;
	MachineState.DeleteState[] deleteState;
	MachineState.InsertState[] insertState;

	static double costCopyCopy = 1, costCopyMismatch = 2, costCopyDelete = 4, costCopyInsert = 4,
			costDeleteCopy = 1, costDeleteMismatch = 2, costDeleteDelete = 4,
			costInsertCopy = 1, costInsertMismatch = 2, costInsertInsert = 4;
	Sequence mSeq;

	ArrayList<MachineState> states = new ArrayList<MachineState>();

	/**
	 * +
	 * Create a profileFSM from three sequence: lflank, rflank, and period
	 *
	 * @param seq
	 * @param repStart
	 * @param repEnd
	 */
	public VNTRpFSM(Sequence lflank, Sequence rSeq, Sequence rflank) {
		SequenceBuilder concat = new SequenceBuilder(lflank, lflank.length() + rSeq.length() + rflank.length());

		concat.append(rSeq);
		concat.append(rflank);

		mSeq = concat.toSequence();

		//start, consume none
		startState = new MachineState.StartState();
		//end have gone through all the hidden state
		endState = new MachineState.EndState();

		states.add(startState);
		states.add(endState);

		MachineState.InsertState insertFirst = new MachineState.InsertState("IS");
		states.add(insertFirst);

		copyState = new MachineState.CopyState[mSeq.length()];
		insertState = new MachineState.InsertState[mSeq.length()];
		deleteState = new MachineState.DeleteState[mSeq.length()];

		int i = 0;
		copyState[i] = new MachineState.CopyState("C" + i, mSeq.getBase(i));
		deleteState[i] = new MachineState.DeleteState("D" + i);
		insertState[i] = new MachineState.InsertState("I" + i);

		states.add(copyState[i]);
		states.add(insertState[i]);
		states.add(deleteState[i]);


		Transition.CopyTransition cTransition =
				new Transition.CopyTransition(copyState[0], costCopyCopy, costCopyMismatch);
		startState.setCopyTransition(cTransition);

		Transition.DeleteTransition dTransition =
				new Transition.DeleteTransition(deleteState[0], costCopyDelete);
		startState.setDeleteTransition(dTransition);

		Transition.InsertTransition iTransition =
				new Transition.InsertTransition(insertFirst, costCopyInsert);
		startState.setInsertTransition(iTransition);


		for (i = 1; i < mSeq.length(); i++) {
			//Create states
			copyState[i] = new MachineState.CopyState("C" + i, mSeq.getBase(i));
			deleteState[i] = new MachineState.DeleteState("D" + i);
			insertState[i] = new MachineState.InsertState("I" + i);

			//link the previous copy state
			cTransition = new Transition.CopyTransition(copyState[i], costCopyCopy, costCopyMismatch);
			copyState[i - 1].setCopyTransition(cTransition);

			dTransition = new Transition.DeleteTransition(deleteState[i], costCopyDelete);
			copyState[i - 1].setDeleteTransition(dTransition);

			iTransition = new Transition.InsertTransition(insertState[i - 1], costCopyInsert);
			copyState[i - 1].setInsertTransition(iTransition);

			//link the previous insert state
			cTransition = new Transition.CopyTransition(copyState[i], costInsertCopy, costInsertMismatch);
			insertState[i - 1].setCopyTransition(cTransition);

			iTransition = new Transition.InsertTransition(insertState[i - 1], costInsertInsert);
			insertState[i - 1].setInsertTransition(iTransition);

			//link the previous delete state
			cTransition = new Transition.CopyTransition(copyState[i], costDeleteCopy, costDeleteMismatch);
			deleteState[i - 1].setCopyTransition(cTransition);

			dTransition = new Transition.DeleteTransition(deleteState[i], costDeleteDelete);
			deleteState[i - 1].setDeleteTransition(dTransition);

		}
		i = mSeq.length() - 1;

		//link the previous copy state
		//copyState[i].copyTransition = null;
		//copyState[i].deleteTransition = null;
		iTransition = new Transition.InsertTransition(insertState[i], costCopyInsert);
		copyState[i].setInsertTransition(iTransition);

		//link the previous insert state
		//insertState[i].copyTransition   = null;

		iTransition = new Transition.InsertTransition(insertState[i], costInsertInsert);
		insertState[i].setInsertTransition(iTransition);
		//insertState[i-1].deleteTransition = null;

		//link the previous delete state
		//deleteState[i].copyTransition   = null;
		//deleteState[i].deleteTransition = null;
		//deleteState[i-1].insertTransition = null;

		//free links to end state
		copyState[i].setEndTransition(new Transition.EndTransition(endState));
		insertState[i].setEndTransition(new Transition.EndTransition(endState));
		deleteState[i].setEndTransition(new Transition.EndTransition(endState));
	}


	public Emission align(Sequence genSeq) {
		//First find a possible path, not nesseary the best. ei the null hypothesis
		//basically go through all the copy states until either the profile
		//sequence (mSeq) or the generated sequence (genSeq) is exhausted. The go
		//to through the deletion/insertion to finish off
		MachineState currentState = startState;
		double bestCost = 0;
		int currentIndexGen = 0;
		while (currentState != endState) {
			//first check if the whole genSeq has been generated

			if (currentIndexGen < genSeq.length()) {
				//do copy as much as I can
				Transition.CopyTransition cTransition = currentState.getCopyTransition();
				if (cTransition != null) {
					bestCost += cTransition.emissionCost(genSeq.getBase(currentIndexGen));
					currentIndexGen++;
					currentState = cTransition.getCopyState();
					//then the next iteration
				} else {
					//at this point, have gone through the profile
					while (currentIndexGen < genSeq.length()) {
						Transition.InsertTransition iTransition = currentState.getInsertTransition();
						bestCost += iTransition.emissionCost(genSeq.getBase(currentIndexGen));
						currentIndexGen++;
						currentState = iTransition.getInsertState();
					}//assert gone through the gene
					Transition.EndTransition eTransition = currentState.getEndTransition();
					if (eTransition == null) {
						throw new RuntimeException("Something very very wrong, shouldn't get here 1");
					}
					bestCost += eTransition.getCost();
					currentState = eTransition.getEndState();//which is also the endState and hence can terminate loop
				}
			} else {
				//need to move to deletion the rest of the profile
				Transition.DeleteTransition dTransition = null;
				while ((dTransition = currentState.getDeleteTransition()) != null) {
					bestCost += dTransition.getCost();
					currentState = dTransition.getDeleteState();
				}
				//assert" dTransition == null meaning the currentState must go to End
				Transition.EndTransition eTransition = currentState.getEndTransition();
				if (eTransition == null) {
					throw new RuntimeException("Something very very wrong, shouldn't get here 2");
				}
				bestCost += eTransition.getCost();
				currentState = eTransition.getEndState();//which is also the endState and hence can terminate loop
			}
		}


		//Transition
		//MachineState nextCopy = startState.copyTransition.toState;
		/*******************************************************************
         Emission retEmission = new Emission(states[0], mSeq.length() - 1, genSeq.length() - 1);
         retEmission.myCost = genSeq.length() * (insEmissionCost + 4);


         Emission currentEmission, finalEmission;//current pointer and last pointer on the linked-list
         currentEmission = finalEmission = new Emission(states[0], -1, -1);
         currentEmission.myCost = 0;

         HashMap<String, Emission> hash = new HashMap<String, Emission>();

         while (currentEmission != null) {//linked list not exhausted
         if (currentEmission.gPos >= genSeq.length() - 1) {
         //done generating genSeq
         if (currentEmission.myCost < retEmission.myCost) {
         retEmission = currentEmission;
         }
         } else if (currentEmission.myCost < retEmission.myCost) {
         String hashKey;
         Emission nextEmission;
         double cost;

         //1. consider deletion if profile has something to offer
         if (currentEmission.mPos + 1 < mSeq.length() && currentEmission.toState.delState != null) {
         cost = currentEmission.myCost + currentEmission.toState.delCost;

         hashKey = Emission.hashKey(currentEmission.toState.delState.name, currentEmission.mPos + 1, currentEmission.gPos);
         nextEmission = hash.get(hashKey);
         if (nextEmission == null) {
         nextEmission = new Emission(currentEmission.toState.delState, currentEmission.mPos + 1, currentEmission.gPos);
         nextEmission.type = EmissionType.DELETION;

         nextEmission.myCost = cost;
         hash.put(hashKey, nextEmission);

         finalEmission.next = nextEmission;
         finalEmission = nextEmission;
         nextEmission.bwdEmission = currentEmission;
         } else {
         if (nextEmission.myCost > cost) {
         nextEmission.myCost = cost;
         nextEmission.bwdEmission = currentEmission;
         nextEmission.type = EmissionType.DELETION;
         }
         }//else - if nextstate != null
         }//if mPos

         //2. insertion
         if (currentEmission.gPos + 1 < genSeq.length() && currentEmission.toState.insState != null) {
         cost = currentEmission.myCost + currentEmission.toState.insCost + this.insEmissionCost;

         hashKey = Emission.hashKey(currentEmission.toState.insState.name, currentEmission.mPos, currentEmission.gPos + 1);
         nextEmission = hash.get(hashKey);
         if (nextEmission == null) {
         nextEmission = new Emission(currentEmission.toState.insState, currentEmission.mPos, currentEmission.gPos + 1);
         nextEmission.myCost = cost;
         nextEmission.type = EmissionType.INSERTION;
         hash.put(hashKey, nextEmission);

         finalEmission.next = nextEmission;
         finalEmission = nextEmission;
         nextEmission.bwdEmission = currentEmission;

         } else {
         if (nextEmission.myCost > cost) {
         nextEmission.myCost = cost;
         nextEmission.bwdEmission = currentEmission;
         nextEmission.type = EmissionType.INSERTION;
         }
         }
         }

         //3.Match
         if (currentEmission.gPos + 1 < genSeq.length() && currentEmission.mPos + 1 < mSeq.length()) {
         EmissionType type = EmissionType.COPY;

         if (mSeq.getBase(currentEmission.mPos + 1) == genSeq.getBase(currentEmission.gPos + 1)) {
         cost = currentEmission.myCost + currentEmission.toState.matchCost + currentEmission.toState.copyCost;
         } else {
         cost = currentEmission.myCost + currentEmission.toState.matchCost + currentEmission.toState.changeCost + this.changeEmissionCost;
         type = EmissionType.MUTATE;
         }

         if (cost < retEmission.myCost) {
         hashKey = Emission.hashKey(currentEmission.toState.matchState.name, currentEmission.mPos + 1, currentEmission.gPos + 1);
         nextEmission = hash.get(hashKey);
         if (nextEmission == null) {
         nextEmission = new Emission(currentEmission.toState.matchState, currentEmission.mPos + 1, currentEmission.gPos + 1);
         nextEmission.type = type;

         nextEmission.myCost = cost;
         hash.put(hashKey, nextEmission);
         finalEmission.next = nextEmission;
         finalEmission = nextEmission;
         nextEmission.bwdEmission = currentEmission;
         } else {
         if (nextEmission.myCost > cost) {
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

         hash.remove(Emission.hashKey(currentEmission.toState.name, currentEmission.mPos, currentEmission.gPos));
         currentEmission.next = null;

         currentEmission = tmp;
         }

         Logging.info("Hash = " + hash.size());

         //timer.systemInfo();
         //Runtime.getRuntime().gc();
         //timer.systemInfo();

         //timer.mark("toc");

         return retEmission;
         /*******************************************************************/
		return null;


	}

}
