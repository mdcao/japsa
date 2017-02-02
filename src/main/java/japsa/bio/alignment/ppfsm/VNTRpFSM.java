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

import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;

/**
 * @author minhduc
 *
 */
public class VNTRpFSM{
	MachineState start,end;
	MachineState[] copyState, deleteState, insertState;
	double costCopyCopy = 1, costCopyMismatch = 2, costCopyDelete = 4, costCopyInsert = 4, 
			costDeleteCopy = 1, costDeleteMismatch = 2, costDeleteDelete = 4, 
			costInsertCopy = 1, costInsertMismatch = 2, costInsertInsert = 4;
	
	Sequence mSeq;
	/**
	 * Create a profileFSM from three sequence: lflank, rflank, and period
	 * @param seq
	 * @param repStart
	 * @param repEnd
	 */	
	public VNTRpFSM(Sequence lflank, Sequence rSeq, Sequence rflank){
		SequenceBuilder concat = new SequenceBuilder(lflank, lflank.length() + rSeq.length() + rflank.length());

		concat.append(rSeq);		
		concat.append(rflank);

		mSeq = concat.toSequence();

		//start, consume none
		MachineState start = new MachineState("S");
		//end have gone through all the hidden state
		MachineState end = new MachineState("E");		
		MachineState insertFirst = new MachineState("IS");
		
		copyState = new MachineState[mSeq.length()];
		insertState = new MachineState[mSeq.length()];
		deleteState = new MachineState[mSeq.length()];

		int i = 0;		
		copyState  [i] = new MachineState("C" + i);
		deleteState[i] = new MachineState("D" + i);
		insertState[i] = new MachineState("I" + i);

		start.copyTransition = new Transition(copyState[0], costCopyCopy);
		start.mismatchTransition = new Transition(copyState[0], costCopyMismatch);
		start.deleteTransition = new Transition(deleteState[0], costCopyDelete);
		start.insertTransition = new Transition(insertFirst, costCopyInsert);
				
		for (i=1; i < mSeq.length();i++){
			//Create states
			copyState  [i] = new MachineState("C" + i);
			deleteState[i] = new MachineState("D" + i);
			insertState[i] = new MachineState("I" + i);
			
			//link the previous copy state
			copyState[i-1].copyTransition = new Transition(copyState[i], costCopyCopy);
			copyState[i-1].mismatchTransition = new Transition(copyState[i], costCopyMismatch);
			copyState[i-1].deleteTransition = new Transition(deleteState[i], costCopyDelete);
			copyState[i-1].insertTransition = new Transition(insertState[i-1], costCopyInsert);
			
			//link the previous insert state			
			insertState[i-1].copyTransition   = new Transition(copyState[i], costInsertCopy);
			insertState[i-1].mismatchTransition   = new Transition(copyState[i], costInsertMismatch);
			insertState[i-1].insertTransition = new Transition(insertState[i-1], costInsertInsert);
			//insertState[i-1].deleteTransition = null;			
			
			//link the previous delete state			
			deleteState[i-1].copyTransition   = new Transition(copyState[i], costDeleteCopy);
			deleteState[i-1].mismatchTransition   = new Transition(copyState[i], costDeleteMismatch);
			deleteState[i-1].deleteTransition = new Transition(deleteState[i], costDeleteDelete);
			//deleteState[i-1].insertTransition = null;		
								
		}
		i = mSeq.length() - 1;
		
		//link the previous copy state
		//copyState[i].copyTransition = null;//
		//copyState[i].mismatchTransition = null;
		//copyState[i].deleteTransition = null;
		copyState[i].insertTransition = new Transition(insertState[i], costCopyInsert);
		
		//link the previous insert state			
		//insertState[i].copyTransition   = null;
		//insertState[i-1].mismatchTransition  = null;
		insertState[i].insertTransition = new Transition(insertState[i], costInsertInsert);
		//insertState[i-1].deleteTransition = null;			
		
		//link the previous delete state			
		//deleteState[i].copyTransition   = null;
		//deleteState[i].mismatchTransition = null;
		//deleteState[i].deleteTransition = null;
		//deleteState[i-1].insertTransition = null;
		
		//free links to end state
		copyState[i].addTransition(end, 0);		
		insertState[i].addTransition(end, 0);
		deleteState[i].addTransition(end, 0);
		
	}
	
	
	

}
