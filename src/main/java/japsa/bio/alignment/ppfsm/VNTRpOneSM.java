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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import japsa.bio.alignment.ppfsm.state.MachineState;
import japsa.bio.alignment.ppfsm.transition.Transition;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.JapsaMath;

/**
 * @author minhduc
 */
public class VNTRpOneSM extends ProfilePFSM{
	/**
	 * 
	 * Create a profileFSM from three sequence: lflank, rflank, and period. This is one state machine
	 * (that is copy/insert/delete merged into one for each profile position)
	 *
	 * @param seq
	 * @param repStart
	 * @param repEnd
	 */
	public VNTRpOneSM(Sequence seqLeft, Sequence seqRep, Sequence seqRight) {
		super(new MachineState("S"), new MachineState("E"));
		
		System.out.println("costCopyCopy = " + costCopyCopy);
		System.out.println("costCopyMismatch = " + costCopyMismatch);
		System.out.println("costCopyDelete = " + costCopyDelete);
		System.out.println("costCopyInsert = " + costCopyInsert);		

		System.out.println("costDeleteCopy = " + costDeleteCopy);
		System.out.println("costDeleteMismatch = " + costDeleteMismatch);
		System.out.println("costDeleteDelete = " + costDeleteDelete);		

		System.out.println("costInsertCopy = " + costInsertCopy);
		System.out.println("costInsertMismatch = " + costInsertMismatch);
		System.out.println("costInsertInsert = " + costInsertInsert);		

		
		MachineState startRepeat =  new MachineState("SR");//which is also the end of left
		MachineState endRepeat =  new MachineState("ER");

		//Do the left flank
		MachineState[] nodeLeft   = new MachineState[seqLeft.length()];

		//Now create states
		int i = 0;
		nodeLeft[i] = new MachineState("Lf_" + i, seqLeft.getBase(i));

		Transition tran = new Transition.CopyTransition(nodeLeft[0], costCopyCopy, costCopyMismatch);		
		startState.addTransition(tran);

		tran = new Transition.DeleteTransition(nodeLeft[0], costCopyDelete);
		startState.addTransition(tran);		

		tran = new Transition.InsertTransition(startState, costCopyInsert);
		startState.addTransition(tran);		

		for (i = 1; i < seqLeft.length(); i++) {
			//Create states
			nodeLeft[i] = new MachineState("LfC_" + i, seqLeft.getBase(i));

			//link the previous copy state
			tran = new Transition.CopyTransition(nodeLeft[i], costCopyCopy, costCopyMismatch);
			nodeLeft[i - 1].addTransition(tran);

			tran = new Transition.DeleteTransition(nodeLeft[i], costCopyDelete);
			nodeLeft[i - 1].addTransition(tran);

			tran = new Transition.InsertTransition(nodeLeft[i - 1], costCopyInsert);
			nodeLeft[i - 1].addTransition(tran);
		}		

		i = seqLeft.length() - 1;

		//link the previous copy state
		tran = new Transition.InsertTransition(nodeLeft[i], costCopyInsert);
		nodeLeft[i].addTransition(tran);
		
		tran = new Transition.FreeTransition(startRepeat, 0);
		tran.setIterIncrease(1);
		nodeLeft[i].addTransition(tran);		
				
		//do the repeat
		MachineState[] nodeRep   = new MachineState[seqRep.length()];
		
		i = 0;
		nodeRep[i] = new MachineState("ReC_" + i, seqRep.getBase(i));

		tran = new Transition.CopyTransition(nodeRep[0], costCopyCopy, costCopyMismatch);		
		startRepeat.addTransition(tran);

		tran = new Transition.DeleteTransition(nodeRep[0], costCopyDelete);
		startRepeat.addTransition(tran);

		tran = new Transition.InsertTransition(startRepeat, costCopyInsert);
		startRepeat.addTransition(tran);
		
		for (i = 1; i < seqRep.length(); i++) {
			//Create states
			nodeRep[i] = new MachineState("ReC_" + i, seqRep.getBase(i));

			//link the previous copy state
			tran = new Transition.CopyTransition(nodeRep[i], costCopyCopy, costCopyMismatch);
			nodeRep[i - 1].addTransition(tran);

			tran = new Transition.DeleteTransition(nodeRep[i], costCopyDelete);
			nodeRep[i - 1].addTransition(tran);

			tran = new Transition.InsertTransition(nodeRep[i - 1], costCopyInsert);
			nodeRep[i - 1].addTransition(tran);
		}		

		i = seqRep.length() - 1;

		//link the previous copy state
		tran = new Transition.InsertTransition(nodeRep[i], costCopyInsert);
		nodeRep[i].addTransition(tran);
		
		tran = new Transition.FreeTransition(endRepeat, 0);			
		nodeRep[i].addTransition(tran);
		
		//do the right
		MachineState[] nodeRight   = new MachineState[seqRight.length()];
		
		i = 0;
		nodeRight[i] = new MachineState("RfC_" + i, seqRight.getBase(i));		
		
		tran = new Transition.CopyTransition(nodeRight[i], costCopyCopy, costCopyMismatch);		
		endRepeat.addTransition(tran);

		tran = new Transition.DeleteTransition(nodeRight[i], costCopyDelete);
		endRepeat.addTransition(tran);

		tran = new Transition.InsertTransition(endRepeat, costCopyInsert);
		endRepeat.addTransition(tran);

		for (i = 1; i < seqRight.length(); i++) {
			//Create states
			nodeRight[i] = new MachineState("RfC_" + i, seqRight.getBase(i));

			//link the previous copy state
			tran = new Transition.CopyTransition(nodeRight[i], costCopyCopy, costCopyMismatch);
			nodeRight[i - 1].addTransition(tran);

			tran = new Transition.DeleteTransition(nodeRight[i], costCopyDelete);
			nodeRight[i - 1].addTransition(tran);

			tran = new Transition.InsertTransition(nodeRight[i - 1], costCopyInsert);
			nodeRight[i - 1].addTransition(tran);			
		}		

		i = seqRight.length() - 1;

		//link the previous copy state
		tran = new Transition.InsertTransition(nodeRight[i], costCopyInsert);
		nodeRight[i].addTransition(tran);
		
		tran = new Transition.FreeTransition(endState, 0);			
		nodeRight[i].addTransition(tran);
		
		////////////////////////////////////////////////////////////////////
		tran = new Transition.FreeTransition(startRepeat, 0);
		tran.setIterIncrease(1);
		endRepeat.addTransition(tran);
	}
	

	public static void main(String [] args) throws IOException{
		SequenceReader reader = SequenceReader.getReader(args[0]);
		Alphabet dnaAlphabet = Alphabet.DNA16();
		Sequence leftSeq = reader.nextSequence(dnaAlphabet);
		Sequence repSeq = reader.nextSequence(dnaAlphabet);
		Sequence rightSeq = reader.nextSequence(dnaAlphabet);

		VNTRpOneSM pFMS = new VNTRpOneSM(leftSeq, repSeq, rightSeq);
		Sequence seq;
		while ((seq = reader.nextSequence(dnaAlphabet)) != null){
			pFMS.align(seq);
		}
	}
}
