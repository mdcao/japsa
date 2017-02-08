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

import japsa.bio.alignment.ppfsm.state.MachineState;
import japsa.bio.alignment.ppfsm.transition.Transition;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

/**
 * @author minhduc
 */
public class VNTRpThreeSM extends ProfilePFSM{

	
	/**
	 * 
	 * Create a profileFSM from three sequence: lflank, rflank, and period
	 *
	 * @param seq
	 * @param repStart
	 * @param repEnd
	 */
	public VNTRpThreeSM(Sequence seqLeft, Sequence seqRep, Sequence seqRight) {
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
		MachineState[] copyLeft   = new MachineState[seqLeft.length()];
		MachineState[] insertLeft = new MachineState[seqLeft.length()];
		MachineState[] deleteLeft = new MachineState[seqLeft.length()];

		//Now create states		
		MachineState insertFirstLeft = new MachineState("LfI");

		int i = 0;
		copyLeft[i] = new MachineState("LfC_" + i, seqLeft.getBase(i));		
		deleteLeft[i] = new MachineState("LfD_" + i);
		insertLeft[i] = new MachineState("LfI_" + i);

		Transition tran = new Transition.CopyTransition(copyLeft[0], costCopyCopy, costCopyMismatch);		
		startState.addTransition(tran);

		tran = new Transition.DeleteTransition(deleteLeft[0], costCopyDelete);
		startState.addTransition(tran);

		tran = new Transition.InsertTransition(insertFirstLeft, costCopyInsert);
		startState.addTransition(tran);

		tran = new Transition.CopyTransition(copyLeft[i], costInsertCopy, costInsertMismatch);
		insertFirstLeft.addTransition(tran);

		tran = new Transition.InsertTransition(insertFirstLeft, costInsertInsert);
		insertFirstLeft.addTransition(tran);

		for (i = 1; i < seqLeft.length(); i++) {
			//Create states
			copyLeft[i] = new MachineState("LfC_" + i, seqLeft.getBase(i));
			deleteLeft[i] = new MachineState("LfD_" + i);
			insertLeft[i] = new MachineState("LfI_" + i);

			//link the previous copy state
			tran = new Transition.CopyTransition(copyLeft[i], costCopyCopy, costCopyMismatch);
			copyLeft[i - 1].addTransition(tran);

			tran = new Transition.DeleteTransition(deleteLeft[i], costCopyDelete);
			copyLeft[i - 1].addTransition(tran);

			tran = new Transition.InsertTransition(insertLeft[i - 1], costCopyInsert);
			copyLeft[i - 1].addTransition(tran);

			//link the previous insert state
			tran = new Transition.CopyTransition(copyLeft[i], costInsertCopy, costInsertMismatch);
			insertLeft[i - 1].addTransition(tran);

			tran = new Transition.InsertTransition(insertLeft[i - 1], costInsertInsert);
			insertLeft[i - 1].addTransition(tran);

			//link the previous delete state
			tran = new Transition.CopyTransition(copyLeft[i], costDeleteCopy, costDeleteMismatch);
			deleteLeft[i - 1].addTransition(tran);

			tran = new Transition.DeleteTransition(deleteLeft[i], costDeleteDelete);
			deleteLeft[i - 1].addTransition(tran);
		}		

		i = seqLeft.length() - 1;

		//link the previous copy state
		tran = new Transition.InsertTransition(insertLeft[i], costCopyInsert);
		copyLeft[i].addTransition(tran);

		tran = new Transition.InsertTransition(insertLeft[i], costInsertInsert);
		insertLeft[i].addTransition(tran);

		tran = new Transition.FreeTransition(startRepeat, 0);
		tran.setIterIncrease(1);
		copyLeft[i].addTransition(tran);		

		tran = new Transition.FreeTransition(startRepeat, 0);
		tran.setIterIncrease(1);
		insertLeft[i].addTransition(tran);

		tran = new Transition.FreeTransition(startRepeat, 0);
		tran.setIterIncrease(1);
		deleteLeft[i].addTransition(tran);

		//do the repeat
		MachineState[] copyRep   = new MachineState[seqRep.length()];
		MachineState[] insertRep = new MachineState[seqRep.length()];
		MachineState[] deleteRep = new MachineState[seqRep.length()];

		//Now create states		
		MachineState insertFirstRep = new MachineState("ReI");

		i = 0;
		copyRep[i] = new MachineState("ReC_" + i, seqRep.getBase(i));		
		deleteRep[i] = new MachineState("ReD_" + i);
		insertRep[i] = new MachineState("ReI_" + i);

		tran = new Transition.CopyTransition(copyRep[0], costCopyCopy, costCopyMismatch);		
		startRepeat.addTransition(tran);

		tran = new Transition.DeleteTransition(deleteRep[0], costCopyDelete);
		startRepeat.addTransition(tran);

		tran = new Transition.InsertTransition(insertFirstRep, costCopyInsert);
		startRepeat.addTransition(tran);

		tran = new Transition.CopyTransition(copyRep[i], costInsertCopy, costInsertMismatch);
		insertFirstRep.addTransition(tran);

		tran = new Transition.InsertTransition(insertFirstRep, costInsertInsert);
		insertFirstRep.addTransition(tran);

		for (i = 1; i < seqRep.length(); i++) {
			//Create states
			copyRep[i] = new MachineState("ReC_" + i, seqRep.getBase(i));
			deleteRep[i] = new MachineState("ReD_" + i);
			insertRep[i] = new MachineState("ReI_" + i);

			//link the previous copy state
			tran = new Transition.CopyTransition(copyRep[i], costCopyCopy, costCopyMismatch);
			copyRep[i - 1].addTransition(tran);

			tran = new Transition.DeleteTransition(deleteRep[i], costCopyDelete);
			copyRep[i - 1].addTransition(tran);

			tran = new Transition.InsertTransition(insertRep[i - 1], costCopyInsert);
			copyRep[i - 1].addTransition(tran);

			//link the previous insert state
			tran = new Transition.CopyTransition(copyRep[i], costInsertCopy, costInsertMismatch);
			insertRep[i - 1].addTransition(tran);

			tran = new Transition.InsertTransition(insertRep[i - 1], costInsertInsert);
			insertRep[i - 1].addTransition(tran);

			//link the previous delete state
			tran = new Transition.CopyTransition(copyRep[i], costDeleteCopy, costDeleteMismatch);
			deleteRep[i - 1].addTransition(tran);

			tran = new Transition.DeleteTransition(deleteRep[i], costDeleteDelete);
			deleteRep[i - 1].addTransition(tran);
		}		

		i = seqRep.length() - 1;

		//link the previous copy state
		tran = new Transition.InsertTransition(insertRep[i], costCopyInsert);
		copyRep[i].addTransition(tran);

		tran = new Transition.InsertTransition(insertRep[i], costInsertInsert);
		insertRep[i].addTransition(tran);

		tran = new Transition.FreeTransition(endRepeat, 0);			
		copyRep[i].addTransition(tran);		

		tran = new Transition.FreeTransition(endRepeat, 0);			
		insertRep[i].addTransition(tran);

		tran = new Transition.FreeTransition(endRepeat, 0);			
		deleteRep[i].addTransition(tran);

		//do the right
		MachineState[] copyRight   = new MachineState[seqRight.length()];
		MachineState[] insertRight = new MachineState[seqRight.length()];
		MachineState[] deleteRight = new MachineState[seqRight.length()];

		//Now create states		
		MachineState insertFirstRight = new MachineState("RfI");

		i = 0;
		copyRight[i] = new MachineState("RfC_" + i, seqRight.getBase(i));		
		deleteRight[i] = new MachineState("RfD_" + i);
		insertRight[i] = new MachineState("RfI_" + i);

		tran = new Transition.CopyTransition(copyRight[0], costCopyCopy, costCopyMismatch);		
		endRepeat.addTransition(tran);

		tran = new Transition.DeleteTransition(deleteRight[0], costCopyDelete);
		endRepeat.addTransition(tran);

		tran = new Transition.InsertTransition(insertFirstRight, costCopyInsert);
		endRepeat.addTransition(tran);

		tran = new Transition.CopyTransition(copyRight[i], costInsertCopy, costInsertMismatch);
		insertFirstRight.addTransition(tran);

		tran = new Transition.InsertTransition(insertFirstRight, costInsertInsert);
		insertFirstRight.addTransition(tran);

		for (i = 1; i < seqRight.length(); i++) {
			//Create states
			copyRight[i] = new MachineState("RfC_" + i, seqRight.getBase(i));
			deleteRight[i] = new MachineState("RfD_" + i);
			insertRight[i] = new MachineState("RfI_" + i);

			//link the previous copy state
			tran = new Transition.CopyTransition(copyRight[i], costCopyCopy, costCopyMismatch);
			copyRight[i - 1].addTransition(tran);

			tran = new Transition.DeleteTransition(deleteRight[i], costCopyDelete);
			copyRight[i - 1].addTransition(tran);

			tran = new Transition.InsertTransition(insertRight[i - 1], costCopyInsert);
			copyRight[i - 1].addTransition(tran);

			//link the previous insert state
			tran = new Transition.CopyTransition(copyRight[i], costInsertCopy, costInsertMismatch);
			insertRight[i - 1].addTransition(tran);

			tran = new Transition.InsertTransition(insertRight[i - 1], costInsertInsert);
			insertRight[i - 1].addTransition(tran);

			//link the previous delete state
			tran = new Transition.CopyTransition(copyRight[i], costDeleteCopy, costDeleteMismatch);
			deleteRight[i - 1].addTransition(tran);

			tran = new Transition.DeleteTransition(deleteRight[i], costDeleteDelete);
			deleteRight[i - 1].addTransition(tran);
		}		

		i = seqRight.length() - 1;

		//link the previous copy state
		tran = new Transition.InsertTransition(insertRight[i], costCopyInsert);
		copyRight[i].addTransition(tran);

		tran = new Transition.InsertTransition(insertRight[i], costInsertInsert);
		insertRight[i].addTransition(tran);

		tran = new Transition.FreeTransition(endState, 0);			
		copyRight[i].addTransition(tran);		

		tran = new Transition.FreeTransition(endState, 0);			
		insertRight[i].addTransition(tran);

		tran = new Transition.FreeTransition(endState, 0);			
		deleteRight[i].addTransition(tran);

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

		VNTRpThreeSM pFMS = new VNTRpThreeSM(leftSeq, repSeq, rightSeq);
		Sequence seq;
		while ((seq = reader.nextSequence(dnaAlphabet)) != null){
			pFMS.align(seq);
		}
	}
}
