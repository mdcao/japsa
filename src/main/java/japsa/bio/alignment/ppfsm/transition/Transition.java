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
package japsa.bio.alignment.ppfsm.transition;

import japsa.bio.alignment.ppfsm.state.MachineState;

/**
 * Transition between states
 * @author minhduc
 *
 */
public abstract class Transition{

	/**
	 * The probability of this transition
	 */
	double cost = 0;

	/**
	 * The count of this transition
	 */
	//int count = 1;	
	int iterIncrease = 0;
	protected final int numGenerate;
	MachineState toState;

	public Transition(MachineState state, int numGenerate){
		toState = state;
		this.numGenerate = numGenerate;
	}

	/**
	 * @return the cost
	 */
	public double getCost() {
		return cost;
	}

	public abstract double emissionCost(byte base);

	/**
	 * @return the iteration
	 */
	public int getIterIncrease() {
		return iterIncrease;
	}

	/**
	 * @param iteration the iteration to set
	 */
	public void setIterIncrease(int it) {
		this.iterIncrease = it;
	}

	/**
	 * @param cost the cost to set
	 */
	public void setCost(double cost) {
		this.cost = cost;
	}

	public MachineState getState() {
		return toState;
	}


	public static class CopyTransition extends Transition{		
		double copyCost;
		double mismatchCost;
		//MachineState.CopyState toState;

		public CopyTransition(MachineState state, double copyCost, double mismatchCost) {
			super(state,1);			
			this.copyCost = copyCost;
			this.mismatchCost = mismatchCost;
		}

		public double emissionCost(byte base) {			
			if (base == toState.getBase())
				return copyCost;
			else
				return mismatchCost;
		}
	}

	public static class DeleteTransition extends Transition{		
		double deleteCost;		

		public DeleteTransition(MachineState state, double cost) {
			super(state,0);
			this.deleteCost = cost;
		}

		public double emissionCost(byte base) {			
			return deleteCost;
		}		
	}

	public static class InsertTransition extends Transition{		
		double insertCost;

		public InsertTransition(MachineState state, double cost) {
			super(state,1);			
			this.insertCost = cost;
		}

		public double emissionCost(byte base) {			
			return insertCost + 2;
		}
				
	}

	public static class FreeTransition extends Transition{
		
		public FreeTransition(MachineState state, double cost) {
			super(state,0);
			this.cost = cost;
		}

		public double emissionCost(byte base) {			
			return cost;
		}
	}
}