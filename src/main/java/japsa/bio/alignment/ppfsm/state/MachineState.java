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

package japsa.bio.alignment.ppfsm.state;

import java.util.ArrayList;
import japsa.bio.alignment.ppfsm.transition.Transition;


/**
 * Implementation of a machine state for the profile FSM
 * @author minhduc
 *
 */

public class MachineState {
	//protected Transition.CopyTransition   copyTransition = null;
	//protected Transition.DeleteTransition deleteTransition = null; 
	//protected Transition.InsertTransition insertTransition = null;
	//protected Transition.EndTransition    endTransition = null;

	ArrayList<Transition> transitions = new ArrayList<Transition>();	
	private String name;
	private byte base;

	public MachineState(String name) {
		this.setName(name);
	}
	
	public MachineState(String name, byte base) {
		this(name);
		setBase(base);
	}
	/**
	 * Add a non-standard transition
	 * @param toState
	 * @param cost
	 */
	public void addTransition(Transition tran){		
		transitions.add(tran);
	}

	/**
	 * Get non-standard transition
	 * @return
	 */
	public ArrayList<Transition> getTransitions(){
		return transitions;
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}
	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * @return the name
	 */
	public byte getBase() {
		return base;
	}
	/**
	 * @param name the name to set
	 */
	public void setBase(byte base) {
		this.base = base;
	}



	//	/**
	//	 * Copy transition
	//	 * @author minhduc
	//	 *
	//	 */
	//	public static class CopyState extends MachineState{
	//		byte base;
	//		
	//		public CopyState(String name, byte base) {
	//			super(name);
	//			this.base = base;
	//		}		
	//		/**
	//		 * @return the base
	//		 */
	//		public byte getBase() {
	//			return base;
	//		}
	//		/**
	//		 * @param base the base to set
	//		 */
	//		public void setBase(byte base) {
	//			this.base = base;
	//		}
	//		/**
	//		 * @param copyTransition the copyTransition to set
	//		 */
	//		public void setCopyTransition(Transition.CopyTransition copyTransition) {
	//			this.copyTransition = copyTransition;
	//		}
	//		/**
	//		 * @param deleteTransition the deleteTransition to set
	//		 */
	//		public void setDeleteTransition(Transition.DeleteTransition deleteTransition) {
	//			this.deleteTransition = deleteTransition;
	//		}
	//		/**
	//		 * @param insertTransition the insertTransition to set
	//		 */
	//		public void setInsertTransition(Transition.InsertTransition insertTransition) {
	//			this.insertTransition = insertTransition;
	//		}
	//	}
	//
	//	/**
	//	 * Insert transition
	//	 * @author minhduc
	//	 *
	//	 */
	//	public static class InsertState extends MachineState{		
	//		public InsertState(String name) {
	//			super(name);			
	//		}		
	//		/**
	//		 * @param copyTransition the copyTransition to set
	//		 */
	//		public void setCopyTransition(Transition.CopyTransition copyTransition) {
	//			this.copyTransition = copyTransition;
	//		}		
	//		/**
	//		 * @param insertTransition the insertTransition to set
	//		 */
	//		public void setInsertTransition(Transition.InsertTransition insertTransition) {
	//			this.insertTransition = insertTransition;
	//		}
	//	}
	//	
	//	/**
	//	 * Copy transition
	//	 * @author minhduc
	//	 *
	//	 */
	//	public static class DeleteState extends MachineState{
	//		public DeleteState(String name) {
	//			super(name);		
	//		}		
	//		
	//		/**
	//		 * @param copyTransition the copyTransition to set
	//		 */
	//		public void setCopyTransition(Transition.CopyTransition copyTransition) {
	//			this.copyTransition = copyTransition;
	//		}
	//		/**
	//		 * @param deleteTransition the deleteTransition to set
	//		 */
	//		public void setDeleteTransition(Transition.DeleteTransition deleteTransition) {
	//			this.deleteTransition = deleteTransition;
	//		}		
	//	}
	//	
	//	public static class StartState extends MachineState{
	//		public StartState() {
	//			super("S");		
	//		}		
	//		/**
	//		 * @param copyTransition the copyTransition to set
	//		 */
	//		public void setCopyTransition(Transition.CopyTransition copyTransition) {
	//			this.copyTransition = copyTransition;
	//		}
	//		/**
	//		 * @param deleteTransition the deleteTransition to set
	//		 */
	//		public void setDeleteTransition(Transition.DeleteTransition deleteTransition) {
	//			this.deleteTransition = deleteTransition;
	//		}
	//
	//		/**
	//		 * @param insertTransition the insertTransition to set
	//		 */
	//		public void setInsertTransition(Transition.InsertTransition insertTransition) {
	//			this.insertTransition = insertTransition;
	//		}
	//	}
	//	//EndState: has no transition
	//	public static class EndState extends MachineState{
	//		public EndState() {
	//			super("E");		
	//		}
	//	}


}