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
 * 02/05/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package misc.pgm;

import java.util.Arrays;

/**
 * Representing a variable in a probabilistic graphical model.
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 * 
 */
public class DiscreteVariable extends Variable {

	/**
	 * Cardinality of the variable
	 */
	private int card;
	private String[] states;

	public DiscreteVariable(String name, int card) {
		super(name);

		this.card = card;
		states = new String[card];
		// default the state names as 0 1 2 ...
		for (int i = 0; i < card; i++)
			states[i] = i + "";
	}

	/**
	 * Create a discrete variable with a predefined state names
	 * 
	 * @param name
	 * @param statesName
	 */
	public DiscreteVariable(String name, String[] statesName) {
		super(name);
		this.card = statesName.length;
		states = Arrays.copyOf(statesName, card);
	}

	public void setState(int stateNumber, String state) {
		states[stateNumber] = state;
	}

	public String getState(int stateNumber) {
		return states[stateNumber];
	}

	public int cardinality() {
		return card;
	}

	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
	}
}
