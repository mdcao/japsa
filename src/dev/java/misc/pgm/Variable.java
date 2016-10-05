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

import java.util.ArrayList;
import java.util.HashMap;

/**
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 * 
 */
public class Variable {
	private static int CURRENT_INDEX = 0;
	// the set of all variables, provide fast access to variables
	private static ArrayList<Variable> varList = new ArrayList<Variable>();
	private static HashMap<String, Variable> varMap = new HashMap<String, Variable>();
	/**
	 * The index of a variable. Each variable will have an unique index
	 */

	private int varIndex = 0;
	private String name;

	/**
	 * Creating a new variable, check if a variable of the same name exists
	 * varIndex is assigned in the order of variable created.
	 * 
	 * @param name
	 */
	public Variable(String name) {
		this.name = name;
		if (varMap.get(name) != null) {
			throw new RuntimeException("Variable named " + name + " exists");
		}

		varIndex = CURRENT_INDEX++;
		varList.add(this);
		varMap.put(name, this);
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	public String toString() {
		return name;
	}

	/**
	 * @param name
	 *            the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Return the index of this variable
	 * 
	 * @return
	 */
	public int getIndex() {
		// Note that there is no setIndex method
		return varIndex;
	}

	public static Variable getVariable(int index) {
		return varList.get(index);
	}

	public static Variable getVariable(String varName) {
		return varMap.get(varName);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
}
