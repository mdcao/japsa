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

package japsadev.misc.pgm;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.BitSet;
import java.util.Iterator;

/**
 * Implement a factor to be used in representing discrete probability
 * distributions in probabilistic graphical models.
 * 
 * A factor represents a distribution, not necessarily normalised. It can be
 * used to represent a joint distribution or a conditional distribution.
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 * 
 */
public class Factor implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 2710221433807184868L;

	private FactorScope scope;
	private double[] values;

	public Factor(FactorScope aScope) {
		this.scope = aScope;

		int numberValues = 1;

		for (int i = 0; i < scope.size(); i++)
			numberValues *= scope.getVariable(i).cardinality();

		values = new double[numberValues];

	}

	/**
	 * Create a factor for a list of variables. Note that the order of variable
	 * in the scope is not neccesarily the order specified in the varList.
	 * 
	 * @param varList
	 */
	public Factor(DiscreteVariable[] varList) {
		this(new FactorScope(varList));
	}

	/**
	 * Convert an assignment of values to index of the factor
	 * 
	 * @return
	 */
	private int assignment2Index(int[] assignment) {
		int ret = 0;
		int size = 1;
		for (int i = scope.size() - 1; i >= 0; i--) {
			ret += size * assignment[i];
			size *= scope.getVariable(i).cardinality();
		}
		return ret;
	}

	/**
	 * Turn the index of the factor to a list of assignment
	 * 
	 * @param index
	 * @return
	 */
	private int[] index2Assignment(int index) {
		int[] assignment = new int[scope.size()];

		for (int i = scope.size() - 1; i >= 0; i--) {
			assignment[i] = index % scope.getVariable(i).cardinality();
			index /= scope.getVariable(i).cardinality();
		}
		return assignment;
	}

	// private double getValue(int index){
	// return values[index];
	// }

	private double getValue(int[] assignment) {
		return values[assignment2Index(assignment)];
	}

	/**
	 * Return the product of this factor with another
	 * 
	 * @param f
	 * @return
	 */
	public Factor multiply(Factor f) {
		return product(this, f);
	}

	/**
	 * Return factor product of the two factors
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static Factor product(Factor a, Factor b) {

		// First find the unions of factor scopes: set the bitset
		FactorScope newScope = FactorScope.union(a.scope, b.scope);
		Factor prod = new Factor(newScope);

		AssignmentIterator iter = prod.iterator();
		// //////////////////////////////////////////////////////////////////////
		// FIXME: the fastest algorithm would be
		// 1. Get the intersection of the two factors a and b
		// 2. Get an iterator of the intersection
		// 3. For each assignment of the intersection,
		// 3a. Get an iterator of assignment of a that includes the intersection
		// assignment
		// 3b. For each assignment of a (that includes the intersection
		// assignment)
		// 3bi. For each assignment of b that includes the intersection
		// assignment
		// 3bii. Get the producnt assignment, and put in the values
		// /////////////////////////////////////////////////////////////////////

		int index = 0;
		while (iter.hasNext()) {
			int[] assignment = iter.next();
			// Get assignment for A
			int[] aAssignment = prod.assignmentReduce(assignment, a.scope);
			int[] bAssignment = prod.assignmentReduce(assignment, b.scope);

			// get assignment for B

			// assert index = assignment2Index(assignment)
			prod.values[index] = a.getValue(aAssignment)
					* b.getValue(bAssignment);

			index++;
		}
		return prod;
	}

	/**
	 * Reduce an assignment to another assignment of a scope
	 * 
	 * @param assignment
	 * @param oScopes
	 * @return
	 */
	private int[] assignmentReduce(int[] assignment, FactorScope aScope) {
		int[] oAssignment = new int[aScope.size()];

		int index = 0;
		for (int i = 0; i < this.scope.size(); i++) {
			if (aScope.getVariable(index) == this.scope.getVariable(i)) {
				oAssignment[index] = assignment[i];
				index++;
				if (index >= oAssignment.length)
					break;
			}
		}
		return oAssignment;
	}

	/**
	 * Marginalise the factor to another scope
	 * 
	 * @param oScopes
	 * @return
	 */
	public Factor marginalise(FactorScope aScope) {
		Factor aFactor = new Factor(aScope);
		int myIndex = 0;

		AssignmentIterator iter = this.iterator();
		while (iter.hasNext()) {
			int[] assignment = iter.next();
			int[] reducedAssignment = assignmentReduce(assignment, aScope);
			int oIndex = aFactor.assignment2Index(reducedAssignment);
			// System.out.println(myIndex + "  " + oIndex);
			aFactor.values[oIndex] += values[myIndex];
			myIndex++;
		}

		return aFactor;
	}

	AssignmentIterator iterator() {
		return new AssignmentIterator();
	}

	public void print(PrintStream out) {
		out.print("         ");
		for (int i = 0; i < scope.size(); i++) {
			out.printf("%s  ",
					(scope.getVariable(i) + "    ").subSequence(0, 3));
		}
		out.println();
		for (int i = 0; i < values.length; i++) {
			int[] assignment = index2Assignment(i);
			out.printf("%4d : ", i);
			for (int j = 0; j < assignment.length; j++)
				out.printf("%4d ", assignment[j]);
			out.printf(" =  %8.6f  \n", values[i]);
		}
	}

	public void putValue(DiscreteVariable[] vList, int[] valList, double val) {
		int[] assignment = makeAssignment(scope, vList, valList);
		values[assignment2Index(assignment)] = val;
	}

	public void putValue(DiscreteVariable[] vList, String[] strList, double val) {
		int[] assignment = makeAssignment(scope, vList, strList);
		values[assignment2Index(assignment)] = val;
	}

	public static int[] makeAssignment(FactorScope scope,
			DiscreteVariable[] vList, int[] valList) {
		int[] assignment = new int[scope.size()];

		for (int i = 0; i < assignment.length; i++) {
			for (int j = 0; j < vList.length; j++) {
				if (scope.getVariable(i) == vList[j])
					assignment[i] = valList[j];
			}
		}
		return assignment;
	}

	public static int[] makeAssignment(FactorScope scope,
			DiscreteVariable[] vList, String[] strList) {
		int[] valList = new int[strList.length];

		for (int i = 0; i < vList.length; i++) {
			for (int j = 0; j < vList[i].cardinality(); j++) {
				if (vList[i].getState(j).equals(strList[i]))
					valList[i] = j;
			}
		}

		return makeAssignment(scope, vList, valList);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		/**********************************************************
		 * System.out.println("Hello World");
		 * 
		 * DiscreteVariable vC = new DiscreteVariable("coherence",new String[]
		 * {"true", "false"}); DiscreteVariable vD = new
		 * DiscreteVariable("difficulty",new String[] {"easy", "moderate",
		 * "hard"});
		 * 
		 * DiscreteVariable vI = new DiscreteVariable("intelligence",new
		 * String[] {"true", "false"}); DiscreteVariable vG = new
		 * DiscreteVariable("grade",new String[] {"true", "false"});
		 * DiscreteVariable vS = new DiscreteVariable("sat",new String[]
		 * {"true", "false"}); DiscreteVariable vL = new
		 * DiscreteVariable("letter",new String[] {"true", "false"});
		 * DiscreteVariable vJ = new DiscreteVariable("job",new String[]
		 * {"true", "false"}); DiscreteVariable vH = new
		 * DiscreteVariable("happy",new String[] {"true", "false"});
		 * 
		 * Factor fC = new Factor(new FactorScope (vC)); Factor fDcC = new
		 * Factor(new DiscreteVariable[] {vC, vD});
		 * 
		 * double [] pC = {.4,.6}; fC.values = pC; double [] pDcC = {.3, .3, .4,
		 * .5, .35, .15}; fDcC.values = pDcC;
		 * 
		 * Factor fCD = product(fC, fDcC);
		 * 
		 * Factor fD = fCD.marginalise(new FactorScope (vD));
		 * 
		 * //System.out.println("P(c)"); //fC.print(System.out);
		 * 
		 * // System.out.println("P(d|c)"); // fDcC.print(System.out);
		 * 
		 * 
		 * //System.out.println("P(dc)"); //fCD.print(System.out);
		 * 
		 * 
		 * System.out.println("P(d)");
		 * 
		 * fD.print(System.out);
		 * 
		 * 
		 * DiscreteVariable l0 = new DiscreteVariable("l0", new String[]
		 * {"true", "false"}); DiscreteVariable r0 = new DiscreteVariable("r0",
		 * new String[] {"true", "false"});
		 * 
		 * DiscreteVariable ri = r0;
		 * 
		 * Factor fL0 = new Factor(new FactorScope (l0)); Factor fR0 = new
		 * Factor(new FactorScope (r0));
		 * 
		 * fL0.values[0] = 0.4; fL0.values[1] = 0.6;
		 * 
		 * fR0.values[0] = 0.7; fR0.values[1] = 0.3;
		 * 
		 * //Factor fLi = fL0, fRi = fR0;
		 * 
		 * //fLi.print(System.out); //fRi.print(System.out);
		 * 
		 * int n = 25; if (args.length > 0) n = Integer.parseInt(args[0]);
		 * 
		 * Factor jdt = product(fL0, fR0); for (int i = 1; i <= n; i++){
		 * DiscreteVariable vli = new DiscreteVariable("l"+(i), new String[]
		 * {"true", "false"}); DiscreteVariable vri = new
		 * DiscreteVariable("r"+(i), new String[] {"true", "false"});
		 * 
		 * //Create CPT DiscreteVariable[] vlList = new DiscreteVariable[] {vli,
		 * (DiscreteVariable) Variable.getVariable(("l"+(i-1))),
		 * (DiscreteVariable) Variable.getVariable(("r"+(i-1)))};
		 * DiscreteVariable[] vrList = new DiscreteVariable[] {vri,
		 * (DiscreteVariable) Variable.getVariable(("l"+(i-1))),
		 * (DiscreteVariable) Variable.getVariable(("r"+(i-1)))};
		 * 
		 * Factor condLi = new Factor(vlList); Factor condRi = new
		 * Factor(vrList);
		 * 
		 * //li.put(new Boolean[]{true,true}, 0.25); condLi.putValue(vlList, new
		 * String[]{"true","true","true"}, 0.25); condLi.putValue(vlList, new
		 * String[]{"false","true","true"}, 0.75);
		 * 
		 * //li.put(new Boolean[]{true,false}, 0.35); condLi.putValue(vlList,
		 * new String[]{"true","true","false"}, 0.35); condLi.putValue(vlList,
		 * new String[]{"false","true","false"}, 0.65);
		 * 
		 * //li.put(new Boolean[]{false,true}, 0.45); condLi.putValue(vlList,
		 * new String[]{"true","false","true"}, 0.45); condLi.putValue(vlList,
		 * new String[]{"false","false","true"}, 0.55);
		 * 
		 * //li.put(new Boolean[]{false,false}, 0.55); condLi.putValue(vlList,
		 * new String[]{"true","false","false"}, 0.55); condLi.putValue(vlList,
		 * new String[]{"false","false","false"}, 0.45);
		 * 
		 * 
		 * //ri.put(new Boolean[]{true,true}, 0.65); condRi.putValue(vrList, new
		 * String[]{"true","true","true"}, 0.65); condRi.putValue(vrList, new
		 * String[]{"false","true","true"}, 0.35);
		 * 
		 * //ri.put(new Boolean[]{true,false}, 0.5); condRi.putValue(vrList, new
		 * String[]{"true","true","false"}, 0.5); condRi.putValue(vrList, new
		 * String[]{"false","true","false"}, 0.5);
		 * 
		 * //ri.put(new Boolean[]{false,true}, 0.4); condRi.putValue(vrList, new
		 * String[]{"true","false","true"}, 0.4); condRi.putValue(vrList, new
		 * String[]{"false","false","true"}, 0.6);
		 * 
		 * //ri.put(new Boolean[]{false,false}, 0.3); condRi.putValue(vrList,
		 * new String[]{"true","false","false"}, 0.3); condRi.putValue(vrList,
		 * new String[]{"false","false","false"}, 0.7);
		 * 
		 * //Factor prod = product(fLi, fRi); //System.out.printf(
		 * "=================== P(l%d, r%d) ==========================\n"
		 * ,i-1,i-1); //prod.print(System.out); //Factor newfL =
		 * product(product(fLi, fLi), fRi).marginalise(new FactorScope(vli));
		 * 
		 * //Factor newfL = product(condLi, prod); //newfL =
		 * newfL.marginalise(new FactorScope(vli));
		 * 
		 * //Factor newfR = product(condRi, prod); //newfR =
		 * newfR.marginalise(new FactorScope(vri));
		 * 
		 * //fLi = newfL; //fRi = newfR;
		 * 
		 * //fLi.print(System.out); //fRi.print(System.out); jdt = product(jdt,
		 * condLi); jdt = product(jdt, condRi);
		 * 
		 * jdt = jdt.marginalise(new FactorScope(new
		 * DiscreteVariable[]{vli,vri})); ri = vri;
		 * 
		 * jdt.marginalise(new FactorScope(ri)).print(System.out); }
		 * //fRi.print(System.out);
		 * /**********************************************************
		 */

	}

	class AssignmentIterator implements Iterator<int[]> {

		int currentIndex = 0;

		AssignmentIterator() {

		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Iterator#hasNext()
		 */
		@Override
		public boolean hasNext() {
			return currentIndex < values.length;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Iterator#next()
		 */
		@Override
		public int[] next() {
			// FIXME: this is a slow implementation of the iterator
			int[] assignment = index2Assignment(currentIndex);
			currentIndex++;

			return assignment;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Iterator#remove()
		 */
		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}
	}
}

/**
 * Implement the factor scope as a set of variables, and the variables are kept
 * in order.
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 * 
 */
class FactorScope {
	DiscreteVariable[] varList;

	/**
	 * Create the scope factor of one variable
	 * 
	 * @param var
	 */
	FactorScope(DiscreteVariable var) {
		varList = new DiscreteVariable[] { var };
	}

	FactorScope(DiscreteVariable[] vList) {
		varList = new DiscreteVariable[vList.length];
		BitSet bSet = new BitSet();
		for (int i = 0; i < vList.length; i++) {
			bSet.set(vList[i].getIndex());
		}
		int index = varList.length - 1;

		for (int i = bSet.length(); i >= 0; i--) {
			if (bSet.get(i)) {
				varList[index] = (DiscreteVariable) Variable.getVariable(i);
				index--;
			}
		}
	}

	private FactorScope(int size) {
		varList = new DiscreteVariable[size];
	}

	public int size() {
		return varList.length;
	}

	public DiscreteVariable getVariable(int i) {
		return varList[i];
	}

	static FactorScope union(FactorScope a, FactorScope b) {
		BitSet bset = new BitSet();
		for (int i = 0; i < a.size(); i++) {
			bset.set(a.varList[i].getIndex());
		}

		for (int i = 0; i < b.size(); i++) {
			bset.set(b.varList[i].getIndex());
		}

		FactorScope uScope = new FactorScope(bset.cardinality());

		int idx = uScope.size() - 1;
		for (int i = bset.length() - 1; i >= 0; i--) {
			if (bset.get(i)) {
				uScope.varList[idx] = (DiscreteVariable) Variable
						.getVariable(i);
				idx--;
			}
		}
		return uScope;
	}
}
