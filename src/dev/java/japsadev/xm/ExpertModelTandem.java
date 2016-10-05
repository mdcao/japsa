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
 * 10/01/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.xm;

import java.util.Iterator;
import java.util.LinkedList;

import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.util.CommandLine;
import japsa.util.JapsaMath;
import japsa.xm.ExpertModel;
import japsa.xm.expert.AdaptiveMarkovExpert;
import japsa.xm.expert.Expert;
import japsa.xm.expert.MarkovExpert;
import japsadev.xm.expert.*;

public class ExpertModelTandem {
	// protected CombinationExpert baseEx;
	MarkovExpert markovEx;
	AdaptiveMarkovExpert adapMarkovEx;//
	// = new MarkovExpert(null,2);

	// Weighted probability distribution of the next character over all models.
	double[] finalD;// ,baseD;

	// W'ghted probability distribution of the next character over both markovs.
	double[] markovD;

	int currentInd;

	// The cost incurred by the expert for not being certain of the correct
	// symbol.
	double[] markovCost = null;
	// double[][] tandemCosts; // tandems of length 3,4,5

	int maxPeriod = 9;

	protected double repeatPrior;
	protected double repeatPriorProb;

	public static double threshold;
	public static int minLen = 10;
	public static JapsaAnnotation anno;

	// Two expert seeds
	// public TandemExpert tandemEx;
	LinkedList<TandemExpert> panel;

	/*
	 * Parameters of the algorithm
	 */

	public ExpertModelTandem(Alphabet alphabet, int context,
			double listenThreshold) {
		this(alphabet, context, listenThreshold, 9);
	}

	public ExpertModelTandem(Alphabet alphabet, int context,
			double listenThreshold, int max) {
		super();

		// this.hashSize = hashSize;
		Expert.setParams(alphabet.size(), context);
		this.repeatPrior = listenThreshold * context;
		repeatPriorProb = JapsaMath.exp2(-repeatPrior);

		finalD = new double[alphabet.size()];
		markovD = new double[alphabet.size()];// Combination of the markov

		maxPeriod = max;
		threshold = listenThreshold;
	}

	public void printParams() {
		System.out.println("Parameters:" + "\nContext          : "
				+ Expert.CONTEXT_LENGTH + "\nListen Threshold : " + repeatPrior
				/ Expert.CONTEXT_LENGTH + "bps");
	}

	protected void initialiseCommon(Sequence seq) {

		anno = new JapsaAnnotation(seq);
		// Background knowledge/based knowledge
		markovEx = new MarkovExpert(2, 100);
		adapMarkovEx = new AdaptiveMarkovExpert(1, 256, 100);

		markovCost = new double[seq.length()];

		for (int i = 0; i < maxPeriod; i++) {
			panel.add(new TandemExpert(seq, i + 1, Expert.CONTEXT_LENGTH));
		}
	}

	/*****************************************************************/

	protected void preCoding() {
		// Get ranking of experts
		double baseRateProb = markovEx.posteriorProb();
		double finalSum = markovEx.posteriorProb()
				+ adapMarkovEx.posteriorProb();

		for (byte a = 0; a < Expert.alphabet().size(); a++) {
			// baseD[a] =
			markovD[a] = finalD[a] = markovEx.posteriorProb()
					* markovEx.probability(a) + adapMarkovEx.posteriorProb()
					* adapMarkovEx.probability(a);

			markovD[a] /= finalSum;
		}

		// Listen to no more than expert limit experts
		Iterator<TandemExpert> panelIter = panel.iterator();

		while (panelIter.hasNext()) {// Inv: head.next = ptr
			// double rate = ptr.rate();//msg leng over an history
			TandemExpert ptr = panelIter.next();
			double score = ptr.posteriorProb();
			// Only listen if the offset/palindrome is sufficiently good
			if (score > baseRateProb / repeatPriorProb) {
				for (byte a = 0; a < Expert.alphabet().size(); a++) {
					finalD[a] += ptr.probability(a) * score;
				}
				// ptr.weight = myCost;
				finalSum += score;
			}
			// ptr = ptr.getNext();
		}

		// Normalise
		for (byte a = 0; a < Expert.alphabet().size(); a++) {
			finalD[a] /= finalSum;
		}
	}

	/**
	 * Update all experts at this positition
	 * 
	 * @param seqArray
	 * @param i
	 * @param sid
	 */
	protected void updateExperts(byte c) {
		markovEx.update(c);
		adapMarkovEx.update(c);

		/********************* Update Exp ********************/
		Iterator<TandemExpert> panelIter = panel.iterator();

		while (panelIter.hasNext()) {
			TandemExpert ptr = panelIter.next();
			ptr.update(c);// This actually the prob
			// ptr = ptr.getNext();
		}
	}

	protected void updateTandemExpertsAtPos(int pos) {

		Iterator<TandemExpert> panelIter = panel.iterator();
		while (panelIter.hasNext()) {
			panelIter.next().updatePos(pos);
		}
	}

	protected void postCoding(Sequence seqArray) {
	}

	/** Get the cost array for the tandem expert. */
	public void calculateTandemCosts(Sequence seq) {
		initialiseCommon(seq);

		for (currentInd = 0; currentInd < seq.length(); currentInd++) {
			updateTandemExpertsAtPos(currentInd);
		}
	}

	/**
	 * Done only after calculate Tandem cost
	 */
	public void filtering() {
	}

	// double gain;
	public double[] encode(Sequence seqArray) {

		initialiseCommon(seqArray);

		// Get the sequence to be encode
		// byte japsa.seq[] = seqArray.toBytes();
		double[] costs = new double[seqArray.length()];

		// double markovTotal = 0;
		// double totalCost = 0.0;
		// go thru the sequence and encode each charactor
		for (currentInd = 0; currentInd < seqArray.length(); currentInd++) {
			// Compute the probability distribution
			preCoding();

			int actual = seqArray.getBase(currentInd);
			// double cost = -JapsaMath.log2(finalD[actual]);

			markovCost[currentInd] = -JapsaMath.log2(markovD[actual]);

			// costs[currentInd] = cost;
			// totalCost += cost;

			// markovTotal += markovCost[currentInd];
			updateExperts(seqArray.getBase(currentInd));
		}

		// gain = (totalCost - markovTotal) / seqArray.length();

		// System.out.print(gain + "  " + (markovTotal / seqArray.length()) +
		// "   ");

		return costs;
	}

	public static CommandLine prepareCmd() {
		CommandLine cmdLine = new CommandLine();
		cmdLine.addInt("context", 15, "Length of the context");
		cmdLine.addDouble("threshold", 0.15, "Listen threshold");

		return cmdLine;
	}

	public static void main11(String[] args) throws Exception {

		// Get params from users

		CommandLine cmdLine = prepareCmd();

		args = cmdLine.parseLine(args);
		System.out.println(ExpertModel.version());

		if (args == null || args.length <= 0) {
			System.err
					.println("Usage: java CommandLine [options]  file1 file2 ...\n"
							+ cmdLine.usageMessage() + "\n");
			System.exit(1);
		}

		ExpertModelTandem tModel = new ExpertModelTandem(Alphabet.DNA4(), 5,
				cmdLine.getDoubleVal("threshold"));

		// Print out all params
		tModel.printParams();

		// cmdLine.printOptions();

		FastaReader fReader = new FastaReader(args[0]);

		Sequence seq;

		while ((seq = fReader.nextSequence(Alphabet.DNA4())) != null) {
			System.out.print("\n" + seq.getName() + " :  ");

			// long start = System.currentTimeMillis();
			// seqHash[1] = japsa.seq;
			// double[] costs =
			tModel.encode(seq);// ,args[1]);

			// long time = (System.currentTimeMillis() - start);

			// System.out.printf("Compress in %d ms\n", time);
			// IOTools.writeDoubleSequence("cost.info", "Cost", costs);
			// IOTools.writeDoubleSequence("markov.info", "Markov Cost",
			// tModel.markovCost);
			// for (int idx = 0; idx < tModel.tandemCosts.length; idx++) {
			// IOTools.writeDoubleSequence("tandem" + (idx + 2) + ".info",
			// "Cost of tandem " + (idx + 2), tModel.tandemCosts[idx]);
			// }
			// System.out.println(tModel.gain);

			System.out
					.println("=============================================================================");

			/*************************************************************************/

		}
		fReader.close();
	}

}
