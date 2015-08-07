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
 * 10/01/2012 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsa.xm;

import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.LinkedList;

import japsa.seq.AbstractSequence;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.tools.xm.ExpertModelCmd;
import japsa.util.IntIterator;
import japsa.util.MyBitSet;
import japsa.util.JapsaMath;
import japsa.xm.expert.AdaptiveMarkovExpert;
import japsa.xm.expert.CombinationExpert;
import japsa.xm.expert.Expert;
import japsa.xm.expert.MarkovExpert;
import japsa.xm.expert.RepeatExpert;
import japsa.xm.hash.GappedHashtable;
import japsa.xm.hash.MyBinaryHashtable;
import japsa.xm.hash.MyHashtable;
import japsa.xm.hash.PatternStore;

import com.colloquial.arithcode.ArithDecoder;
import com.colloquial.arithcode.ArithEncoder;
import com.colloquial.arithcode.BitInput;

public class ExpertModel {
	public static String VERSION = "V 3.0";
	private int checkPoint = 1000000;

	/**
	 * @return the checkPoint
	 */
	public int getCheckPoint() {
		return checkPoint;
	}

	/**
	 * @param checkPoint the checkPoint to set
	 */
	public void setCheckPoint(int checkPoint) {
		this.checkPoint = checkPoint;
	}

	// 'global' Variables
	protected int expertCount;// The number of experts currently empliyed
	protected int[] accLengths;// lengths of sequences (including context)

	// protected CombinationExpert baseEx;//Any fixed expert are in this linked
	// list
	MarkovExpert markovEx;
	AdaptiveMarkovExpert adapMarkovEx;//
	// = new MarkovExpert(null,2);
	double[] finalD;// ,baseD;
	double[] markovD;
	public String hashType = "hash";
	int currentInd;

	public static double[] markovCost = null;

	protected CombinationExpert repEx;

	protected double repeatPrior;
	protected double repeatPriorProb;

	// protected Distribution finalDist;
	protected int posSize = 0;
	protected int lengthSeqs = 0;// Total lengths of all sequences
	protected int total = 1024 * 1024;
	MyBitSet bitSet, pBitSet;// To store if an particular offset has been stored

	// Two expert seeds
	public RepeatExpert offSetSeed; // = new OffsetCountExpert(null,1,null);// a
									// dummy one
	public RepeatExpert palinSeed;// = new PalindromeCountExpert(null,1,null);//
									// a dummy one
	/*
	 * Parameters of the algorithm
	 */
	protected int expertsLimit;// Limit of experts
	protected int chances;// Number of chances given to each expert
	PatternStore myHash; // The hashtable
	int hashSize;
	
	boolean selfRep = false;
	boolean binaryHash = false;
	
	LinkedList<RepeatExpert> panel;
	Alphabet alphabet;

	public ExpertModel(int hashSize, Alphabet alphabet , int context,
			int expertsLimit, double listenThreshold, int chances,
			boolean binaryHash) {
		super();
		Expert.setAlphabet(alphabet);
		this.alphabet = alphabet;

		// this.hashSize = hashSize;
		Expert.setContext(context);
		this.expertsLimit = expertsLimit;
		this.chances = chances;
		this.repeatPrior = listenThreshold * context;
		repeatPriorProb = JapsaMath.exp2(-repeatPrior);

		// thExp = JapsaMath.exp2(listenThreshold);
		// this.listenThreshold = Math.pow(2,listenThreshold);

		// finalDist = new Distribution(Expert.alphabet().size());

		finalD = new double[Expert.alphabet().size()];
		markovD = new double[Expert.alphabet().size()];// Combination of the markov
													// experts
		this.hashSize = hashSize;
		this.binaryHash = binaryHash;
	}

	public void setBinaryHash(boolean useBinary) {
		binaryHash = useBinary;
	}

	public void setSelfRep(boolean rep) {
		selfRep = rep;
	}

	public void setHashType(String hashT){
		// Normal hash
		if ("hash".equals(hashT))
			this.hashType = "hash";
		// Suffix array
		else if ("sfa".equals(hashT))
			this.hashType = "" + "sfa";
		// suffix tree
		else if ("sft".equals(hashT))
			this.hashType = "" + "sft";
		// gapped
		else {
			for (int i = 0; i < hashT.length(); i++) {
				char c = hashT.charAt(i);
				if (c != '0' && c != '1') {
					throw new RuntimeException("Unknown hash type : " + hashT);
				}
			} // assert: hashT ok

			int i = hashT.length() - 1;
			while (i >= 0 && hashT.charAt(i) == '0')
				i--;
			if (i < 0)
				throw new RuntimeException("Unknown hash type : " + hashT);
			this.hashType = hashT.substring(hashT.indexOf('1'), i + 1);
			System.out.println(hashT.indexOf('1') + " vs " + i);

			this.hashSize = hashType.length();

		}
	}

	public void printParams() {
		String hashName = "Hashtable";

		if ("hash".equals(hashType)) {

		} else if ("sfa".equals(hashType))
			hashName = "PrefixArray";
		else if ("sft".equals(hashType)) {
			hashName = "PrefixTree";
		} else {
			hashName = "Gapped(" + hashType + ")";
		}

		System.out.println("Parameters:" + "\nHash size        : " + hashSize
				+ "\nExpert Limit     : " + expertsLimit
				+ "\nContext          : " + Expert.CONTEXT_LENGTH
				+ "\nListen Threshold : " + repeatPrior / Expert.CONTEXT_LENGTH
				+ "bps" + "\nChances          : " + chances
				+ "\nBinaryHash       : " + binaryHash
				+ "\nHashType         : " + hashName + "\nExpert Type      : "
				+ offSetSeed.getClass());
	}

	// protected abstract void initilise(BioCompSequence[] seqArray);//should
	// store

	protected void initiliseGappedHash(AbstractSequence[] seqArray) {
		/*************************************************************************/
		myHash = new GappedHashtable(hashType);
		/*************************************************************************/

		store(seqArray);
	}

	protected void initiliseHash(AbstractSequence[] seqArray) {
		/*************************************************************************/
		if (binaryHash) {
			myHash = new MyBinaryHashtable(hashSize, (int) Math.ceil(JapsaMath.log2(alphabet.size())));
		} else
			myHash = new MyHashtable(hashSize, (int) Math.ceil(JapsaMath.log2(alphabet.size())));
		/*************************************************************************/

		store(seqArray);
	}
	
	public void store(AbstractSequence[] seqArray) {
		// Store all back ground sequence in the hash
		for (int sid = 0; sid < seqArray.length - 1; sid++) {
			for (int i = 0; i < seqArray[sid].length(); i++) {
				myHash.nextKey(seqArray[sid].symbolAt(i));
				myHash.putCurrentValue((sid << posSize) + i);
			}
		}
	}

	protected void initilise_optimise(AbstractSequence[] seqArray) {
		// Common stuff
		lengthSeqs = 0;
		accLengths = new int[seqArray.length + 1];
		accLengths[0] = 0;
		for (int i = 0; i < seqArray.length; i++) {
			lengthSeqs += seqArray[i].length();
			accLengths[i + 1] = lengthSeqs;
		}

		expertCount = 0;

		int numSeqs = 1;

		if (seqArray.length > 1)
			numSeqs = ((int) JapsaMath.log2(seqArray.length - 1)) + 1;

		posSize = 31 - numSeqs;

		// Background knowledge/based knowledge

		markovEx = new MarkovExpert(2);
		adapMarkovEx = new AdaptiveMarkovExpert(1,256);
		//markovEx.setNext(adapMarkovEx);

		
		repEx = new CombinationExpert();
		panel = new LinkedList<RepeatExpert>();
		
		MyHashtable hash = new MyHashtable(hashSize, (int) Math.ceil(JapsaMath.log2(alphabet.size())));
		/*************************************************************************/
		
		System.out.println("Run first pass");

		for (int sid = 0; sid < seqArray.length; sid++) {
			for (int i = 0; i < seqArray[sid].length(); i++) {
				hash.nextKey(seqArray[sid].symbolAt(i));
				hash.putCurrentValue_psuedo((sid << posSize) + i);
			}
		}

		hash.reinitialise_optimise();
		System.out.println("Finish first pass");

		myHash = hash;
		store(seqArray);

		checkPoint(0);

		bitSet = new MyBitSet(lengthSeqs);
		pBitSet = new MyBitSet(lengthSeqs * 2);

	}

	protected void initiliseCommon(AbstractSequence[] seqArray) {
		// Common stuff
		lengthSeqs = 0;
		accLengths = new int[seqArray.length + 1];
		accLengths[0] = 0;
		for (int i = 0; i < seqArray.length; i++) {
			lengthSeqs += seqArray[i].length();
			accLengths[i + 1] = lengthSeqs;
		}

		expertCount = 0;

		int numSeqs = 1;

		if (seqArray.length > 1)
			numSeqs = ((int) JapsaMath.log2(seqArray.length - 1)) + 1;

		posSize = 31 - numSeqs;
		// Background knowledge/based knowledge

		markovEx = new MarkovExpert(2);
		adapMarkovEx = new AdaptiveMarkovExpert(1, 256);
		
		panel = new LinkedList<RepeatExpert>();
		repEx = new CombinationExpert();

		if ("hash".equals(this.hashType))
			initiliseHash(seqArray);
		else {
			initiliseGappedHash(seqArray);
		}

		bitSet = new MyBitSet(lengthSeqs);
		pBitSet = new MyBitSet(lengthSeqs * 2);

	}

	/*****************************************************************/
	/**
	 * Precoding when the next symbol is known (used in encoding)
	 * @param nextSym
	 */
	protected void preCoding(byte nextSym) {
		// double baseRateProb = (markovEx.rateProb() + adapMarkovEx.rateProb())
		// / 2; //MsgLenth of markov expert
		double baseRateProb = markovEx.posteriorProb();

		double finalSum = baseRateProb + adapMarkovEx.posteriorProb();

		markovD[nextSym] = finalD[nextSym] = markovEx.posteriorProb()
				* markovEx.probability(nextSym) + adapMarkovEx.posteriorProb()
				* adapMarkovEx.probability(nextSym);
		markovD[nextSym] /= finalSum;


		// Initialise expert prediction and total prediction
		repEx.getCombDistribution().setWeights(0.0);
		double repSum = 0.0;

		Iterator<RepeatExpert> panelIter = panel.iterator();		
		while (panelIter.hasNext()) {
			RepeatExpert ptr = panelIter.next();

			double score = ptr.posteriorProb(); // costToWeight(rate);

			if (score > baseRateProb / repeatPriorProb){
				repEx.getCombDistribution().addWeight(nextSym,	score * ptr.probability(nextSym));
				
				repSum += score;
				ptr.resetCounter();				
			} else {
				ptr.incrementCounter();

				if (ptr.getCounter() > chances) {
					resignExpert(ptr);// ptr.resign();
					expertCount--;
					panelIter.remove();
				}
			}
		}

		if (repSum != 0.0) {// There is repeat			
			repEx.getCombDistribution().scale(1.0 / repSum);

			double repExScore = repEx.posteriorProb() * repeatPriorProb;

			// Combine 2 experts
			for (int a = 0; a < Expert.alphabet().size(); a++) {
				finalD[a] += repEx.getCombDistribution().getWeight(a)
						* repExScore;
			}
			finalSum += repExScore;
		} else {
			repEx.getCombDistribution().setWeights(1.0 / Expert.alphabet().size());
		}

		for (int a = 0; a < Expert.alphabet().size(); a++) {
			finalD[a] /= finalSum;
		}
	}


	protected void preCoding() {
		// double baseRateProb = (markovEx.rateProb() + adapMarkovEx.rateProb())
		// / 2; //MsgLenth of markov expert
		double baseRateProb = markovEx.posteriorProb();

		double finalSum = markovEx.posteriorProb() + adapMarkovEx.posteriorProb();

		for (int a = 0; a < Expert.alphabet().size(); a++) {
			// baseD[a] =
			markovD[a] = finalD[a] = markovEx.posteriorProb()
					* markovEx.probability(a) + adapMarkovEx.posteriorProb()
					* adapMarkovEx.probability(a);

			markovD[a] /= finalSum;
		}

		// Initialise expert prediction and total prediction
		repEx.getCombDistribution().setWeights(0.0);
		double repSum = 0.0;

		Iterator<RepeatExpert> panelIter = panel.iterator();		
		while (panelIter.hasNext()) {
			RepeatExpert ptr = panelIter.next();

			double score = ptr.posteriorProb(); // costToWeight(rate);

			if (score > baseRateProb / repeatPriorProb){														
				for (int a = 0; a < Expert.alphabet().size(); a++) {
					repEx.getCombDistribution().addWeight(a,
							score * ptr.probability(a));
				}
				
				repSum += score;
				ptr.resetCounter();				
			} else {
				ptr.incrementCounter();

				if (ptr.getCounter() > chances) {
					resignExpert(ptr);// ptr.resign();
					expertCount--;
					panelIter.remove();
				}
			}
		}

		if (repSum != 0.0) {// There is repeat			
			repEx.getCombDistribution().scale(1.0 / repSum);

			double repExScore = repEx.posteriorProb() * repeatPriorProb;

			// Combine 2 experts
			for (int a = 0; a < Expert.alphabet().size(); a++) {
				finalD[a] += repEx.getCombDistribution().getWeight(a)
						* repExScore;
			}
			finalSum += repExScore;
		} else {
			repEx.getCombDistribution().setWeights(1.0 / Expert.alphabet().size());
		}

		for (int a = 0; a < Expert.alphabet().size(); a++) {
			finalD[a] /= finalSum;
		}
	}

	protected void resurrectExpert(RepeatExpert e, AbstractSequence bs, int pos, int past) {
		e.resurrect(bs, pos, past);
	}

	/**
	 * Update all experts at this positition
	 * 
	 * @param seqArray
	 * @param i
	 * @param sid
	 */
	protected void updateExperts(int c) {
		markovEx.update(c);
		adapMarkovEx.update(c);

		repEx.update(c);
		/********************* Update Exp ********************/

		Iterator<RepeatExpert> panelIter = panel.iterator();
		while (panelIter.hasNext()) {
			RepeatExpert ptr = panelIter.next();
			double thisCost = ptr.update(c);// This actually the prob
			if (thisCost < 0) {// A cost should be a positive number
				resignExpert(ptr);// ptr.resign();
				expertCount--;
				panelIter.remove();
			}
		}
	}

	void resignExpert(RepeatExpert p) {
		p.resign();
	}

	protected void postCoding(AbstractSequence[] seqArray, int sid) {
		/**********************************************************************/
		// Move to next key
		// if (currentInd < seqArray[sid].length() - 1)
		myHash.nextKey(seqArray[sid].symbolAt(currentInd));

		// Add good diagonals.
		if (currentInd >= hashSize && currentInd < seqArray[sid].length() - 1) {
			// assert n >= 0
			// IntIterator iter = myHash.copyIterator();
			IntIterator iter = myHash.iterator();

			while (expertsLimit > expertCount && iter.hasNext()) {
				int position = iter.next();
				RepeatExpert e = null;

				if (position > 0) {// Offset expert
					int id = position >> posSize;
					int pos = position % (1 << posSize);

					if (pos <= hashSize) {
						continue;
					}

					if (pos > hashSize && // Have enough for resurrect
							!bitSet.get(currentInd + accLengths[sid]
									- accLengths[id] - pos)// Not in there
							&& pos < seqArray[id].length() - 3) {// have some
																	// thing to
																	// predict
						e = offSetSeed.duplicate(seqArray[id], pos,
								bitSet);
						e.setID(currentInd + accLengths[sid] - accLengths[id]
								- pos);
					}
				} else {// Palindrome expert
					position = -position;
					// position = position - hashSize;
					int id = position >> posSize;
					int pos = position % (1 << posSize);

					if (pos > hashSize && // Have enough for resurrect
							(!pBitSet.get(currentInd + accLengths[sid]
									+ accLengths[id] + pos))// ?
							&& (pos + 3 < seqArray[id].length())) {// Have
																			// something
																			// to
																			// predict
						e = palinSeed.duplicate(seqArray[id], pos
								- hashSize + 1, pBitSet);
						e.setID(currentInd + accLengths[sid] + accLengths[id]
								+ pos);
					}
				}
				// Add this expert in only if an identical expert not in the
				// list
				if (e != null) {
					resurrectExpert(e, seqArray[sid], currentInd,
							hashSize);
					// e.resurrect(seqArray[sid].toBytes(),currentInd,hashSize);
					panel.add(e);
					//e.setNext(repEx.getNext());
					//repEx.setNext(e);
					expertCount++;
				}
			}
		}
		// System.out.println(i + "   " + expertCount);
		/********************** Store current *********************/
		if (selfRep)
			myHash.putCurrentValue((sid << posSize) + currentInd);
		/******************************************************************/
	}

	protected void checkPoint(int steps) {
		System.out.print("Reach milestone " + steps + " : ");
		myHash.printSummary();

		System.out.println("  Memory availabe "
				+ Runtime.getRuntime().freeMemory() + "     "
				+ Runtime.getRuntime().totalMemory());
		System.gc();
		System.out.println("  Memory availabe "
				+ Runtime.getRuntime().freeMemory() + "     "
				+ Runtime.getRuntime().totalMemory());
	}


	// Decode a sequence
	/*****************************************************************************/
	public void decode(AbstractSequence[] seqArray, File encodedFile)
			throws IOException {
		FileInputStream fileIn = new FileInputStream(encodedFile);
		int length = (new DataInputStream(fileIn)).readInt();

		// Get the sequence to be encode
		seqArray[seqArray.length - 1] = new Sequence(Alphabet.DNA4(), length);

		initiliseCommon(seqArray);

		int sid = seqArray.length - 1;
		AbstractSequence seq = seqArray[sid];

		ArithDecoder decoder = new ArithDecoder(new BitInput(fileIn));
		
		currentInd = 0;
		while (!decoder.endOfStream()) {
			int mid = decoder.getCurrentSymbolCount(total);
			if (mid >= total)
				break;
			preCoding();

			int actual = 0;
			double accu = 0;

			// Will later implement using binary search, for now just linear
			// search

			while (mid >= (int) ((accu + finalD[actual]) * total)) {
				accu += finalD[actual];
				actual++;
			}

			decoder.removeSymbolFromStream((int) (accu * total),
					(int) ((accu + finalD[actual]) * total), total);
			seq.setSymbol(currentInd,  actual);
	
			updateExperts(actual);
			postCoding(seqArray, sid);

			currentInd++;
			if (currentInd >= seq.length()) {
				break;
			}
		}
		decoder.close();
	}

	/**********************************************************************/
	public File realEncode(AbstractSequence[] seqArray, String filename) {

		try {
			initiliseCommon(seqArray);

			// Get the sequence to be encode
			int sid = seqArray.length - 1;
			AbstractSequence seq = seqArray[seqArray.length - 1];
			
			File file = new File(filename);
			FileOutputStream fileOut = new FileOutputStream(file);

			// Write the length of the sequence, this is to make encoding easier
			int len = seq.length();

			(new DataOutputStream(fileOut)).writeInt(len);

			ArithEncoder encoder = new ArithEncoder(fileOut);
			// initilise(seqHash);
			for (currentInd = 0; currentInd < seq.length(); currentInd++) {
				preCoding();
				int actual = seq.symbolAt(currentInd);

				/**********************************************************************/
				double accu = 0;
				for (int j = 0; j < actual; ++j) {
					accu += finalD[j];
				}
				int low = (int) (accu * total);
				int high = (int) ((accu + finalD[actual]) * total);
				encoder.encode(low, high, total);
				/**********************************************************************/

				updateExperts(actual);
				postCoding(seqArray, sid);
			}
			encoder.close();

			return file;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	/**********************************************************************/
	public double[] encode(AbstractSequence[] seqArray) {

		initiliseCommon(seqArray);

		// Get the sequence to be encode
		AbstractSequence seq = seqArray[seqArray.length - 1];
		double[] costs = new double[seq.length()];
		// double [] markovCosts = new double [japsa.seq.length];
		int sid = seqArray.length - 1;

		double totalCost = 0.0;

		// go thru the sequence and encode each charactor
		for (currentInd = 0; currentInd < seq.length(); currentInd++) {
			// Compute the probability distribution
			preCoding();

			int actual = seq.symbolAt(currentInd);
			double cost = -JapsaMath.log2(finalD[actual]);

			costs[currentInd] = cost;
			totalCost += cost;


			updateExperts(seq.symbolAt(currentInd));
			postCoding(seqArray, sid);

			/*******************************************************************/
			// #CHECKPOINT_BEGIN
			if ((currentInd + 1) % checkPoint == 0) {//
				checkPoint(currentInd + 1);
				System.out.println("    Current comp = "
						+ (totalCost / currentInd) + "   ("
						+ new java.util.Date() + ")");

			}
			// #CHECKPOINT_END
			/*******************************************************************/
		}
		return costs;
	}

	public double encode1(AbstractSequence[] seqArray) {
		initiliseCommon(seqArray);

		// Get the sequence to be encodes		
		int sid = seqArray.length - 1;
		AbstractSequence seq = seqArray[seqArray.length - 1]; 

		double totalCost = 0.0;

		// go thru the sequence and encode each charactor
		for (currentInd = 0; currentInd < seq.length(); currentInd++) {
			// Compute the probability distribution
			preCoding();

			int actual = seq.symbolAt(currentInd);
			// double cost = -JapsaMath.log2(finalDist.getWeight(actual));
			double cost = -JapsaMath.log2(finalD[actual]);

			totalCost += cost;

			updateExperts(actual);
			postCoding(seqArray, sid);

			/*******************************************************************/
			// #CHECKPOINT_BEGIN
			if ((currentInd + 1) % checkPoint == 0) {//
				checkPoint(currentInd + 1);
				System.out.println("    Current comp = "
						+ (totalCost / currentInd) + "   ("
						+ new java.util.Date() + ")");

			}
			// #CHECKPOINT_END
			/*******************************************************************/
		}

		return totalCost / seq.length();
	}

	public double encode_optimise(AbstractSequence[] seqArray) throws Exception {

		this.initilise_optimise(seqArray);
		// initiliseCommon(seqArray);

		// Get the sequence to be encode		
		int sid = seqArray.length - 1;
		AbstractSequence seq = seqArray[seqArray.length - 1]; 

		double totalCost = 0.0;

		// go thru the sequence and encode each charactor
		for (currentInd = 0; currentInd < seq.length(); currentInd++) {
			// Compute the probability distribution
			preCoding();

			int actual =  seq.symbolAt(currentInd);
			totalCost -= JapsaMath.log2(finalD[actual]);

			updateExperts(actual);
			postCoding(seqArray, sid);

			/*******************************************************************/
			// #CHECKPOINT_BEGIN
			if ((currentInd + 1) % checkPoint == 0) {//
				checkPoint(currentInd + 1);
				System.out.println("    Current comp = "
						+ (totalCost / currentInd) + "   ("
						+ new java.util.Date() + ")");

			}
			// #CHECKPOINT_END
			/*******************************************************************/
		}

		return (totalCost / seq.length());
	}

	public void encode(AbstractSequence[] seqArray, String infoFile, String markovFile)
			throws Exception {

		initiliseCommon(seqArray);

		// Get the sequence to be encode
		int sid = seqArray.length - 1;
		AbstractSequence seq = seqArray[seqArray.length - 1]; 

		double totalCost = 0.0, totalMarkovCost = 0.0, cost;
		PrintStream infoPs = new PrintStream(new BufferedOutputStream(new FileOutputStream(infoFile)));
		PrintStream markovPs = null;
		
		if (markovFile != null) {
			markovPs = new PrintStream(new BufferedOutputStream(new FileOutputStream(markovFile)));
			markovPs.println("#Information content produced using Markov Model by the eXpert Model (XM,DCC'07, doi:10.1109/DCC.2007.7) ");
		}

		infoPs.println("#Information content produced by the eXpert Model (XM,DCC'07, doi:10.1109/DCC.2007.7) ");

		// go thru the sequence and encode each charactor
		for (currentInd = 0; currentInd < seq.length(); currentInd++) {
			// Compute the probability distribution
			preCoding();

			int actual = seq.symbolAt(currentInd);
			// double cost = -JapsaMath.log2(finalDist.getWeight(actual));
			cost = -JapsaMath.log2(finalD[actual]);
			totalCost += cost;
			infoPs.println(cost);

			// For markov cost
			if (markovPs != null) {
				cost = -JapsaMath.log2(markovD[actual]);
				totalMarkovCost += cost;
				markovPs.println(cost);
			}

			if (cost <= 0) {
				System.err.println(cost + "\t" + finalD[actual]);
			}

			updateExperts(actual);
			postCoding(seqArray, sid);

			/*******************************************************************/
			// #CHECKPOINT_BEGIN
			if ((currentInd + 1) % checkPoint == 0) {//
				checkPoint(currentInd + 1);
				System.out.println("    Current comp = "
						+ (totalCost / currentInd) + "   ("
						+ new java.util.Date() + ")");

			}
			// #CHECKPOINT_END
			/*******************************************************************/
		}
		infoPs.close();
		if (markovPs != null)
			markovPs.close();

		System.out.println(totalCost / seq.length() + "    " + totalMarkovCost
				/ seq.length());
	}

	/******************************************************************************/

	public static String version() {
		return "The eXpert-Model (XM) for Compression of DNA Sequences "
				+ VERSION
				+ "\n  Minh Duc Cao, T. I. Dix, L. Allison, C. Mears."
				+ "\n  A simple statistical algorithm for biological sequence compression"
				+ "\n  IEEE Data Compression Conf., Snowbird, Utah, 2007, [doi:10.1109/DCC.2007.7]\n";
	}

	public static void main(String[] args) throws Exception {
		ExpertModelCmd.main(args);
	}

}
