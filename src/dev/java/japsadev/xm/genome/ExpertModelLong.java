/*****************************************************************************
 * Copyright (c) 2010 Minh Duc Cao, Monash University.  All rights reserved. *
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
 * 3. Neither the name of Monash University nor the names of its contributors*
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

package japsadev.xm.genome;

import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.JapsaMath;
import japsadev.xm.genome.MyBigHashtableLong.MyBigHashLongIterator;

import java.io.*;

import com.colloquial.arithcode.*;

public class ExpertModelLong {
	public boolean optimise = false;

	public static String VERSION = "V 2.0";

	// #DEBUG_BEGIN
	static public int DEBUG_LEVEL = 2;
	// #DEBUG_END

	static final long TWO_MAX_INT = 1l << 32;

	// #CHECKPOINT_BEGIN
	protected int checkPoint = 1000000;
	// #CHECKPOINT_END

	// 'global' Variables
	protected int expertCount;// The number of experts currently empliyed
	protected long[] accLengths;// lengths of sequences (including context)

	// protected CombinationExpert baseEx;//Any fixed expert are in this linked
	// list
	MarkovExpertLong markovEx;
	AdaptiveMarkovExpertLong adapMarkovEx;//
	// = new MarkovExpert(null,2);
	double[] finalD;// ,baseD;
	double[] markovD;
	public String hashType = "hash";
	long currentInd;

	public static double[] markovCost = null;

	protected CombinationExpertLong repEx;
	protected RepeatExpertLong restHead = null;

	protected double repeatPrior;
	protected double repeatPriorProb;

	// protected Distribution finalDist;
	// protected int posSize = 0;
	protected long lengthSeqs = 0;// Total lengths of all sequences
	protected int total = 1024 * 1024;
	MyBitSetLong bitSet, pBitSet;// To store if an particular offset has been
									// stored

	// Two expert seeds
	// public RepeatExpertLong offSetSeed; //= new
	// OffsetCountExpert(null,1,null);// a dummy one
	// public RepeatExpertLong palinSeed;// = new
	// PalindromeCountExpert(null,1,null);// a dummy one
	/*
	 * Parameters of the algorithm
	 */
	protected int expertsLimit;// Limit of experts
	protected int chances;// Number of chances given to each expert
	// PatternStore myHash; //The hashtable

	// HASH MyHashtableLong myHash; //The hashtable

	MyBigHashtableLong myHash; // The hashtable

	int hashSize;

	boolean selfRep = false;
	boolean binaryHash = false;

	public ExpertModelLong(int hashSize, int alphabetSize, int context,
			int expertsLimit, double listenThreshold, int chances) {
		super();

		// this.hashSize = hashSize;
		ExpertLong.setParams(alphabetSize, context);
		this.expertsLimit = expertsLimit;
		this.chances = chances;
		this.repeatPrior = listenThreshold * context;
		repeatPriorProb = JapsaMath.exp2(-repeatPrior);

		// thExp = JapsaMath.exp2(listenThreshold);
		// this.listenThreshold = Math.pow(2,listenThreshold);

		// finalDist = new Distribution(Expert.ALPHABET_SIZE);

		finalD = new double[ExpertLong.ALPHABET_SIZE];
		markovD = new double[ExpertLong.ALPHABET_SIZE];// Combination of the
														// markov experts
		this.hashSize = hashSize;
		// this.binaryHash = binaryHash;
	}

	public void setBinaryHash(boolean useBinary) {
		binaryHash = useBinary;
	}

	public void setSelfRep(boolean rep) {
		selfRep = rep;
	}

	public void setHashType(String hashT) throws Exception {
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
					throw new Exception("Unknown hash type : " + hashT);
				}
			} // assert: hashT ok

			int i = hashT.length() - 1;
			while (i >= 0 && hashT.charAt(i) == '0')
				i--;
			if (i < 0)
				throw new Exception("Unknown hash type : " + hashT);
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
				+ "\nContext          : " + ExpertLong.CONTEXT_LENGTH
				+ "\nListen Threshold : " + repeatPrior
				/ ExpertLong.CONTEXT_LENGTH + "bps" + "\nChances          : "
				+ chances + "\nBinaryHash       : " + binaryHash
				+ "\nHashType         : " + hashName + "\nExpert Type      : "
				+ restHead.getClass());
	}

	// protected abstract void initilise(BioCompSequence[] seqArray);//should
	// store
	/*************************************************************************
	 * protected void initiliseGappedHash(BioCompSequence[] seqArray){ myHash =
	 * new GappedHashtable(hashType); store(seqArray); } /
	 *************************************************************************/
	protected void initiliseHash(GenomeSequence[] seqArray) {

		// if (binaryHash){
		// myHash = new MyBinaryHashtable(hashSize);
		// }else

		// HASH myHash = new MyHashtableLong(hashSize);
		myHash = new MyBigHashtableLong(seqArray.length, hashSize);
		store(seqArray);
	}

	/*************************************************************************
	 * protected void initilisePrefixArray(BioCompSequence[] seqArray){
	 * 
	 * //#TIME_BEGIN long buildingTime=System.currentTimeMillis();
	 * System.out.println(" #Building Suffix Array"); //#TIME_END
	 * 
	 * PrefixArrayAbstract sArray;
	 * 
	 * if (binaryHash) sArray = new PrefixArrayBinary(seqArray[0].toBytes());
	 * else sArray = new PrefixArray(seqArray[0].toBytes());
	 * 
	 * myHash = sArray; sArray.setHashSize(hashSize);
	 * 
	 * sArray.setSeqToMatch(seqArray[1].toBytes()); //#TIME_BEGIN
	 * System.out.printf
	 * (" #Suffix Array (%s) built in %d ms\n",sArray.getClass()
	 * ,(System.currentTimeMillis() - buildingTime)); //#TIME_END }
	 * /************
	 * *************************************************************
	 * 
	 * 
	 * private void initilisePrefixTree(BioCompSequence[] seqArray){ //myHash =
	 * new SuffixTree(seqArray[0].toBytes()); SuffixTree tree = new
	 * SuffixTree(seqArray[0].toBytes()); tree.setMinHash(hashSize); myHash =
	 * tree; tree.setIncrementLeaf(1);
	 * 
	 * store(seqArray); System.out.println("Tree build done");
	 * tree.setIncrementLeaf(0); tree.setRefSeq(seqArray[seqArray.length -
	 * 1].toBytes()); }
	 * 
	 * /
	 *************************************************************************/

	public void store(GenomeSequence[] seqArray) {
		// Store all back ground sequence in the hash
		for (int sid = 0; sid < seqArray.length - 1; sid++) {
			for (long i = 0; i < seqArray[sid].getLength(); i++) {
				myHash.nextKey(seqArray[sid].getBase(i));
				// HASH myHash.putCurrentValue((int)((sid << posSize) + i));
				myHash.putCurrentValue(sid, (int) (i));
			}
		}
	}

	protected void initilise_optimise(GenomeSequence[] seqArray) {
		// Common stuff
		lengthSeqs = 0;

		// Compute accumulated length
		accLengths = new long[seqArray.length + 1];
		accLengths[0] = 0;
		for (int i = 0; i < seqArray.length; i++) {
			lengthSeqs += seqArray[i].getLength();
			accLengths[i + 1] = lengthSeqs;
		}

		expertCount = 0;

		// Background knowledge/based knowledge
		markovEx = new MarkovExpertLong(null, 2);
		adapMarkovEx = new AdaptiveMarkovExpertLong(null, 1);
		markovEx.setNext(adapMarkovEx);

		repEx = new CombinationExpertLong();

		// HASH MyHashtableLong hash = new MyHashtableLong(hashSize);
		MyBigHashtableLong hash = new MyBigHashtableLong(seqArray.length,
				hashSize);
		/*************************************************************************/

		System.out.println("Run first pass");
		// In this pass, count the needed cell on the hash table

		for (int sid = 0; sid < seqArray.length; sid++) {
			for (long i = 0; i < seqArray[sid].getLength(); i++) {
				hash.nextKey(seqArray[sid].getBase(i));
				// HASH hash.putCurrentValue_psuedo((int) ((sid << posSize) +
				// i));
				hash.putCurrentValue_psuedo(sid, (int) (i));
			}
		}
		hash.printMemoryNeeded();

		hash.reinitialise_optimise();
		System.out.println("Finish first pass");

		myHash = hash;
		store(seqArray);

		checkPoint(0);

		bitSet = new MyBitSetLong(lengthSeqs);
		pBitSet = new MyBitSetLong(lengthSeqs * 2);

	}

	/*************************************************************************/
	protected void initiliseCommon(GenomeSequence[] seqArray) {
		// Common stuff
		lengthSeqs = 0;

		// Compute accumulated length
		accLengths = new long[seqArray.length + 1];
		accLengths[0] = 0;
		for (int i = 0; i < seqArray.length; i++) {
			lengthSeqs += seqArray[i].getLength();
			accLengths[i + 1] = lengthSeqs;
		}

		expertCount = 0;

		// Background knowledge/based knowledge
		markovEx = new MarkovExpertLong(null, 2);
		adapMarkovEx = new AdaptiveMarkovExpertLong(null, 1);
		markovEx.setNext(adapMarkovEx);

		repEx = new CombinationExpertLong();

		// MyHashtable hash = new MyHashtable(hashSize);
		/*************************************************************************/
		// Initilise hashtable
		initiliseHash(seqArray);

		// if ("sft".equals(this.hashType))
		// initilisePrefixTree(seqArray);
		// else if ("sfa".equals(this.hashType))
		// initilisePrefixArray(seqArray);
		// else if ("hash".equals(this.hashType))
		// initiliseHash(seqArray);
		// else{//gappped
		// initiliseGappedHash(seqArray);
		// }

		bitSet = new MyBitSetLong(lengthSeqs);
		pBitSet = new MyBitSetLong(lengthSeqs * 2);
	}

	/*****************************************************************/
	protected void preCoding() {
		double baseRateProb = markovEx.rateProb();

		for (byte a = 0; a < ExpertLong.ALPHABET_SIZE; a++) {
			markovD[a] = finalD[a] = markovEx.rateProb()
					* markovEx.probability(a) + adapMarkovEx.rateProb()
					* adapMarkovEx.probability(a);
		}

		// Initialise expert prediction and total prediction
		this.repEx.getCombDistribution().setWeights(0.0);

		double repSum = 0.0;

		// Listen to no more than expert limit experts
		ExpertLong ptr = repEx.getNext();// go thro all rep experts
		ExpertLong head = repEx;
		while (ptr != null) {// Inv: head.next = ptr
			// double rate = ptr.rate();//msg leng over an history
			double score = ptr.rateProb(); // costToWeight(rate);
			// Only listen if the offset/palindrome is sufficiently good
			// if (rate < baseRate - repeatPrior){// * Expert.CONTEXT_LENGTH){
			if (score > baseRateProb / repeatPriorProb) {// *
															// Expert.CONTEXT_LENGTH){
				for (byte a = 0; a < ExpertLong.ALPHABET_SIZE; a++) {
					repEx.getCombDistribution().addWeight(a,
							score * ptr.probability(a));
				}
				// ptr.weight = myCost;
				repSum += score;
				ptr.resetCounter();
				// goodRep.add(ptr);
			} else {
				// ptr.weight = 0.0;
				ptr.incrementCounter();

				if (ptr.getCounter() > chances) {
					// remove ptr
					resignExpert(ptr);// ptr.resign();

					expertCount--;

					head.setNext(ptr.getNext());
					ExpertLong retired = ptr;

					ptr = ptr.getNext();
					retired.next = restHead;
					restHead = (RepeatExpertLong) retired;

					continue;
				}
			}
			head = ptr;
			ptr = ptr.getNext();
		}

		if (repSum != 0.0) {// There is repeat
			repEx.incrementCount();
			repEx.getCombDistribution().scale(1.0 / repSum);

			// double repExRate = repEx.rate() +
			// repeatPrior;//JapsaMath.log2(repEx.getPrior()) ;
			// double repExScore = costToWeight(repExRate);

			double repExScore = repEx.rateProb() * repeatPriorProb;

			// Combine 2 experts
			for (byte a = 0; a < ExpertLong.ALPHABET_SIZE; a++) {
				// finalDist.addWeight(a,
				// repEx.getCombDistribution().getWeight(a) * repExScore);
				finalD[a] += repEx.getCombDistribution().getWeight(a)
						* repExScore;
			}
			// finalSum += repExScore;
		} else {
			repEx.getCombDistribution().setWeights(
					1.0 / ExpertLong.ALPHABET_SIZE);
		}

		double fSum = 0;
		for (byte a = 0; a < ExpertLong.ALPHABET_SIZE; a++) {
			// finalDist.addWeight(a, repEx.getCombDistribution().getWeight(a) *
			// repExScore);
			fSum += finalD[a];
		}

		for (byte a = 0; a < ExpertLong.ALPHABET_SIZE; a++) {
			// finalDist.addWeight(a, repEx.getCombDistribution().getWeight(a) *
			// repExScore);
			finalD[a] /= fSum;
		}
	}

	protected void resurrectExpert(RepeatExpertLong e, GenomeSequence bs,
			long pos, int past) {
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
		// Update base experts
		// Expert

		markovEx.update(c);
		adapMarkovEx.update(c);

		// Update rep experts
		// baseEx.update(c);
		repEx.update(c);
		/********************* Update Exp ********************/
		ExpertLong ptr = repEx.getNext();
		ExpertLong head = repEx;

		ExpertLong retired;

		while (ptr != null) {
			double thisCost = ptr.update(c);// This actually the prob
			if (thisCost > 0) {// A cost should be a positive number
				head = ptr;
				ptr = ptr.getNext();
			} else {// resign
				resignExpert(ptr);// ptr.resign();
				expertCount--;
				head.setNext(ptr.getNext());
				retired = ptr;
				ptr = ptr.getNext();

				// put retired into the rest list
				retired.next = this.restHead;
				this.restHead = (RepeatExpertLong) retired;
			}
		}
	}

	void resignExpert(ExpertLong p) {
		p.resign();
	}

	// static int id_zero = 0;
	protected void postCoding(GenomeSequence[] seqArray, int sid) {
		/**********************************************************************/
		// Move to next key
		// if (currentInd < seqArray[sid].length() - 1)
		myHash.nextKey(seqArray[sid].getBase(currentInd));

		// Add good diagonals.
		if (currentInd >= hashSize
				&& currentInd < seqArray[sid].getLength() - 1) {
			// assert n >= 0
			// IntIterator iter = myHash.copyIterator();
			// HASH MyHashLongIterator iter = myHash.getLongIterator();
			MyBigHashLongIterator iter = myHash.iterator();

			// #DEBUG_BEGIN
			if (DEBUG_LEVEL > 3) {
				System.out.printf("Find : %10d%10d\n", iter.sizeAvailable(),
						expertCount);
			}
			// #DEBUG_END

			while (expertsLimit > expertCount && iter.hasNext()) {
				int position = iter.next();
				RepeatExpertLong e = null;
				int id = iter.sid;

				if (!iter.isPalin) {// Offset expert
					long pos = position;
					if (pos < 0) {
						pos = TWO_MAX_INT + pos;
					}

					if (pos <= hashSize) {
						continue;
					}

					if (pos > hashSize && // Have enough for resurrect
							!bitSet.get(currentInd + accLengths[sid]
									- accLengths[id] - pos)// Not in there
							&& pos < seqArray[id].getLength() - 3) {// have some
																	// thing to
																	// predict
						// get the next free
						e = this.restHead;
						this.restHead = (RepeatExpertLong) this.restHead.next;
						e.reuseExpert(seqArray[id], pos, bitSet,
								RepeatExpertLong.COPY_TYPE);
						// e = offSetSeed.duplicate(seqArray[id], posSrc, bitSet);
						e.setID(currentInd + accLengths[sid] - accLengths[id]
								- pos);
					}
				} else {// Palindrome expert
					long pos = position;
					if (pos < 0) {
						pos = TWO_MAX_INT + pos;
					}

					if (pos > hashSize
							&& (!pBitSet.get(currentInd + accLengths[sid]
									+ accLengths[id] + pos))// ?
							&& (pos + 3 < seqArray[id].getLength())) {// Have
																		// something
																		// to
																		// predict
						e = this.restHead;
						this.restHead = (RepeatExpertLong) this.restHead.next;
						e.reuseExpert(seqArray[id], pos - hashSize + 1,
								pBitSet, RepeatExpertLong.PALIN_TYPE);
						// e = palinSeed.duplicate(seqArray[id], posSrc - hashSize
						// + 1,pBitSet);
						e.setID(currentInd + accLengths[sid] + accLengths[id]
								+ pos);
					}
				}
				// Add this expert in only if an identical expert not in the
				// list
				if (e != null) {
					resurrectExpert(e, seqArray[sid], currentInd, hashSize);
					// e.resurrect(seqArray[sid].toBytes(),currentInd,hashSize);
					e.setNext(repEx.getNext());
					repEx.setNext(e);
					expertCount++;
				}
			}
		}
		// System.out.println(i + "   " + expertCount);
		/********************** Store current *********************/
		if (selfRep)
			// HASH myHash.putCurrentValue( (int) (currentInd ));
			myHash.putCurrentValue(sid, (int) (currentInd));
		// myHash.putCurrentValue( (int) ((sid << posSize) + currentInd));
		/******************************************************************/
	}

	protected void checkPoint(long steps) {
		System.out.print("Reach milestone " + steps + " : ");
		myHash.printSummary();

		System.out.println("  Memory availabe "
				+ Runtime.getRuntime().freeMemory() + "     "
				+ Runtime.getRuntime().totalMemory());
		// System.gc();
		System.out.println("  Memory availabe "
				+ Runtime.getRuntime().freeMemory() + "     "
				+ Runtime.getRuntime().totalMemory());
	}

	// public static double costToWeight(double cost){
	// return JapsaMath.exp2(-cost);
	// }

	// Decode a sequence
	/*****************************************************************************/
	public void decode(GenomeSequence[] seqArray, File encodedFile)
			throws IOException {
		FileInputStream fileIn = new FileInputStream(encodedFile);
		long length = (new DataInputStream(fileIn)).readLong();

		// Get the sequence to be encode
		int sid = seqArray.length - 1;
		GenomeSequence seq = new GenomeSequence(length);
		seqArray[seqArray.length - 1] = seq;

		initiliseCommon(seqArray);

		ArithDecoder decoder = new ArithDecoder(new BitInput(fileIn));

		currentInd = 0;
		double accu = 0;
		while (!decoder.endOfStream()) {
			int mid = decoder.getCurrentSymbolCount(total);
			if (mid >= total)
				break;
			preCoding();

			int actual = 0;
			// double
			accu = 0;

			// Will later implement using binary search, for now just linear
			// search
			try {
				while (mid >= (int) ((accu + finalD[actual]) * total)) {
					accu += finalD[actual];
					actual++;
				}

				decoder.removeSymbolFromStream((int) (accu * total),
						(int) ((accu + finalD[actual]) * total), total);
				seq.putBase(currentInd, actual);
			} catch (Exception e) {
				e.printStackTrace();
				System.out.println(actual + "  " + mid + "  " + total + "  at "
						+ currentInd + "  acc = " + accu);
				double nAcc = 0;
				for (int y = 0; y < finalD.length; y++) {
					nAcc += finalD[y];
					System.out.println(finalD[y] + "   " + nAcc + " "
							+ (nAcc * total) + "   " + ((int) (nAcc * total)));
				}
				System.exit(-1);
			}
			// if (seqDec[i] != japsa.seq[i]){
			// System.out.println("Wrong at the possition " + i + " expected " +
			// japsa.seq[i] + " but see " + actual);
			// System.out.println(mid + "   " + accu * total + "  " + (accu +
			// final_distribution[actual]) * total);
			// return false;
			// }

			updateExperts(actual);

			postCoding(seqArray, sid);

			currentInd++;
			if (currentInd >= seq.getLength()) {
				break;
				// return seqDec;
			}
			/*******************************************************************/
			// #CHECKPOINT_BEGIN
			if ((currentInd + 1) % checkPoint == 0) {//
				checkPoint(currentInd + 1);
				System.out.println("    Current decoding at ("
						+ new java.util.Date() + ")");
			}
			// #CHECKPOINT_END
			/*******************************************************************/
		}// while
		decoder.close();
	}

	/*******************************************************************/
	public double encode1(GenomeSequence[] seqArray) {
		long s1 = System.currentTimeMillis();
		if (optimise)
			this.initilise_optimise(seqArray);
		else
			initiliseCommon(seqArray);

		// initiliseCommon(seqArray);

		long s2 = System.currentTimeMillis();
		System.out.println((s2 - s1));

		// Get the sequence to be encode
		int sid = seqArray.length - 1;
		GenomeSequence seq = seqArray[sid];// seqArray.length -1];

		double totalCost = 0.0;

		// go thru the sequence and encode each charactor
		for (currentInd = 0; currentInd < seq.getLength(); currentInd++) {
			// Compute the probability distribution
			preCoding();

			int actual = seq.getBase(currentInd);
			double cost = -JapsaMath.log2(finalD[actual]);

			totalCost += cost;

			updateExperts(seqArray[sid].getBase(currentInd));
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
		long s3 = System.currentTimeMillis();
		System.out.println((s3 - s2));

		return totalCost / seq.getLength();
	}

	public File realEncode(GenomeSequence[] seqArray, String filename) {

		try {
			if (optimise)
				this.initilise_optimise(seqArray);
			else
				initiliseCommon(seqArray);

			// Get the sequence to be encode
			int sid = seqArray.length - 1;
			GenomeSequence seq = seqArray[sid];
			double totalCost = 0.0;

			File file = new File(filename);
			FileOutputStream fileOut = new FileOutputStream(file);

			// Write the length of the sequence, this is to make encoding easier
			long len = seq.getLength();

			(new DataOutputStream(fileOut)).writeLong(len);

			ArithEncoder encoder = new ArithEncoder(fileOut);
			// initilise(seqHash);
			for (currentInd = 0; currentInd < seq.getLength(); currentInd++) {
				preCoding();
				int actual = seq.getBase(currentInd);

				/**********************************************************************/
				double accu = 0;
				for (int j = 0; j < actual; ++j) {
					accu += finalD[j];
				}
				int low = (int) (accu * total);
				int high = (int) ((accu + finalD[actual]) * total);
				// System.out.println("Encode " + low + "  " + high + "  " +
				// total);
				encoder.encode(low, high, total);
				/**********************************************************************/

				updateExperts(actual);
				postCoding(seqArray, sid);

				totalCost -= JapsaMath.log2(finalD[actual]);
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

			System.out.println(" Theoritical bps "
					+ (totalCost / seq.getLength()));
			encoder.close();

			return file;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	public static String version() {
		return "The eXpert-Model (XM) for Compression of DNA Sequences "
				+ VERSION
				+ "\n   Minh Duc Cao, T. I. Dix, L. Allison, C. Mears."
				+ "\n   A simple statistical algorithm for biological sequence compression"
				+ "\n   IEEE Data Compression Conf., Snowbird, Utah, 2007, [doi:10.1109/DCC.2007.7]\n";
	}

	public static void main(String[] args) {
		try {
			// Get params from users
			CommandLine cmdLine = new CommandLine();
			cmdLine.addInt("hashSize", 11, "Hash size");
			cmdLine.addInt("context", 15, "Length of the context");
			cmdLine.addInt("limit", 200, "Expert Limit");
			cmdLine.addDouble("threshold", 0.15, "Listen threshold");
			cmdLine.addInt("chance", 20, "Chances");
			cmdLine.addString(
					"offsetType",
					"counts",
					"Way of update offset/palindrome expert: possible value count, subs, viaprotein,static");
			cmdLine.addBoolean(
					"optimise",
					false,
					"Running in optimise mode, just report the entropy,recommended for long sequence");

			cmdLine.addString("real", null, "File name of the real compression");
			cmdLine.addString("decode", null, "File name of the encoded");
			cmdLine.addString("output", "decoded",
					"The output file of decoded file");

			// #CHECKPOINT_BEGIN
			cmdLine.addInt("checkPoint", 1000000, "Frequency of check point");
			// #CHECKPOINT_END

			args = cmdLine.parseLine(args);

			int alphabetSize = 4;
			int hashSize = cmdLine.getIntVal("hashSize");
			int context = cmdLine.getIntVal("context");
			int expertsLimit = cmdLine.getIntVal("limit");
			double listenThreshold = cmdLine.getDoubleVal("threshold");
			int chances = cmdLine.getIntVal("chance");

			ExpertModelLong eModel = new ExpertModelLong(hashSize,
					alphabetSize, context, expertsLimit, listenThreshold,
					chances);

			eModel.optimise = cmdLine.getBooleanVal("optimise");

			eModel.setHashType("hash");
			eModel.setSelfRep(true);

			// eModel.offSetSeed = new RepeatCountExpertLong(new
			// GenomeSequence(1),0,null,RepeatExpertLong.COPY_TYPE);
			// eModel.palinSeed = new RepeatCountExpertLong(new
			// GenomeSequence(1),0,null,RepeatExpertLong.PALIN_TYPE);

			eModel.restHead = new RepeatCountExpertLong(new GenomeSequence(1),
					0, null, RepeatExpertLong.COPY_TYPE);
			ExpertLong tail = eModel.restHead;

			for (int x = 0; x < expertsLimit + 1; x++) {
				RepeatCountExpertLong t = new RepeatCountExpertLong(
						new GenomeSequence(1), 0, null,
						RepeatExpertLong.COPY_TYPE);
				tail.setNext(t);
				tail = t;
			}

			System.out.println(ExpertModelLong.version());
			// Print out all params
			eModel.printParams();

			if (cmdLine.getStringVal("decode") == null
					&& (args == null || args.length <= 0)) {
				System.err.println("Usage: java CommandLine [options]\n"
						+ cmdLine.usageMessage() + " file1 file2 ...");
				System.exit(1);
			}
			// cmdLine.printOptions();

			// #TIME_BEGIN
			long timeStart = System.currentTimeMillis();
			System.out.println(" #Reading file(s)");
			// #TIME_END

			// BioCompDNA[] dnaArray;//
			GenomeSequence[] genSeqs;
			if (cmdLine.getStringVal("decode") != null) {
				if (args == null)
					args = new String[0];

				genSeqs = new GenomeSequence[args.length + 1];
			} else {// Read per normal
				genSeqs = new GenomeSequence[args.length];
			}
			for (int i = 0; i < args.length; i++) {
				System.out.print("Read " + args[i] + "...");
				genSeqs[i] = GenomeSequence.guessFormat(args[i]);
				System.out.println(" done (" + genSeqs[i].getLength() + ")");
			}

			// #TIME_BEGIN
			long timeEnd = System.currentTimeMillis();
			System.out.println(" #Read file(s) in " + (timeEnd - timeStart)
					+ "ms ");
			// #TIME_END

			// if (cmdLine.getStringVal("hashType").equals("sft")){
			// Augment all the background sequences
			// GenomeSequence[] augdnaArray = new GenomeSequence[2];
			// augdnaArray[1] = genSeqs[genSeqs.length - 1];
			// TODO: fix the following
			// augdnaArray[0] = new BioCompDNA(new byte[0],"Combine");

			// for (int x = 0; x < dnaArray.length - 1; x++){
			// augdnaArray[0].concatenate(dnaArray[x]);
			// augdnaArray[0].concatenate(dnaArray[x].reverseComplement());
			// }
			// dnaArray = augdnaArray;
			// }else if (cmdLine.getStringVal("hashType").equals("sfa") &&
			// genSeqs.length > 1){
			// Augment all the background sequences for suffix array
			// BioCompDNA[] augdnaArray = new BioCompDNA[2];
			// augdnaArray[1] = dnaArray[dnaArray.length - 1];

			// int combileLength = 0;
			// for (int x = 0; x < dnaArray.length - 1; x++){
			// combileLength += dnaArray[x].length();
			// }

			// Combine all sequences, including the rev comp into 1/
			// in the reverse order. the correct order is retrieved later
			// byte [] combineByte = new byte[combileLength * 2];///Backward and
			// forward
			// int ind = 0;

			// for (int x = 0; x < dnaArray.length - 1; x++){
			// byte[] seqX = dnaArray[x].toBytes();
			// for (int y = 0; y < seqX.length; y++){
			// //
			// combineByte[ind ++] = seqX[y];//(byte)(3 - seqX[y]);//compliment
			// seqX[y];
			// combineByte[combineByte.length - ind] = (byte)(3 -
			// seqX[y]);//seqX[y];//(byte)(3 - seqX[y]);//compliment
			// }
			// }

			// augdnaArray[0] = new BioCompDNA(combineByte,"Combine");
			// dnaArray = augdnaArray;
			// }

			// for (int i = 0; i < dnaArray.length; i++){
			// if (i == dnaArray.length - 1){
			// if (cmdLine.getStringVal("decode") == null){
			// System.out.println("Encode  : " + dnaArray[i] +
			// "\t"+dnaArray[i].length()+"");
			// }
			// }else
			// System.out.println("Context : " + dnaArray[i] +
			// "\t"+dnaArray[i].length()+"");
			// }
			System.out
					.println("----------------------------------------------------------------------");

			long start;
			/*************************************************************************
			 * if (cmdLine.getBooleanVal("optimise")){ //System.gc(); start =
			 * System.currentTimeMillis();
			 * 
			 * double cost = eModel.encode_optimise(genSeqs); long time =
			 * (System.currentTimeMillis() - start); //
			 * System.out.printf("%f bps in %d ms\n",cost,time);
			 * System.out.println(
			 * "======================================================================"
			 * );
			 * 
			 * }else /
			 *************************************************************************/
			if (cmdLine.getStringVal("decode") != null) {
				if (cmdLine.getBooleanVal("optimise")) {
					System.err
							.println("Warn: optimise option cannot be used with decode, disabled");
				}

				String file = cmdLine.getStringVal("decode");
				String output = cmdLine.getStringVal("output");

				System.out.println("Decoding");
				start = System.currentTimeMillis();

				eModel.decode(genSeqs, new File(file));

				System.out.println(" Time decode "
						+ (System.currentTimeMillis() - start) / 1000.0
						+ "seconds\n");
				SequenceOutputStream outPrintStream = SequenceOutputStream
						.makeOutputStream(output);
				genSeqs[genSeqs.length - 1].write(outPrintStream);
				outPrintStream.close();

			}
			/*************************************************************************/
			else if (cmdLine.getStringVal("real") != null) {
				// System.gc();
				System.out.println("Real  encoding");
				start = System.currentTimeMillis();
				File outputFile = eModel.realEncode(genSeqs,
						cmdLine.getStringVal("real"));

				System.out.println(" Time encode "
						+ (System.currentTimeMillis() - start) / 1000.0
						+ "seconds");
				System.out.println(" Encoding cost = "
						+ (outputFile.length() * 8.0)
						/ genSeqs[genSeqs.length - 1].getLength() + "bps");

			} else {// Normal
				// System.gc();
				start = System.currentTimeMillis();
				// seqHash[1] = japsa.seq;
				double costs = eModel.encode1(genSeqs);// ,args[1]);

				long time = (System.currentTimeMillis() - start);
				// System.out.printf(" Comp rate : #%f#\n",total /
				// costs.length);
				System.out.printf("%f bps in %d ms\n", costs, time);
				System.out
						.println("=============================================================================");
			}

		} catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}

}
