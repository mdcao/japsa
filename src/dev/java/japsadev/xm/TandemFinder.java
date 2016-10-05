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
 * 15/03/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.xm;

/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;

import java.util.*;

/**
 * A Tool for finding tandem repeats.
 * 
 * @author Jackson Gatenby (jackson.gatenby@gmail.com)
 * 
 *         Usage:
 * 
 *         <pre>
 * java bio.xm.TandemFinder [reads.fq]
 * </pre>
 * 
 */
public class TandemFinder {

	static ExpertModelTandem tModel;

	// Thresholds for determining if a tandem repeat exists
	public static double REPEAT_THRESHOLD = 1.0;

	/**
	 * Determine if a sequence has a tandem repeat in it.
	 * 
	 * @param seqArray
	 *            the sequence
	 * @return an empty set if no repeat, else a set of the possible lengths of
	 *         the repeat unit.
	 */
	public static Set<Integer> classify(Sequence seqArray) {
		/*
		 * - Do the encoding. - If one of the tandem costs has a low region, it
		 * wins.
		 */
		tModel.calculateTandemCosts(seqArray);

		// In case a count of each successful candidate is needed, use a map.
		// Map<Integer, Integer> candidates = new HashMap<Integer, Integer>();
		Set<Integer> candidates = new HashSet<Integer>();

		/****************************************************************
		 * for (int i = 1; i < seqArray.length(); ++i) { for (int t = 0; t <
		 * tModel.tandemCosts.length; ++t) { if (tModel.tandemCosts[t][i] <
		 * REPEAT_THRESHOLD) { // We're sure about the sequence repeating.
		 * candidates.add(t + 2); } } } /
		 ****************************************************************/

		/*
		 * Debug: *-/ for (int n : candidates) { System.out.print(n + ", "); }
		 * System.out.println(); //
		 */

		// If there is a repeat, multiples of the repeat period will occur
		return candidates;
	}

	@SuppressWarnings("unused")
	private static int mostCommon(Map<Integer, Integer> nums) {
		int best = 0, bestScore = -1;
		for (int n : nums.keySet()) {
			if (nums.get(n) > bestScore) {
				best = n;
				bestScore = nums.get(n);
			}
		}
		return best;
	}

	private static int gcd(Collection<Integer> nums) {
		int res = 0;
		for (int n : nums) {
			res = gcd(res, n);
		}
		return res;
	}

	private static int gcd(int a, int b) {
		while (b != 0) {
			int c = a;
			a = b;
			b = c % b;
		}
		return a;
	}

	public static CommandLine prepareCmd() {
		CommandLine cmdLine = new CommandLine();
		cmdLine.addInt("context", 15, "Length of the context");
		cmdLine.addDouble("threshold", 0.15, "Listen threshold");
		cmdLine.addBoolean("all-periods", false,
				"Show all candidate repeat periods");
		return cmdLine;
	}

	public static void main(String[] args) throws Exception {
		// Get params from users

		CommandLine cmdLine = prepareCmd();

		args = cmdLine.parseLine(args);

		if (args == null || args.length <= 0) {
			System.err
					.println("Usage: java bio.xm.TandemFinder [options] file\n"
							+ cmdLine.usageMessage() + "\n");
			System.exit(1);
		}

		tModel = new ExpertModelTandem(Alphabet.DNA5(), 20,
				cmdLine.getDoubleVal("threshold"));

		// cmdLine.printOptions();

		// SequenceFileReader fReader
		SequenceReader fReader = SequenceReader.getReader(args[0]);

		Sequence s;
		while ((s = fReader.nextSequence(Alphabet.DNA5())) != null) {
			Set<Integer> candidates = classify(s);
			int period = gcd(candidates);
			System.out.print(s.getName() + "\t" + period);
			if (cmdLine.getBooleanVal("all")) {
				// Show all potential candidates
				System.out.print("\t");
				Iterator<Integer> it = candidates.iterator();
				for (int i = 0; i < candidates.size(); ++i) {
					int c = it.next();
					if (i == 0) {
						System.out.print(c);
					} else {
						System.out.print(" " + c);
					}
				}
			}
			System.out.println();
		}
	}
}
