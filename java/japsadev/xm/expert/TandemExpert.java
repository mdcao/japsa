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

package japsadev.xm.expert;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.Sequence;
import japsa.util.JapsaMath;
import japsa.xm.expert.Expert;
import japsadev.xm.ExpertModelTandem;

import java.util.Arrays;

/**
 * Expert predict based on thinking the current base is part of a tandem repeat
 * 
 * @author Minh Duc Cao
 * 
 */
public class TandemExpert extends Expert {
	Sequence seq;
	private int period = 3;
	private int learnWindow = 20;
	private int currentPos = 0;
	private int totalCounts = 4, countRight = 1;

	double totalGain = 2.5, countPredict = 1;
	public double[] tandemCosts;
	public double[] positive;
	int trial = 0;
	int currentRun = 0;
	double currentSum = 0;
	double trialInfo = 0;

	int maxTrial = 3;
	int windowIndex = 0;
	int[] windowCount;

	double sumLastFour = 8;

	public TandemExpert(Sequence seq, int period) {
		this(seq, period, 20);
	}

	public TandemExpert(Sequence seq, int period, int window) {
		super();
		this.seq = seq;
		this.period = period;
		learnWindow = window;

		tandemCosts = new double[seq.length()];
		positive = new double[seq.length()];

		windowCount = new int[learnWindow];
		Arrays.fill(windowCount, -1);

		maxTrial = period * 2 / 3;
		if (maxTrial < 4)
			maxTrial = 4;

		// totalCounts = learnWindow;
		// countRight = totalCounts - 3 * totalCounts / 4;
	}

	@Override
	public double probability(int character) {
		if (currentPos < period)
			return 0.25;

		double prob = (countRight + 0.0) / (totalCounts);
		if (prob <= 0 || prob >= 1) {
			System.err.println("Error " + prob + "  " + totalCounts + "   "
					+ countRight);
			(new Exception()).printStackTrace();
			System.exit(1);
		}

		if (seq.symbolAt(currentPos - period) == character)
			return prob;
		else
			return (1 - prob) / (Expert.alphabet().size() - 1);
	}

	@Override
	public double update(int actual) {
		double prob = probability(actual);

		// update position

		// update counts
		if (currentPos >= period) {
			totalCounts++;
			if (seq.symbolAt(currentPos - period) == actual)
				countRight++;
		}

		if (currentPos >= learnWindow + period) {
			totalCounts--;
			if (seq.symbolAt(currentPos - learnWindow) == seq
					.symbolAt(currentPos - learnWindow - period))
				countRight--;
		}
		currentPos++;

		updateCost(prob);
		return prob;
	}

	public double updatePos(int pos) {
		int actual = seq.symbolAt(pos);
		double prob = probability(actual);

		// if (posSrc == 1000){
		// int x;
		// x = 1;
		// }
		// update counts

		if (currentPos >= period) {
			if (windowCount[windowIndex] == 0) {
				totalCounts--;
			}
			if (windowCount[windowIndex] == 1) {
				totalCounts--;
				countRight--;
			}

			totalCounts++;
			if (seq.symbolAt(currentPos - period) == actual
			// || countRight * 1.6 < totalCounts
			) {
				windowCount[windowIndex] = 1;
				countRight++;
			} else {
				windowCount[windowIndex] = 0;
			}

			windowIndex++;
			if (windowIndex >= windowCount.length)
				windowIndex = 0;
		}

		currentPos++;

		double msgLen = -JapsaMath.log2(prob);
		tandemCosts[pos] = msgLen;

		if (pos >= 4) {
			this.sumLastFour -= tandemCosts[pos - 4];
		} else
			this.sumLastFour -= 2;

		sumLastFour += msgLen;

		if (sumLastFour / 4 < totalGain / countPredict
				- ExpertModelTandem.threshold) {
			trial = 0;
			currentRun++;
			currentSum += msgLen;
			positive[pos] = currentSum / currentRun;
		} else {
			trial++;
			if (trial < maxTrial) {
				currentRun++;
				currentSum += msgLen;
				positive[pos] = currentSum / currentRun;
			} else {
				if (currentRun - trial > ExpertModelTandem.minLen) {
					int start = pos - currentRun - period + 1 - learnWindow;

					if (start < 1)
						start = 1;
					TandemRepeat str = new TandemRepeat(seq.getName(), start,
							pos - maxTrial);
					str.setPeriod(period);
					str.setScore(totalGain / countPredict - currentSum
							/ currentRun);
					str.setUnitNo((double) str.getLength() / period);

					str.setID("M" + str.getStart());

					str.addDesc("@R:" + period);
					// Unit no
					str.addDesc("@N:" + str.getUnitNo());
					// myCost
					str.addDesc("@S:" + str.getScore());

					ExpertModelTandem.anno.add(str);
				}
				positive[pos] = -0.5;
				currentRun = 0;
				currentSum = 0;
			}
		}

		positive[pos] = sumLastFour / 4;

		totalGain += msgLen;
		countPredict++;

		updateCost(prob);
		return prob;
	}

}
