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
 * 04/01/2012 - Minh Duc Cao:                                         
 *  
 ****************************************************************************/

package japsa.util;

import java.util.Random;




/**
 * Implement a distribution
 * 
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 * 
 */
public class Distribution {
	protected double[] dist;

	public Distribution(int n) {
		dist = new double[n];
	}

	public Distribution(double[] aDist) {
		dist = aDist;
	}

	public Distribution(int[] counts) {
		dist = new double[counts.length];
		for (int i = 0; i < dist.length; i++)
			dist[i] = counts[i];
		normalise();
	}

	public int length() {
		return dist.length;
	}

	public double getMean() {
		double sum = 0;
		for (int i = 0; i < dist.length; i++) {
			sum += dist[i];
		}

		return sum / dist.length;
	}

	public void normalise() {
		double sum = 0;
		for (int i = 0; i < dist.length; i++) {
			sum += dist[i];
		}

		for (int i = 0; i < dist.length; i++) {
			dist[i] /= sum;// let exception happen
		}

	}

	public void setWeight(int idx, double value) {
		dist[idx] = value;
	}

	// Set akll weights to value
	public void setWeights(double value) {
		for (int i = 0; i < dist.length; i++)
			dist[i] = value;
	}

	public void addWeight(int idx, double value) {
		dist[idx] += value;
	}

	public void scale(double factor) {
		for (int i = 0; i < dist.length; i++) {
			dist[i] *= factor;
		}

	}

	public double getWeight(int idx) {
		return dist[idx];
	}

	public String toString() {
		StringBuffer sb = new StringBuffer("(");

		for (int i = 0; i < dist.length; i++) {
			sb.append("," + dist[i]);
		}

		return sb.append(")").toString();
	}

	public byte randomGenerate(Random rd) {
		double x = rd.nextDouble();
		for (byte i = 0; i < dist.length; i++) {
			if (x <= dist[i])
				return i;
			else
				x -= dist[i];
		}
		return 0;
	}

	/**
	 * Precond: aDist is normalise
	 * 
	 * @param aDist
	 * @return
	 */
	public double expectedCode() {
		double msgLen = 0.0;
		for (int i = 0; i < dist.length; i++) {
			msgLen -= JapsaMath.log2(dist[i]) * dist[i];
		}
		return msgLen;
	}

}
