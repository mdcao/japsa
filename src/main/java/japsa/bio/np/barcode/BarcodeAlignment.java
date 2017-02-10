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

/****************************************************************************
 *                           Revision History                                
 * 22 Nov 2016 - Minh Duc Cao: Started                                 
 *  
 ****************************************************************************/
package japsa.bio.np.barcode;


import japsa.seq.Sequence;

/**
 * Implement based on jaligner from Ahmed Moustafa
 * See license below
 */

public final class BarcodeAlignment {

	/**
	 * Traceback direction stop
	 */
	public static final byte STOP = 0;
	/**
	 * Traceback direction left
	 */
	public static final byte LEFT = 1;
	/**
	 * Traceback direction diagonal
	 */
	public static final byte DIAGONAL = 2;
	/**
	 * Traceback direction up
	 */
	public static final byte UP = 3;

	public BarcodeAlignment(Sequence s1, Sequence s2) {
		super();
		this.barcodeSequence = s1;
		this.readSequence = s2;	

		m = s1.length() + 1;
		n = s2.length() + 1;

		//Initilise the arrays
		pointers = new byte[m * n];
		sizesOfVerticalGaps = new short[m * n];
		sizesOfHorizontalGaps = new short[m * n];
	}	
	
	Sequence barcodeSequence;
	Sequence readSequence;
	int m,n;
	byte[] pointers;
	short[] sizesOfVerticalGaps;
	short[] sizesOfHorizontalGaps;
	//BLOSSOM62
	//double [][] scores = 
	//	{{4.0,0.0,0.0,0.0},
	//			{0.0,9.0,-3.0,-1.0}, 
	//			{0.0,-3.0,6.0,-2.0},
	//			{0.0,-1.0,-2.0,5.0}
	//	};
	
	//poreFUME's scores
	double openPenalty = 4.7;
	double extendPenalty = 1.6;

	double [][] scores = {
			{  2.7, -4.5, -4.5, -4.5},
			{ -4.5,  2.7, -4.5, -4.5},
			{ -4.5, -4.5,  2.7, -4.5},
			{ -4.5, -4.5, -4.5,  2.7}			
	};
		
	
	private int cellRow;
	/**
	 * Column of the cell
	 */
	private int cellCol;
	/**
	 * Alignment score at this cell
	 */
	private double cellScore;

	
	public void setBarcodeSequence(Sequence seq){
		barcodeSequence = seq;
	}
	
	public void setReadSequence(Sequence seq){
		readSequence = seq;
	}
	
	

	public double align() {		
		// Initializes the boundaries of the traceback matrix to STOP.
		for (int i = 0, k = 0; i < m; i++, k += n) {
			pointers[k] = STOP;
		}
		for (int j = 1; j < n; j++) {
			pointers[j] = STOP;
		}

		for (int i = 0, k = 0; i < m; i++, k += n) {
			for (int j = 0; j < n; j++) {
				sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
			}
		}
		return construct();
	}

	/**
	 * Constructs directions matrix for the traceback
	 * 
	 * @param barcodeSequence
	 *            sequence #1
	 * @param readSequence
	 *            sequence #2
	 * @param scores
	 *            scoring matrix
	 * @param openPenalty
	 *            open gap penalty
	 * @param extendPenalty
	 *            extend gap penalty
	 * @return The cell where the traceback starts.
	 */
	private double construct() {
		//logger.info("Started...");
		//long start = System.currentTimeMillis();

		double f; // score of alignment x1...xi to y1...yi if xi aligns to yi
		double[] g = new double[n]; // score if xi aligns to a gap after yi
		double h; // score if yi aligns to a gap after xi
		double[] v = new double[n]; // best score of alignment x1...xi to
		// y1...yi
		double vDiagonal;

		g[0] = Float.NEGATIVE_INFINITY;
		h = Float.NEGATIVE_INFINITY;
		v[0] = 0;

		for (int j = 1; j < n; j++) {
			g[j] = Float.NEGATIVE_INFINITY;
			v[j] = 0;
		}

		double similarityScore, g1, g2, h1, h2;

		cellScore = Float.NEGATIVE_INFINITY;
		//Cell cell = new Cell();

		for (int i = 1, k = n; i < m; i++, k += n) {
			h = Float.NEGATIVE_INFINITY;
			vDiagonal = v[0];
			for (int j = 1, l = k + 1; j < n; j++, l++) {
				similarityScore = scores[barcodeSequence.getBase(i-1)][readSequence.getBase(j-1)];

				// Fill the matrices
				f = vDiagonal + similarityScore;

				g1 = g[j] - extendPenalty;
				g2 = v[j] - openPenalty;
				if (g1 > g2) {
					g[j] = g1;
					sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
				} else {
					g[j] = g2;
				}

				h1 = h - extendPenalty;
				h2 = v[j - 1] - openPenalty;
				if (h1 > h2) {
					h = h1;
					sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
				} else {
					h = h2;
				}

				vDiagonal = v[j];
				v[j] = maximum(f, g[j], h, 0);

				// Determine the traceback direction
				if (v[j] == 0) {
					pointers[l] = STOP;
				} else if (v[j] == f) {
					pointers[l] = DIAGONAL;
				} else if (v[j] == g[j]) {
					pointers[l] = UP;
				} else {
					pointers[l] = LEFT;
				}

				// Set the traceback start at the current cell i, j and score
				if (v[j] > cellScore) {
					cellRow  = i;
					cellCol = j;
					cellScore = v[j];
					//cell.set(i, j, v[j]);
				}
			}
		}		
		return cellScore;
	}


	/**
	 * Returns the maximum of 4 float numbers.
	 * 
	 * @param a
	 *            float #1
	 * @param b
	 *            float #2
	 * @param c
	 *            float #3
	 * @param d
	 *            float #4
	 * @return The maximum of a, b, c and d.
	 */
	private static double maximum(double a, double b, double c, double d) {
		if (a > b) {
			if (a > c) {
				return a > d ? a : d;
			} else {
				return c > d ? c : d;
			}
		} else if (b > c) {
			return b > d ? b : d;
		} else {
			return c > d ? c : d;
		}
	}

}

/**
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
