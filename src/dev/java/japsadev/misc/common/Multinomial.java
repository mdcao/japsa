/*  
 *  Copyright (c) David Powell <david@drp.id.au>
 *
 * 
 * This file is used by both FuzzyLZ and AlignCompress

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


 */

package japsadev.misc.common;

public class Multinomial {
	/**
	 * Compute the MML estimate for encoding the parameters to the multinomial.
	 * p[] are the probabilites (should sum to 1) N is the number of data items
	 * that will be encoded.
	 * 
	 * Encoding length of the Multinomial parameters are returned in bits.
	 */
	static public double MMLparameter_cost(double[] p, double N) {

		double h = MyMath.factorial(p.length - 1); // h(theta) - prior
													// probabilty density =
													// (K-1)!
		double F = 1.0 / p[0]; // F will be the Fischer = N^(K-1)/(p1*p2*...*pk)
		for (int i = 1; i < p.length - 1; i++) {
			F *= N / p[i];
		}

		double cost = 0.5 * MyMath.log2(1 + F
				/ (h * h * Math.pow(12, p.length - 1)));

		cost += 0.5 * (p.length - 1) * MyMath.log2e;

		return cost;
	}
}
