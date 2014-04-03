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
 * 11/01/2012 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsa.util;


/**
 * Some mathematical operations.
 * 
 * @author Minh Duc Cao
 * 
 */
public class JapsaMath {
	
	//Avoid creating object
	private JapsaMath(){
		
	}

	public static final double loge2 = Math.log(2);	
	public static final double log2e = 1.0 / loge2;	

	public static double log2(double a) {
		return Math.log(a) * log2e;
	}

	public static double exp2(double a) {
		return Math.exp(a * loge2);
	}
	
	public static double sqr(double x){
		return x * x;
	}

	/*************************************************************************
	 * Operations on array
	 ************************************************************************/	
	/**
	 * Compute the mean of an array of double values
	 * @param values
	 * @return The mean of the array
	 */
	public static double mean(double[] values) {
		double sum = 0;
		for (int i = 0; i < values.length; i++)
			sum += values[i];
		if (values.length > 0)
			sum /= values.length;

		return sum;
	}
	
	/**
	 * Scale an array of double so that the lowest value becomes 0 and the 
	 * highest becomes 1.
	 * @param values
	 */
	public static void scale(double[] values) {
		double min = values[0], max = values[0];
		for (int i = 1; i < values.length; i++) {
			if (values[i] < min)
				min = values[i];
			if (values[i] > max)
				max = values[i];
		}
		scale(values, min, max);
	}
	
	/**
	 * Scale an array of numbers so that the value low becomes 0 and high becomes 1. 
	 * @param values
	 * @param low
	 * @param high
	 */
	public static void scale(double[] values, double low, double high) {
		double diff = high - low;
		if (diff == 0)
			return;
		for (int i = 0; i < values.length; i++) {
			values[i] = (values[i] - low) / diff;
		}
	}
	
	/**
	 * Compute the statistics (mean and sd) of a sample
	 * @return
	 */
	public static double [] getStat(double[] inp){
		return getStat(inp, inp.length);	}
	/**
	 * Get mean and standard deviations of the first n elements of array
	 * @param inp
	 * @param len
	 * @return
	 */
	public static double [] getStat(double[] inp, int len){
		double [] stat = new double[2];
		
		if (len > inp.length)
			len = inp.length;
		
		if (len == 0) return stat;
		double sum = 0;
		
		for (int i = 0; i < len; i++)
			sum += inp[i];
		
		stat[0] = sum / len;
		
		sum = 0;
		for (int i = 0; i < len; i++)
			sum += (inp[i] - stat[0]) * (inp[i] - stat[0]);
		
		stat[1] = Math.sqrt(sum / len);
				
		return stat;		
	}

	/*************************************************************************
	 * Math operations for matrices
	 * 
	 *************************************************************************/
	
	
	/**
	 * Multiply two matrices. This may not be the best implementation of matrix
	 * multiplication.
	 * 
	 * @param a
	 * @param b
	 * @return the product of the two matrix
	 */

	public static double[][] timeMatrix(double[][] a, double[][] b) {
		int m = a.length, n = b[0].length;
		double[][] c = new double[m][n];
		for (int x = 0; x < m; x++) {
			for (int y = 0; y < n; y++) {
				c[x][y] = 0;
				for (int t = 0; t < a[0].length; t++) {
					c[x][y] += a[x][t] * b[t][y];
				}
			}
		}
		return c;
	}
	
	
	final static byte [] DigitTens = {
		'0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
		'1', '1', '1', '1', '1', '1', '1', '1', '1', '1',
		'2', '2', '2', '2', '2', '2', '2', '2', '2', '2',
		'3', '3', '3', '3', '3', '3', '3', '3', '3', '3',
		'4', '4', '4', '4', '4', '4', '4', '4', '4', '4',
		'5', '5', '5', '5', '5', '5', '5', '5', '5', '5',
		'6', '6', '6', '6', '6', '6', '6', '6', '6', '6',
		'7', '7', '7', '7', '7', '7', '7', '7', '7', '7',
		'8', '8', '8', '8', '8', '8', '8', '8', '8', '8',
		'9', '9', '9', '9', '9', '9', '9', '9', '9', '9',
		} ; 

	    final static byte [] DigitOnes = { 
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		} ;
	    
	final static byte[] digits = {
		'0' , '1' , '2' , '3' , '4' , '5' ,
		'6' , '7' , '8' , '9' , 'a' , 'b' ,
		'c' , 'd' , 'e' , 'f' , 'g' , 'h' ,
		'i' , 'j' , 'k' , 'l' , 'm' , 'n' ,
		'o' , 'p' , 'q' , 'r' , 's' , 't' ,
		'u' , 'v' , 'w' , 'x' , 'y' , 'z'
	    };
	
	/**
	 * Count the number of bytes require to represent number num in base 10
	 * @param num
	 * @return
	 */
	public static int radixSize(int num){
		int count = 1;
		
		if (num == Integer.MIN_VALUE){
			return 10;//FIXME: need to check this
        }
		
		if (num < 0){
			count = 2;
			num = -num;
		}
		
		while (num >= 10){
			count ++;
			num /= 10;
		}
		
		return count;
	}
	
	/**
	 * Print the number num into space of fm char, right justified. This method
	 * is largely based on the format integer method from the Java Standard 
	 * implementation.
	 *
	 *  @param num The number to print
	 * @param fm The number of spaces to represent the number
	 * @return
	 */
	public static byte[] formatLong(long num, int fm){		
		byte [] buf = new byte[fm];
		
		long q, r;
        int charPos = fm;//the last pos
        byte sign = 0;
        
        if (num == Integer.MIN_VALUE){
        	///
        }
        
        if (num < 0) { 
            sign = '-';
            num = -num;
        }
        
        // Generate two digits per iteration
        while (num >= 65536 && charPos > 0) {
            q = num / 100;
        // really: r = i - (q * 100);
            r = num - ((q << 6) + (q << 5) + (q << 2));
            num = q;
            buf [--charPos] = DigitOnes[(int) r];
            buf [--charPos] = DigitTens[(int) r];
        }
        
        while (num > 0 && charPos > 0) {
        	q = (num * 52429) >>> (16+3);
            r = num - ((q << 3) + (q << 1));  // r = i-(q*10) ...
            buf [--charPos] = digits [(int) r];
            num = q; 
        }
        
        if (sign != 0 && charPos > 0) {
            buf [--charPos] = sign;
        }        
        while (charPos > 0)
        	buf [--charPos] = ' ';
		
		
        return buf;
	}	
	/**
	 * Print the number num into space of fm char, right justified. This method
	 * is largely based on the format integer method from the Java Standard 
	 * implementation.
	 *
	 *  @param num The number to print
	 * @param fm The number of spaces to represent the number
	 * @return
	 */
	public static byte[] formatInt(int num, int fm){		
		byte [] buf = new byte[fm];
		
		int q, r;
        int charPos = fm;//the last pos
        byte sign = 0;
        
        if (num == Integer.MIN_VALUE){
        	///
        }
        
        if (num < 0) { 
            sign = '-';
            num = -num;
        }
        
        // Generate two digits per iteration
        while (num >= 65536 && charPos > 0) {
            q = num / 100;
        // really: r = i - (q * 100);
            r = num - ((q << 6) + (q << 5) + (q << 2));
            num = q;
            buf [--charPos] = DigitOnes[r];
            buf [--charPos] = DigitTens[r];
        }
        
        while (num > 0 && charPos > 0) {
        	q = (num * 52429) >>> (16+3);
            r = num - ((q << 3) + (q << 1));  // r = i-(q*10) ...
            buf [--charPos] = digits [r];
            num = q; 
        }
        
        if (sign != 0 && charPos > 0) {
            buf [--charPos] = sign;
        }        
        while (charPos > 0)
        	buf [--charPos] = ' ';
		
		
        return buf;
	}	
	
	/**
	 * Compute the phred score of a likelihood
	 * Q = -10log_{10} P
	 * Note: This should take the error probability
	 * @param prob
	 * @return
	 */
	public static double prob2phred(double prob){
		return -10.0 * Math.log10(prob);
	}
	
	
	/**
	 *  P = 10^(-Q/10)
	 * @param phred
	 * @return
	 */	
	public static double phred2prob(double phred){
		return Math.pow(10, -phred/10);
	}
}
