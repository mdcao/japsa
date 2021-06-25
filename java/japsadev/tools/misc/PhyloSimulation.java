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

package japsadev.tools.misc;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.util.Distribution;
import japsa.util.JapsaMath;

import java.util.Random;




//import Jama.Matrix;


/**
 * @author Minh Duc Cao www.caominhduc.org
 * Simulation used in XMas
 */
public class PhyloSimulation {
	public static int LENGTH = 10000;// Length = 1000

	static double[] fre = { 0.1, 0.4, 0.4, 0.1 };

	// static double [] fre = {0.25,0.25,0.25,0.25};
	static double alpha = 0.001;

	static double ga = fre[0], gc = fre[1], gg = fre[2], gt = fre[3];
	static double beta = alpha / 2.0;

	public static double[][] unit = { { 1.0, 0.0, 0.0, 0.0 },
			{ 0.0, 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0, 0.0 },
			{ 0.0, 0.0, 0.0, 1.0 } };

	public static double[][] jukes_cantor = {
			{ 1 - 3 * alpha, alpha, alpha, alpha },
			{ alpha, 1 - 3 * alpha, alpha, alpha },
			{ alpha, alpha, 1 - 3 * alpha, alpha },
			{ alpha, alpha, alpha, 1 - 3 * alpha } };

	public static double[][] kimura = {
			{ 1 - alpha - 2 * beta, beta, alpha, beta },
			{ beta, 1 - alpha - 2 * beta, beta, alpha },
			{ alpha, beta, 1 - alpha - 2 * beta, beta },
			{ beta, alpha, beta, 1 - alpha - 2 * beta } };

	public static double[][] equal_input = {
			{ 1 - alpha * (1 - ga), alpha * gc, alpha * gg, alpha * gt },
			{ alpha * ga, 1 - alpha * (1 - gc), alpha * gg, alpha * gt },
			{ alpha * ga, alpha * gc, 1 - alpha * (1 - gg), alpha * gt },
			{ alpha * ga, alpha * gc, alpha * gg, 1 - alpha * (1 - gt) } };

	public static double[][] hky = {
			{ 1 - beta * gc - alpha * gg - beta * gt, beta * gc, alpha * gg,
					beta * gt },
			{ beta * ga, 1 - beta * ga - beta * gg - alpha * gt, beta * gg,
					alpha * gt },
			{ alpha * ga, beta * gc, 1 - alpha * ga - beta * gc - beta * gt,
					beta * gt },
			{ beta * ga, alpha * gc, beta * gg,
					1 - beta * ga - alpha * gc - beta * gg } };

	static double pa = 0.001, pb = 0.003, pc = 0.002, pd = 0.0015, pe = 0.0025,
			pf = 0.0035;
	public static double[][] rev = {
			{ 1 - (pa * gc + pb * gg + pc * gt), pa * gc, pb * gg, pc * gt },
			{ pa * ga, 1 - (pa * ga + pd * gg + pe * gt), pd * gg, pe * gt },
			{ pb * ga, pd * gc, 1 - (pb * ga + pd * gc + pf * gt), pf * gt },
			{ pc * ga, pe * gc, pf * gg, 1 - (pc * ga + pe * gc + pf * gg) } };

	/*********************************************************************************
	 * public static double [][] rateMtx = {{0.9997,0.0001,0.0001,0.0001},
	 * {0.0001,0.9997,0.0001,0.0001}, {0.0001,0.0001,0.9997,0.0001},
	 * {0.0001,0.0001,0.0001,0.9997} };
	 * /****************************************
	 * ***************************************** public static double [][]
	 * rateMtx = {{0.999,0.0002,0.0005,0.0003}, {0.0003,0.999,0.0002,0.0005},
	 * {0.0005,0.0003,0.999,0.0002}, {0.0002,0.0005,0.0003,0.999} };
	 * /***********
	 * **********************************************************************
	 * 
	 * /
	 *********************************************************************************/
	public static double[][] rateMtx = 
		  { { 0.999, 0.0002, 0.0005, 0.0003 },
			{ 0.0005, 0.9985, 0.0004, 0.0006 },
			{ 0.0004, 0.0007, 0.9983, 0.0006 },
			{ 0.0003, 0.0004, 0.0003, 0.999 } };

	/*********************************************************************************
	 * public static double [][] rateMtx = {{0.99,0.003,0.005,0.002},
	 * {0.003,0.99,0.002,0.005}, {0.005,0.002,0.99,0.003},
	 * {0.002,0.006,0.002,0.99} }; /
	 *********************************************************************************/

	/**
	 * matrix a (mxm) times b (mxm) -- a slow implemetation
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static double[][] times(double[][] a, double[][] b) {
		int m = a.length;
		double[][] c = new double[m][m];

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				c[i][j] = 0;
				for (int k = 0; k < m; k++)
					c[i][j] += a[i][k] * b[k][j];
			}
		}
		return c;
	}

	public static double[] times(double[] a, double[][] b) {
		int m = a.length;
		double[] c = new double[m];

		for (int i = 0; i < 1; i++) {
			for (int j = 0; j < m; j++) {
				c[j] = 0;
				for (int k = 0; k < m; k++)
					c[j] += a[k] * b[k][j];
			}
		}
		return c;
	}

	public static double[][] copy(double[][] a) {
		int m = a.length;
		double[][] c = new double[m][m];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				c[i][j] = a[i][j];
			}
		}
		return c;
	}

	public static double entropy(double[] s) {
		double ent = 0;
		for (int i = 0; i < s.length; i++)
			ent -= s[i] * JapsaMath.log2(s[i]);

		return ent;
	}

	public static void generate(int y, int z) throws Exception {

		System.out.printf("%4d %4d %f\n", y, z, y * 1.0 / z);
		int length = 10000;

		double[] aGeneDist = { .2, .3, .3, .2 };
		japsa.util.Distribution dist = new japsa.util.Distribution(aGeneDist);

		byte[] seqX = new byte[length], seqZ = new byte[length], seqY = new byte[length];

		for (int i = 0; i < seqX.length; i++) {
			seqX[i] = dist.randomGenerate(rnd);
		}

		japsa.util.Distribution[] mtY = new japsa.util.Distribution[4], mtZ = new japsa.util.Distribution[4];

		double[][] c = copy(rateMtx);
		for (int i = 0; i < y; i++)
			// y times
			c = times(c, rateMtx);

		for (int i = 0; i < c.length; i++) {
			double s = 0;
			for (int j = 0; j < c[i].length; j++) {
				s += c[i][j];
				System.out.printf("%f   ", c[i][j]);
			}
			System.out.println("   " + s);
		}

		for (int i = 0; i < 4; i++) {
			mtY[i] = new Distribution(c[i]);
		}

		c = copy(rateMtx);
		for (int i = 0; i < z; i++)
			// y times
			c = times(c, rateMtx);

		for (int i = 0; i < c.length; i++) {
			double s = 0;
			for (int j = 0; j < c[i].length; j++) {
				s += c[i][j];
				System.out.printf("%f   ", c[i][j]);
			}
			System.out.println("   " + s);
		}

		for (int i = 0; i < 4; i++) {
			mtZ[i] = new Distribution(c[i]);
		}

		for (int i = 0; i < seqX.length; i++) {
			seqY[i] = mtY[seqX[i]].randomGenerate(rnd);
			seqZ[i] = mtZ[seqX[i]].randomGenerate(rnd);
		}

		Sequence Y = new Sequence(Alphabet.DNA4(), seqY), X = new Sequence(
				Alphabet.DNA4(), seqX), Z = new Sequence(Alphabet.DNA4(), seqZ);

		X.writeFasta("seqX");
		Y.writeFasta("seqY");
		Z.writeFasta("seqZ");

		System.out
				.println("=======================================================");
	}
	

	static Random rnd = new Random(13);

	
	
	public static byte[] genDis(byte[] src, int dis) {
		byte[] target = new byte[src.length];
		// Random rnd = new Random();

		double[][] ct = copy(rateMtx);

		for (int i = 1; i < dis; i++)
			// y times
			ct = times(ct, rateMtx);

		japsa.util.Distribution[] pvGene = new japsa.util.Distribution[4];
		for (int i = 0; i < pvGene.length; i++) {
			pvGene[i] = new japsa.util.Distribution(ct[i]);
		}

		for (int i = 0; i < src.length; i++) {
			target[i] = pvGene[src[i]].randomGenerate(rnd);
		}
		return target;
	}

	public static double disVector2Matrix(double[] v, double[][] mt) {
		double res = 0.0;
		for (int i = 0; i < mt.length; i++) {
			for (int j = 0; j < mt[i].length; j++) {
				res -= v[i] * mt[i][j] * JapsaMath.log2(mt[i][j]);
			}
		}
		return res;
	}

	public static void main4() throws Exception {
		double scale = 1;
		double[] s = { 0.1, 0.2, 0.3, 0.4 };

		// the subs matrix
		double[][] ct = copy(rateMtx);
		for (int i = 0; i < 1800; i++) {// y times
			ct = times(ct, rateMtx);

			// target
			double[] t = times(s, ct);

			// if (i % 100 ==0)
			double I_t = entropy(t);
			double I_ts = disVector2Matrix(s, ct);

			double I_s = entropy(s);
			double I_st = disVector2Matrix(t, ct);

			double dis1 = -JapsaMath.log2((I_t - I_ts) / (I_t + I_ts));

			double dis2 = -JapsaMath.log2((2 - I_st) / (2 + I_st));

			double dis = -JapsaMath.log2((I_s - I_st + I_t - I_ts)
					/ (I_s + I_st + I_t + I_ts));

			System.out.println(i / scale + " " + I_st + " " + I_s + " "
					+ (I_s - I_st) + " " + dis1 + " " + dis2 + " " + dis);

		}
	}

	

	public static double nlogn(double x) {
		return 6 * x * (x - 1) * JapsaMath.log2e / (x + 4 * Math.sqrt(x) + 1);
	}

	public static double calx(double p) {
		return JapsaMath.log2e
				* 2
				* (1 - p)
				* 9
				/ 4
				* ((1 + 3 * p) / (5 + 3 * p + 8 * Math.sqrt(1 + 3 * p)) + (p + 3)
						/ (5 - p + 8 * Math.sqrt(1 - p)));
	}

	public static void test2() {
		for (int t = 1; t < 100; t++) {

			double p = Math.exp(-t * 0.04);
			double d = (1.0 + 3.0 * p) / 4.0;

			double od = (1.0 - p) / 4.0;

			double y = -(d * JapsaMath.log2(d) + 3 * od * JapsaMath.log2(od));

			// double x = MyMath.log2((2 - y)/(2 + y));

			double z = 2 - 2 * Math.pow(p, 2);
			// System.out.println(t + " " + y + " " + z);

			double app = -nlogn(d) - 3 * nlogn(od);
			// (1 + 3 * p)/(5 + 3 * p + 8 * Math.sqrt(1 + 3*p))
			// + (3+p)/(5 - p + 8 * Math.sqrt(1 - p));

			// app = app * 9 * (1 - p)/MyMath.loge2 / 2;

			// System.out.println(y/(1 - p) + "  " + app/(1 - p) + " " + z/(1 -
			// p) +" " + calx(p)/(1 - p));
			System.out.println(y + "  " + app + " " + z + " " + calx(p) + " "
					+ (1 - p));
		}
	}

	public static void testDistance() {

		// rateMtx = jukes_cantor;
		rateMtx = kimura;
		// rateMtx = hky;
		// rateMtx = rev;
		// rateMtx = equal_input;

		double[][] mx = copy(unit);

		double[][] my = copy(unit);
		// double I_s = entropy(fre);

		for (int t = 1; t < 2000; t++) {
			mx = times(mx, rateMtx);

			// time my twice
			my = times(my, rateMtx);
			// my = times(my, rateMtx);

			double[] fx = times(fre, mx);
			double[] fy = times(fre, my);

			double I_x = entropy(fx);
			double I_y = entropy(fy);

			// double I_xs = disVector2Matrix(fre, mx);
			// double I_sx = disVector2Matrix(fx, mx);

			// System.out.println(-MyMath.log2( (I_x + I_s - I_sx - I_xs)/(I_s +
			// I_x) ));

			double I_xy = disVector2Matrix(fy, times(mx, my));

			double I_yx = disVector2Matrix(fx, times(my, mx));

			System.out.println(t +" " + (-JapsaMath.log2((I_x + I_y - I_xy - I_yx)
					/ (I_x + I_y))));

		}
	}
	
	/**
	 * Compare 
	 */
	public static void testFunction() {
		double x = 0.9999;
		double p = 1;
		for (int t = 1; t < 50000; t++) {
			p *= x;
			System.out.println(2 - (-(1 + 3 * p) / 4
					* JapsaMath.log2((1 + 3 * p) / 4) - 3 * (1 - p) / 4
					* JapsaMath.log2((1 - p) / 4)));
		}
	}

	public static void functionFitting() {
		double p;
		//alpha = 0.0001;
		for (int t = 0; t < 5000; t++) {
			double a = -8*t*alpha;
			p = Math.exp(a);

			double y = (1.0 - p) / 4.0;
			double x = (1.0 + 3 * p) / 4.0;
			double r = -3 * y * JapsaMath.log2(y) - x * JapsaMath.log2(x);

			// System.out.println(r);
			System.out.println(t + " " + r / 2 + " " + (1 - p * p)+" " + (1 - p*p - r/2));
		}
	}

	public static void main(String[] args) throws Exception {
		// main2Seq();
		// testPriors(args);

		// threeSome(args);
		// testFunction();
		// main4();
		// test2();
		// mainTemp2();

		 testDistance();

		 //functionFitting();
		 //testFunction();
	}
}
