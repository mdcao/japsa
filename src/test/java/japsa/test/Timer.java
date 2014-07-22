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
 * 01/01/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.test;

/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
public class Timer {

	/**
	 * Timer class, largely based on the class of the same name of Lloyd 
	 * Allison (http://allisons.org/ll)
	 * @param args
	 */
	private final long milliSecs0;
	long last;

	public Timer() {
		milliSecs0 = System.currentTimeMillis(); // msecs since Jan. 1 1970
		last = milliSecs0;
		Runtime runtime = Runtime.getRuntime();
		System.out.println("#processors=" + runtime.availableProcessors()
				+ ", maxMemory=" + runtime.maxMemory() / 1000000.0 + " MB");
	}// constructor

	public void mark(String msg) {
		final long now = System.currentTimeMillis();
		System.out.println(msg + ":" + " increment=" + (now - last) / 1000.0
				+ " sec" + " (total=" + (now - milliSecs0) / 1000.0 + ")");
		last = now;
	}// mark()

	
	public void systemStatus(){
		Runtime s_runtime = Runtime.getRuntime();
		long tMem     = s_runtime.totalMemory(),
		      freeMem =s_runtime.freeMemory();	
		long useMem = tMem - freeMem;
		
		System.out.println("Total memory " + tMem + " bytes (" + (tMem >> 20) + "MB)" );
		System.out.println("Free  memory " + freeMem + " bytes (" + (freeMem >> 20) + "MB)" );
		System.out.println("Used  memory " + useMem + " bytes (" + (useMem >> 20) + "MB)" );		
	}
	
	public static void main(String[] args) {
		System.out.println("-- test Timer.java --");
		Timer t = new Timer();
		t.mark("tick");
		final int n = 10000;
		int[][] a = new int[n][n];
		int sum = 0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				sum += a[i][j];
		System.out.println(sum);
		t.mark("tock");
		System.out.println("-- done --");
	}// main()

}
