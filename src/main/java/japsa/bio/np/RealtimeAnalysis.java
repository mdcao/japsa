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

/*****************************************************************************
 *                           Revision History                                
 * 5 Sep 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsa.bio.np;

import java.util.Date;

import japsa.util.Logging;

/**
 * Real time analysis in a thread that runs in parallel with the thread conecting
 * data.
 *  
 * @author minhduc
 *
 */
public abstract class RealtimeAnalysis implements Runnable {

	private int readPeriod = 0;//Min number of reads before a new analysis
	private int timePeriod = 0;//Min number of mini-seconds before a new analysis	
	private int powerNap = 1000;//sleep time in miniseconds (1 second by default)

	RealtimeAnalysis(){
	}

	private boolean waiting = true;

	long lastTime;//The last time an analysis is done
	long startTime;
	int lastReadNumber;
	String timeNow;

	public void stopWaiting(){	
		Logging.info("All reads received at " + new Date());
		waiting = false;
		//TODO: implement notifies if I am sleeping
	}

	/* (non-Javadoc)
	 * @see java.lang.Runnable#run()
	 */
	@Override
	public void run() {
		startTime = System.currentTimeMillis();
		Logging.info("Start analysing data at " + new Date(startTime));

		try {
			Thread.sleep(timePeriod);
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}

		while (waiting){
			//Need to wait if timing is not right
			long timeSleep = timePeriod - (System.currentTimeMillis() - lastTime);
			if (timeSleep > 0){
				try {
					Logging.info("Not due time, sleep for " + timeSleep/1000.0 + " seconds");
					Thread.sleep(timeSleep);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
				continue;//the next round I wont be here again
			}
			//assert: time satisfied
			int currentRead = getCurrentRead();
			if (currentRead - lastReadNumber < readPeriod){
				try {
					//Logging.info("Not due read (" + currentRead + "), sleep for " + powerNap/1000.0 + " minutes");
					Thread.sleep(powerNap);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
				continue;//the next round I wont be here again
			}
			lastReadNumber = currentRead;
			lastTime = System.currentTimeMillis();
			//assert: read number satisfied
			timeNow = new Date(lastTime).toString();
			analysis();
			Logging.info("RUNTIME\t" + timeNow  + "\t" + (this.lastTime - this.startTime)/1000.0 + "\t" + this.lastReadNumber + "\t" + (System.currentTimeMillis() - lastTime)/1000.0);
		}//while

		//perform the final analysis
		lastTime = System.currentTimeMillis();
		lastReadNumber =  getCurrentRead();
		timeNow = new Date(lastTime).toString();
		analysis();
		Logging.info("RUNTIME\t" + timeNow  + "\t" + (this.lastTime - this.startTime)/1000.0 + "\t" + this.lastReadNumber + "\t" + (System.currentTimeMillis() - lastTime)/1000.0);
		//.. and close it
		close();
		Logging.info("Real time analysis done");
	}

	abstract protected void close();
	abstract protected void analysis();
	abstract protected int getCurrentRead();	

	/**
	 * Set the minimum number of reads before a new analysis 
	 * @param readNumber
	 */
	public void setReadPeriod(int readNumber){
		readPeriod = readNumber;
	}
	/**
	 * Set the minimum time before a new analysis
	 * @param timePeriod: number of miniseconds
	 */
	public void setTimePeriod(int timePeriod){
		this.timePeriod = timePeriod;
	}

	void setPowerNap(int powerNap){
		this.powerNap = powerNap;
	}	
}
