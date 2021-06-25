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
 * 28/02/2016 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsadev.tools.work;


import java.io.BufferedReader;
import java.io.IOException;

import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.DoubleArray;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Compute cross correllation between two datasets in bedgraph format
 * 
 * 
 * @author minhduc
 *
 */
@Deployable(
		scriptName = "jsa.dev.plasmaCrossCor", 
		scriptDesc = "Analysis of plasma sequencing using Cross correlation"
		)
public class PlasmaAnalysisCrossCorrelationCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(PlasmaAnalysisCrossCorrelationCmd.class);

	//CommandLine cmdLine;
	public PlasmaAnalysisCrossCorrelationCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("xFile", null, "File 1",true);		
		addString("yFile", null, "File 2",true);		
		addInt("window", 500, "Window sise");
		addInt("lag", 10, "lag");
		addString("output", null, "Name of the output file",true);

		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		PlasmaAnalysisCrossCorrelationCmd cmdLine = new PlasmaAnalysisCrossCorrelationCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/
		String xFile  = cmdLine.getStringVal("xFile");
		//String controlBam = cmdLine.getStringVal("controlBam");
		String yFile = cmdLine.getStringVal("yFile");
		String output = cmdLine.getStringVal("output");		
		int window = cmdLine.getIntVal("window");
		int lag = cmdLine.getIntVal("lag");

		LOG.info("read file 1");
		BufferedReader bf = SequenceReader.openFile(xFile);
		DoubleArray array = new DoubleArray();

		int start = -1;
		String chr = "";
		String line = "";
		while ( (line = bf.readLine())!= null){
			// skip the bedGraph file header
			if (line.startsWith("track type=")){
				continue;
			}
			
			String [] toks  = line.split("\t");
			array.add(Double.parseDouble(toks[3]));

			if (start < 0){
				start = Integer.parseInt(toks[1]);
				chr = toks[0];
			}
		}
		bf.close();

		double [] x = array.toArray();
		array.clear();

		LOG.info("read file 2");
		bf = SequenceReader.openFile(yFile);		
		while ( (line = bf.readLine())!= null){
			// skip the header
			if (line.startsWith("track type=")){
				continue;
			}
			
			String [] toks  = line.split("\t");
			array.add(Double.parseDouble(toks[3]));			
		}
		bf.close();

		double [] y = array.toArray();

		double [] crr = new double[x.length];
		LOG.info("Run 0");
		cross_correlation(x,y,window,0,crr);		
		for (int i=1; i < lag; i++){
			LOG.info("Run " + i);
			cross_correlation(x,y,window,i,crr);
			LOG.info("Run -" + i);
			cross_correlation(y,x,window,i,crr);		
		}

		LOG.info("Write");
		SequenceOutputStream fCount = SequenceOutputStream.makeOutputStream(output);
		fCount.print("track type=bedGraph\n");		
		char sep = '\t';		
		for (int i = 0; i < crr.length;i++){
			fCount.print(chr);
			fCount.print(sep);
			fCount.print(i + start);
			fCount.print(sep);						
			fCount.print(i+1 + start);
			fCount.print(sep);			
			fCount.print(crr[i]);
			fCount.print('\n');			
		}
		fCount.close();
	}
	
	/**
	 * Compute cross correlation between two arrays of double 
	 * @param x
	 * @param y
	 * @param windows
	 * @param lag
	 * @param results
	 * @return
	 */

	public static double[] cross_correlation(double [] x, double [] y, int windows, int lag, double [] results){
		//maxdelay=20				
		int length=x.length;
		//double [] xcorr = new double[length];
		double xSum = 0, ySum = 0, xSq = 0, ySq = 0;

		//the first windows		
		for (int i = 0; i < windows - 1; i++){
			xSum += x[i];
			xSq  += x[i] * x[i];

			ySum += y[i + lag];
			ySq  += y[i + lag] * y[i + lag];			
		}		


		for (int start = 0; start < length - windows - lag; start++){
			int end = start + windows - 1;
			xSum += x[end];
			xSq  += x[end] * x[end];

			ySum += y[end + lag];
			ySq  += y[end + lag] * y[end + lag];

			//mean
			double mx = xSum/windows;
			double my = ySum/windows;			
			double denom = Math.sqrt((xSq - windows * mx * mx) * (ySq - windows * my * my));

			double sum = 0;
			for (int i = 0; i < windows;i++){
				sum += (x[start + i] - mx) * (y [start + i + lag] - my);
			}

			double xcorr  = sum / denom;//Math.abs(sum / denom);
			if (xcorr > results[start])
				results[start] = xcorr;

			xSum -= x[start];
			xSq  -= x[start] * x[start];

			ySum -= y[start + lag];
			ySq  -= y[start + lag] * y[start + lag];
		}//fo
		return results;
	}
}
