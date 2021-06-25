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
 * 08/04/2012 - Minh Duc Cao: Revised                                        
 * 16/11/2013 - Minh Duc Cao 
 ****************************************************************************/

package japsadev.tools.misc;

import japsa.bio.tr.TandemRepeatVariant;
import japsa.bio.tr.TandemRepeat;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

import java.io.BufferedReader;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * @author minhduc
 * C
 */
public class CollectSTRV {

	/**
	 * Collect variation from many individuals and compute the variability
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		
		System.err.println("Process " + (args.length) + " files");
		BufferedReader [] bin = new BufferedReader[args.length];
		double qual = 0.5;
		
		
		SequenceOutputStream outS = SequenceOutputStream.makeOutputStream("output.strvar");
		SequenceOutputStream outbed = SequenceOutputStream.makeOutputStream("output.bed");
		outS.print("#H:chr\tstart\tend\tperiod\tvar\tconfidence\n");
		
		DescriptiveStatistics stats = new DescriptiveStatistics();
		
		for (int i = 0; i < bin.length;i++){
			System.err.println(args[i]);
			bin[i] = SequenceReader.openFile(args[i]);
		}
		
		String line = bin[0].readLine();
		String [] headers = line.trim().substring(3). split("\t");		
		
		for (int i = 1; i < bin.length;i++){
			line = bin[i].readLine();
		}
		
		
		while (true){
			int ssss = -1;
			stats.clear();
			TandemRepeatVariant aSTRV = null;
			for (int i = 0; i < bin.length;i++){
				line = bin[i].readLine();
				if (line == null)
					break;
				
				aSTRV = TandemRepeatVariant.read(line, headers);
				
				if (ssss < 0)
					ssss = aSTRV.getStart();
				
				if (ssss != aSTRV.getStart()){
					System.err.println("ERROR " + aSTRV);
					System.exit(1);
				}
				
				if (aSTRV.getConfidence() > qual){
					stats.addValue(aSTRV.getVar() * aSTRV.getPeriod());
					stats.addValue(aSTRV.getVar2() * aSTRV.getPeriod());
				}	
			}//for
			if (aSTRV == null)
				break;//while
			double var = 0, conf = 0;
			if (stats.getN() > 0){
				var =  stats.getStandardDeviation() ;
				conf = stats.getN();
			}
			outS.print(aSTRV.getChr()+"\t"+aSTRV.getStart()+"\t"+aSTRV.getEnd()+"\t" +aSTRV.getPeriod() + 
					"\t" + var + "\t"  + conf + "\n");
			
			TandemRepeat str = aSTRV.getTandemRepeat();
			str.setScore(var);
			str.writeBED(outbed);			
			
		}//while
		outS.close();
		outbed.close();
	}
}
