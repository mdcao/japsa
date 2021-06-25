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
 *  
 ****************************************************************************/

package japsadev.tools.misc;

import japsa.bio.tr.TandemRepeatVariant;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

import java.io.BufferedReader;


/**
 * Combine strv from various tools for an individual
 * 
 * @author minhduc
 *  
 */

public class CombineSTRV {

	public static void main(String[] args) throws Exception {		
		System.err.println("Process " + (args.length - 1) + " files");
		BufferedReader [] bin = new BufferedReader[args.length -1 ];
		int lineNo = 1;
		double qual = 0.5;
		
		SequenceOutputStream outS =SequenceOutputStream.makeOutputStream(args[0] + ".strvc");
		outS.print("#H:chr\tstart\tend\tperiod\tvar\tvar2\tconfidence\n");
		
		for (int i = 0; i < bin.length;i++)
			bin[i] = SequenceReader.openFile(args[i+1]);
		
		String line = bin[0].readLine();
		String [] headers = line.trim().substring(3). split("\t");		
		
		for (int i = 1; i < bin.length;i++){
			line = bin[i].readLine();
		}
		
		while (true){
			line = bin[0].readLine();
			if (line == null)		
				break;
			
			lineNo ++;			
						
			TandemRepeatVariant strv = TandemRepeatVariant.read(line, headers);
			if (strv.getVar() > strv.getVar2()){
				strv.swapVar(); 
			}
			if (strv.getConfidence() > 1)
				strv.setConfidence(1);
			
			int count = 0;
			
			if (strv.getConfidence() > qual)
				count = 1;
			
			for (int i = 1; i < bin.length;i++){
				line = bin[i].readLine();
				TandemRepeatVariant aSTRV = TandemRepeatVariant.read(line, headers);
				if (strv.getStart() != aSTRV.getStart() || strv.getEnd() != aSTRV.getEnd())
					throw new RuntimeException("Error at line " + lineNo);
				
				if (aSTRV.getConfidence() > 0.5){
					if (aSTRV.getVar() > aSTRV.getVar2()){
						aSTRV.swapVar(); 
					}
					if (aSTRV.getConfidence() > 1)
						aSTRV.setConfidence(1);
					
					strv.setVar(strv.getVar() + aSTRV.getVar());
					strv.setVar2(strv.getVar2() + aSTRV.getVar2());
					strv.setConfidence(strv.getConfidence() + aSTRV.getConfidence());
					
					count ++;
				}	
			}//for
			if (count > 0){
				strv.setVar(strv.getVar() / count);
				strv.setVar2(strv.getVar2() / count);
				strv.setConfidence(strv.getConfidence() / count);
			}
			outS.print(strv.getChr()+"\t"+strv.getStart()+"\t"+strv.getEnd()+"\t" +strv.getTandemRepeat().getPeriod() + "\t" + strv.getVar()+"\t"+strv.getVar2()+"\t" + strv.getConfidence() + "\n");
			
		}//while
		outS.close();
	}
}
