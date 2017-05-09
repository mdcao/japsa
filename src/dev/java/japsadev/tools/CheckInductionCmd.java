/*
 * Copyright (c) 2017  Minh Duc Cao (minhduc.cao@gmail.com).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the names of the institutions nor the names of the contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/*                           Revision History                                
 * 11/01/2012 - Minh Duc Cao: Revised 
 * 01/01/2013 - Minh Duc Cao, revised                                       
 ****************************************************************************/

package japsadev.tools;

import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.IOException;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.mrsainduction",
	scriptDesc = "Sample script description"
	)
public class CheckInductionCmd extends CommandLine{
	public CheckInductionCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addStdInputFile();
		//addBoolean("reverse",false,"Reverse sort order");
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new CheckInductionCmd();
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		
		String input = cmdLine.getStringVal("input");
		
		BufferedReader reader = SequenceReader.openFile(input);
		String line = "";
		while ((line = reader.readLine())!=null){
			
			if (line.startsWith("sbatch")){
				String [] toks = line.split(" ");
				System.out.println(toks[4]);
			}else if (line.startsWith("chr")){
				line = line.trim();
				String [] toks = line.split(":");
				System.out.print(line + "\t" + toks[2]);
				toks = toks[1].split("-");
				System.out.println(Double.parseDouble(toks[1]) - Double.parseDouble(toks[0]));
			}			
		}		
		reader.close();
	}	
}
/*RST*



 
  
  
  
  
  
  
  
  
*RST*/
  
