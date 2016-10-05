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
 * 28/05/2014 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsadev.tools.misc;


import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.HTSUtilities;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.dev.myn50", 
	scriptDesc = "compute n50 of a group"
	)
public class GetN50 extends CommandLine{
	//CommandLine cmdLine;
	public GetN50(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("input", null, "Name of the file",true);

		addStdHelp();
	}
	public static void main(String [] args) throws IOException, InterruptedException{
		GetN50 cmdTool = new GetN50 ();
		args = cmdTool.stdParseLine(args);

		/**********************************************************************/
		String input = cmdTool.getStringVal("input");
		processDB(input);
		//ArrayList<Sequence> seqs = SequenceReader.readAll(input, Alphabet.DNA());			
		//double n50 = HTSUtilities.n50(seqs);
		//System.out.println(n50 + "\t" + seqs.size());
	}


	private static void processDB(String dbFile) throws IOException, InterruptedException{		
		BufferedReader bf = SequenceReader.openFile(dbFile);
		ExecutorService executor = Executors.newFixedThreadPool(16);
				
		String line = "";
		int index = 0;
		while ( (line = bf.readLine())!=null){
			if (line.startsWith("#"))
				continue;
			
			index ++;
			line = line.trim();			
			executor.execute(new CompressSingle (index, line));
		}
		bf.close();		
		executor.shutdown();
		boolean finished = executor.awaitTermination(3, TimeUnit.DAYS);
		Logging.info("" + finished);
	}

	static class CompressSingle implements Runnable {
		int mIndex = 0;
		String mLine = "";
		public CompressSingle (int index, String line){
			mIndex = index;
			mLine = line;
		}
		public void run() {
			String [] toks = mLine.split("\t");			
			String fnaFile = toks[5];
			ArrayList<Sequence> seqs = null;
			try {
				seqs = SequenceReader.readAll(fnaFile, Alphabet.DNA());
			} catch (IOException e) {
				e.printStackTrace();
			}			
			double n50 = HTSUtilities.n50(seqs);
			synchronized(System.out){
				System.out.println(mLine + "\t" + n50 + "\t" + seqs.size() + "\t" + mIndex);
			}
		}
	}	
}
