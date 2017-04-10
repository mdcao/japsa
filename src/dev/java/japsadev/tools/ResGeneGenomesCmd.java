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
 * 07/09/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsadev.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import japsa.bio.amra.ResistanceGeneDB;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.tools.amra.Genomes2ResistanceGeneCmd;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;


/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.dev.refseq2res", 
	scriptDesc = "Extract resistance classes from sequences"
	)
public class ResGeneGenomesCmd  extends CommandLine{
	//CommandLine cmdLine;
	public ResGeneGenomesCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());						

		addString("db", "-", "name of db file");		
		addString("resDB", null, "Name of the resistance gene database",true);
		addString("output", "-", "name of output file");

		addDouble("identity", 0.9, "Minimum identity");
		addDouble("coverage", 0.8, "Minimum coverage of gene");

		addInt("thread",1, "Number of threads");
		
		addStdHelp();
	}

	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws Exception 
	 * @throws OutOfMemoryError 
	 */
	public static void main(String[] args) 
		throws	IOException, InterruptedException{

		CommandLine cmdLine = new ResGeneGenomesCmd();
		args = cmdLine.stdParseLine(args);		
		
		String db   = cmdLine.getStringVal("db");
		String resDBPath = cmdLine.getStringVal("resDB");
		
		String output   = cmdLine.getStringVal("output");

		gIdentity = cmdLine.getDoubleVal("identity");
		gCoverage = cmdLine.getDoubleVal("coverage");

		int thread  = cmdLine.getIntVal("thread");


		gResDB = new ResistanceGeneDB(resDBPath);
		gSos = SequenceOutputStream.makeOutputStream(output);
		
		gSos.print("#strainID\tstrainName\tstrainType\tclasses\n");
		
		if (!processDB(db, thread)){
			Logging.error("Job queue for too long");
		}
		gSos.close();
	}

	static double gIdentity, gCoverage;
	static SequenceOutputStream gSos;
	static ResistanceGeneDB gResDB;

	private static boolean processDB(String dbFile, int threadNumber) throws IOException, InterruptedException{		
		BufferedReader bf = SequenceReader.openFile(dbFile);		
		String line = "";

		ExecutorService executor = Executors.newFixedThreadPool(threadNumber);

		while ( (line = bf.readLine())!=null){
			if (line.startsWith("#"))
				continue;

			String [] toks = line.trim().split("\t");
			String strainID = toks[4];

			String fnaFile = toks[5];
			
			double n50 = Double.parseDouble(toks[7]);

			if (n50 < 100000){
				Logging.info(strainID + " Ignored because of low n50 " + n50);
				continue;//while
			}

			if (toks.length < 14){
				Logging.info(strainID + " Ignored because of malform " + toks.length);
				continue;//while
			}

			if (!toks[13].equals("0")){
				Logging.info(strainID + " Ignored because of not good ST " + toks[13]);
				continue;//while
			}

			String ST = toks[12];

			String organismName = toks[2];			
			if (!organismName.equals(toks[1] + " " + toks[3])){
				organismName += " " + toks[3];
			}			
			//organismName += "_ST" + ST;

			//Remove all the weird chars
			String strainName = organismName.replaceAll(" ", "_");
			strainName = strainName.replaceAll("/", "_");
			strainName = strainName.replaceAll("'", "_");
			strainName = strainName.replaceAll("\"", "_");
			strainName = strainName.replaceAll(";", "_");	
			strainName = strainName.replaceAll(":", "_");
			strainName = strainName.replaceAll("__*", "_");//Make sure no double hyphen
			/*****************************************************/				
			executor.execute(new ResGenome (strainID, fnaFile, ST, strainName));			
		}
		bf.close();
		executor.shutdown();
		boolean finished = executor.awaitTermination(7, TimeUnit.DAYS);	
		return finished;
	}

	static class ResGenome implements Runnable{

		String strainID;
		String strainName;
		String strainType;
		String fnaFile;		

		ResGenome(String strainID, String fnaFile, String strainType, String organismName){
			//mapClasses = shareMap;
			this.strainID = strainID;
			this.fnaFile = fnaFile;
			this.strainName = organismName;
			this.strainType = strainType;

		}

		/* (non-Javadoc)
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run() {
			try{
				ArrayList<Sequence> seqs = SequenceReader.readAll(fnaFile, Alphabet.DNA());
				HashSet<String> dgClasses = Genomes2ResistanceGeneCmd.blastn(seqs, gResDB, gIdentity, gCoverage);
								
				synchronized(gSos){
					gSos.print(strainID + "\t" + strainName + "\t" + strainType + "\t");
					for (String c:dgClasses){
						gSos.print(c + ", ");
					}
					gSos.println();					
				}

			}catch(IOException e){
				e.printStackTrace();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

	}

}

