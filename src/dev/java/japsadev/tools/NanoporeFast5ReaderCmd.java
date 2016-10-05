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
 * 24 Sep 2016 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsadev.tools;

import japsa.seq.SequenceOutputStream;
import japsa.seq.nanopore.Fast5DetailReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(	
	scriptName = "jsa.np.openfast5", 
	scriptDesc = "Extract data from a fast5 file. Still under development",
	seeAlso = "jsa.np.npreader"
	)
public class NanoporeFast5ReaderCmd extends CommandLine{	
	public NanoporeFast5ReaderCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("input", null,"Name of the fast5 file");
		addString("output", "-","Name of the output file, - for stdout");				
		addString("type", "fastq", 
				"Type of data to be extracted:" 
					+ "\nfastq: sequence read in fastq format"
					+ "\nevents: get events"
					+ "\nmodels: get models"
					+ "\nkeys:   list all keys"
				);		
		addStdHelp();		
	}

	public static void main(String[] args) throws Exception{		
		CommandLine cmdLine = new NanoporeFast5ReaderCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String input = cmdLine.getStringVal("input");
		String output = cmdLine.getStringVal("output");		
		String type   = cmdLine.getStringVal("type");
		try{
			//NanoporeReader npReader = new NanoporeReader(input);
			SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);
			//npReader.readData();
			//npReader.readKeys();
			//npReader.close();
			
			//NanoporeReader.readKeys(input, sos, true);
			
			//if (type.equals("fastq"))
			//	npReader.readFastq(fileList, minLength, sos, stats,number);
			//else if (type.equals("events"))
			//	readEvents(fileList, sos, stats);
			//else if (type.equals("models"))
			//	readModels(fileList, sos, stats);
			//else if (type.equals("keys"))
			//	readKeys(fileList, sos, stats);
			
			
			
			Fast5DetailReader f5Reader = new Fast5DetailReader(input);
			f5Reader.getRawEvent();						
			sos.print("Sampling rate " + f5Reader.getSamplingRate() + " " + f5Reader.getChannelNumber() + " " + f5Reader.getStartTime() );
			f5Reader.readData();
			
			
			f5Reader.close();
			sos.close();

			
			
			//npReader.readFastq(pFolderName);
		//	npReader.getEvents()
			
			//npReader.close();
		}catch (Exception e){
			throw e;
		}

	}//main
	
	
}
/******************************************************************

public static void main(String[] args) throws OutOfMemoryError, Exception {
	/*********************** Setting up script ****************************
	Deployable annotation = NanoporeReader.class.getAnnotation(Deployable.class);
	CommandLine cmdLine = new CommandLine("\nUsage: "
		+ annotation.scriptName() + " [options] f1.fast5 f2.fast5 ...",
		annotation.scriptDesc());

	cmdLine.addString("output", "-",
		"Name of the output file, -  for stdout");
	
	cmdLine.addInt("minLength", 0, 
		"Minimum sequence length");

	cmdLine.addBoolean("stats", false, "Compute statistics of reads");
	cmdLine.addBoolean("number", false, "Add a unique number to read name");
	cmdLine.addString("f5list",null, "File containing list of fast5 files, one file per line");		

	cmdLine.addStdHelp();		
	args = cmdLine.stdParseLine(args);
	/**********************************************************************

	
	String output = cmdLine.getStringVal("output");
	String f5list = cmdLine.getStringVal("f5list");
	int minLength  = cmdLine.getIntVal("minLength");
	boolean stats  = cmdLine.getBooleanVal("stats");
	boolean number  = cmdLine.getBooleanVal("number");

	ArrayList<String> fileList = new ArrayList<String>();
	if (f5list != null){
		BufferedReader bf = SequenceReader.openFile(f5list);
		String line;
		while ((line = bf.readLine()) != null){
			fileList.add(line.trim());
		}
		bf.close();
	}

	for (int i = 0; i < args.length; i++){
		fileList.add(args[i]);
	}
	

		//int maxLength = 0, minLength = Integer.MAX_VALUE;
}//main
/**********************************************************************/