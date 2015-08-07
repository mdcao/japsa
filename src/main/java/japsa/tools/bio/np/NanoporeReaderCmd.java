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
 * 7 Aug 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsa.tools.bio.np;

import org.jfree.data.time.TimeTableXYDataset;

import japsa.seq.nanopore.NanoporeReaderStream;
import japsa.seq.nanopore.NanoporeReaderWindow;
import japsa.util.CommandLine;
import japsa.util.JapsaException;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(	
	scriptName = "jsa.np.f5reader", 
	scriptDesc = 
	"Extract Oxford Nanopore sequencing data from FAST5 files, perform an "
	+ "initial analysis of the date and stream them to realtime analysis pipelines"
	)
public class NanoporeReaderCmd extends CommandLine{	
	public NanoporeReaderCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addBoolean("GUI", false,"Run with a Graphical User Interface");
		addBoolean("realtime", false,"Run the program in real-time mode, i.e., keep waiting for new data from Metrichor agent");
		addString("folder", null,"The folder containing base-called reads");
		addBoolean("fail", false,"Get sequence reads from fail folder");
		addString("output", "-","Name of the output file, - for stdout");
		addString("streams", null,"Stream output to some servers, format \"IP:port,IP:port\" (no spaces)");
		addString("format", "fastq","Format of sequence reads (fastq or fasta)");
		addInt("minLength", 0,"Minimum read length");
		addBoolean("number", false,"Add a unique number to read name");
		addBoolean("stats", false,"Generate a report of read statistics");
		addBoolean("time", false,"Extract the sequencing time of each read -- only work with Metrichor > 1.12");		
			
		
		addStdHelp();		
	}
	
	public static void main(String[] args) throws OutOfMemoryError, Exception {		
		CommandLine cmdLine = new NanoporeReaderCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String output = cmdLine.getStringVal("output");		
		String folder = cmdLine.getStringVal("folder");		
		int minLength  = cmdLine.getIntVal("minLength");
		boolean stats  = cmdLine.getBooleanVal("stats");
		boolean number  = cmdLine.getBooleanVal("number");
		boolean time  = cmdLine.getBooleanVal("time");
		boolean GUI  = cmdLine.getBooleanVal("GUI");
		boolean realtime  = cmdLine.getBooleanVal("realtime");
		boolean fail  = cmdLine.getBooleanVal("fail");
		String format = cmdLine.getStringVal("format");		
		String streamServers = cmdLine.getStringVal("streams");

		//String pFolderName = cmdLine.getStringVal("pFolderName");
		//String f5list = cmdLine.getStringVal("f5list");
		//int interval = cmdLine.getIntVal("interval");//in second		
		//int age = cmdLine.getIntVal("age") * 1000;//in second
		int age = 20 * 1000;//cmdLine.getIntVal("age") * 1000;//in second
		int interval = 30;
		String pFolderName = null;

		if (!GUI && folder == null){// && f5list == null){
			Logging.exit("Download folder need to be specified", 1);
		}

		NanoporeReaderStream reader = new NanoporeReaderStream();

		reader.getTime = time;
		reader.stats = stats;
		reader.number = number;
		reader.minLength = minLength;

		reader.interval = interval;
		reader.age = age;

		//reader.f5List = f5list;
		reader.folder = folder;		
		reader.doFail = fail;
		reader.output = output;
		reader.format = format.toLowerCase();
		reader.realtime = realtime;
		reader.streamServers = streamServers;
		NanoporeReaderWindow mGUI = null;

		if (GUI){
			reader.realtime = true;
			System.setProperty("java.awt.headless", "false");
			reader.stats = true;//GUI implies stats
			reader.ready = false;//wait for the command from GUI

			TimeTableXYDataset dataset = new TimeTableXYDataset();
			mGUI = new NanoporeReaderWindow(reader,dataset);

			while (!reader.ready){
				Logging.info("NOT READY");
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {					
					e.printStackTrace();
				}			
			}
			Logging.info("GO");

			new Thread(mGUI).start();
		}else{
			String msg = reader.prepareIO();
			if (msg != null){
				Logging.exit(msg, 1);
			}
		}
		//reader need to wait until ready to go

		//reader.sos = SequenceOutputStream.makeOutputStream(reader.output);
		try{
			reader.readFastq(pFolderName);
		}catch (JapsaException e){
			System.err.println(e.getMessage());
			e.getStackTrace();
			if (mGUI != null)
				mGUI.interupt(e);
		}catch (Exception e){
			throw e;
		}finally{		
			reader.close();
		}

	}//main

}
