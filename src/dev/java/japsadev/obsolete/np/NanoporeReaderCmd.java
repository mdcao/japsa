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
package japsadev.obsolete.np;

import org.jfree.data.time.TimeTableXYDataset;
import japsa.util.CommandLine;
import japsa.util.JapsaException;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */
@Deployable(	
	scriptName = "jsa.np.npreader", 
	scriptDesc = "Extract and stream Oxford Nanopore sequencing data in real-time",
	seeAlso = "jsa.np.filter, jsa.util.streamServer, jsa.util.streamClient,jsa.np.rtSpeciesTyping, jsa.np.rtStrainTyping, jsa.np.rtResistGenes"
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
		//addString("group", "","Group of base-called to be extracted ()");
		addInt("minLength", 1,"Minimum read length");
		addBoolean("number", false,"Add a unique number to read name");
		addBoolean("stats", false,"Generate a report of read statistics");
		//addBoolean("time", false,"Extract the sequencing time of each read -- only work with Metrichor > 1.12");		


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
		//boolean time  = cmdLine.getBooleanVal("time");
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

		//reader.getTime = time;
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


/*RST*
-------------------------------------------------------------------------
*npReader*: real-time conversion and analysis of Nanopore sequencing data
-------------------------------------------------------------------------

*npReader* (jsa.np.npreader) is a program that extracts Oxford Nanopore
sequencing data from FAST5 files, performs an initial analysis of the date and
streams them to real-time analysis pipelines. These pipelines can run on the
same computer or on computing clouds/high performance clusters.

npReader is included in the `Japsa package <http://mdcao.github.io/japsa/>`_.
It requires
`JAVA HDF5 INTERFACE (JHI5) library <https://www.hdfgroup.org/products/java/JNI/jhi5/index.html>`_
to be installed prior to setting up Japsa. Details of installation as follows:

**On Windows/Mac**

1. Download and install HDF-View from
https://www.hdfgroup.org/products/java/release/download.html.
Note the folder that the JHI library is installed, e.g.,
*C:\\Program Files\\HDF_Group\\HDFView\\2.11.0\\lib*

2. Follow the instructions to install Japsa on
http://japsa.readthedocs.org/en/latest/install.html.
Upon prompting for "Path to HDF library", enter the above path.

**On Linux**

You can either install the JHI5 library by downloading the software from
*https://www.hdfgroup.org/products/java/JNI/jhi5/index.html* or from your
Linux distribution software repository, such as::

   sudo apt-get install libjhdf5-jni

The library is typically installed to *#/usr/lib/jni*. Enter this path when
prompted for "Path to HDF library" during installation of Japsa.

HDF-View (https://www.hdfgroup.org/products/java/release/download.html) also 
contains the neccessary library. Please install HDF-2.10.1 instead of the 
latest version.

<usage>

~~~~~~~~~~~~~~
Usage examples
~~~~~~~~~~~~~~

A summary of npReader usage can be obtained by invoking the --help option::

   jsa.np.npreader --help

The simplest way to run *npReader* in GUI mode is by typing::

   jsa.np.npreader -GUI -realtime

and specify various options in the GUI. All of these options can be specified
from the command line::

   jsa.np.npreader -GUI -realtime -folder c:\Downloads\ -fail -output myrun.fastq --minLength 200 --stats

npReader can run natively on a Windows laptop that runs the Metrichor agent. It
can stream sequence data to multiple analysis pipelines on the same computer
and/or on high performance clusters and computing clouds.

Start several analysis pipelines on some remote machines. Such a pipeline can
be to count how many reads aligned to chromosomes A and B::

   jsa.util.streamServer --port 3456 | \
   bwa mem -t 8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 index - | \
   awk -F "\t" 'BEGIN{A=0;B=0;N++} NF>4 \
       {if ($3=="chrA") A++; if ($3=="chrB") B++; \
        if (NR %100==0) \
          {print "At " NR " reads, " A " aligned to chr A; " B " aligned to chr B"} \
       }'  

In this pipeline, the *jsa.util.streamServer* program receives stream data
from *npReader* and forwards to *bwa*, which aligns the data to a reference
and in turn streams the alignment in sam format to the awk program to perform
a simple analysis of counting reads aligned to chrA and chrB.

The Japsa package contains several real-time analysis (jsa.np.speciesTyping,
jsa.np.geneStrainTyping, jsa.np.resistGenes). They can be used to set up
analysis pipelines, such as::

   jsa.util.streamServer --port 3457 \
   bwa mem -t 8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 index - | \   
   jsa.np.speciesTyping  -bam - --index speciesIndex -output output.dat

Once these pipelines are ready, npReader can start streaming data off the
MinION and the Metrichor agent to these pipelines::

   jsa.np.npreader -realtime -folder c:\Downloads\ -fail -output myrun.fastq \
      --minLength 200 --streams server1IP:3456,server2IP:3457

One can run *npReader* on a computing cloud if the download folder (containing
base-called data) can be mounted to the cloud. In such case, npReader can
direct stream data to the pipelines without the need of
*jsa.util.streamServer*::

   jsa.np.npreader -realtime -folder c:\Downloads\ -fail -output - | \
   bwa mem -t 8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 index - | \
   jsa.np.speciesTyping  -bam - --index speciesIndex -output output.dat

Japsa also provides *jsa.np.filter*, a tool to bin sequence data in groups of
the user's liking. Like any other streamline tools, jsa.np.filter can run
behind *jsa.util.streamServer* on a remote machine, or can get data directly
from npReader via pipe::

   jsa.np.npreader -realtime -folder c:\Downloads\ -fail -output - | \
   jsa.np.filter -input - -lenMin 2000 --qualMin 10 -output goodreads.fq

One can also use *tee* to group data into different bins *in real-time* with
*jsa.np.filter*::

   jsa.np.npreader -realtime -folder c:\Downloads\ -fail -output - | \   
   tee >(jsa.np.filter -input - -lenMax 2000 -output 0k2k.fq) \ 
   >(jsa.np.filter -lenMin 2000 -lenMax 4000 -input - -output 2k4k.fq) \ 
   >(jsa.np.filter -lenMin 4000 -lenMax 6000 -input - -output 4k6k.fq) \
   >(jsa.np.filter -lenMin 6000 -input - -output 6k.fq) \
   > all.fq

These bins can also be piped/streamed to different analysis pipelines as above.

*RST*/
