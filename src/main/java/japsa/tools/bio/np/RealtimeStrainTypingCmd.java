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


import java.io.IOException;

import japsa.bio.np.RealtimeStrainTyping;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 *
 */

@Deployable(
	scriptName = "jsa.np.rtStrainTyping", 
	scriptDesc = "Realtime strain typing using Nanopore sequencing data",
	seeAlso = "jsa.np.f5reader, jsa.np.rtSpeciesTyping, jsa.np.rtResistGenes, jsa.util.streamServer, jsa.util.streamClient")
public class RealtimeStrainTypingCmd extends CommandLine{	
	public RealtimeStrainTypingCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("geneDB", null,  " Path to the gene database",true);		
		addString("bamFile", null,  "The bam file",true);

		addDouble("qual", 0,  "Minimum alignment quality");
		addBoolean("twodonly", false,  "Use only two dimentional reads");

		//addInt("top", 10,  "The number of top strains");		
		addInt("read", 50,  "Minimum number of reads between analyses");		
		addInt("time", 30,   "Minimum number of seconds between analyses");

		addString("output", "output.dat",  "Output file");

		addStdHelp();		
	} 

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException, InterruptedException{
		CommandLine cmdLine = new RealtimeStrainTypingCmd();		
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/

		String output = cmdLine.getStringVal("output");					
		String bamFile = cmdLine.getStringVal("bamFile");
		String geneDB = cmdLine.getStringVal("geneDB");
		int read       = cmdLine.getIntVal("read");
		int time       = cmdLine.getIntVal("time");		
		double qual      = cmdLine.getDoubleVal("qual");		
		boolean twoOnly      = cmdLine.getBooleanVal("twodonly");

		RealtimeStrainTyping paTyping = new RealtimeStrainTyping(read, time, geneDB,  output);
		paTyping.setMinQual(qual);	
		paTyping.setTwoOnly(twoOnly);
		paTyping.typing(bamFile);
	}
}


/*RST*
--------------------------------------------------------------------------------
*jsa.np.rtStrainTyping*: Bacterial strain typing with Oxford Nanopore sequencing
--------------------------------------------------------------------------------

*jsa.np.rtStrainTyping* strain types a bacterial sample from Oxford Nanopore
sequencing in real-time. It reads data in SAM/BAM format from a file or from
a stream and identifies the genes in the samples. Based on the patterns of 
gene presence, it makes an inference of the strain, together with the confidence
interval of 95%.

We provide the gene databases for three bacterial species  K. pneumoniae, 
E. coli and S. aureus on  http://genomicsresearch.org/public/researcher/npAnalysis/StrainTyping.tar.gz 
(or https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/StrainTyping.tar.gz).
Obtain them by::

   wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/StrainTyping.tar.gz
   tar zxvf StrainTyping.tar.gz

which will generate three folders for the three species.

<usage> 

~~~~~~~~~~~~~~
Usage examples
~~~~~~~~~~~~~~
If there is a sam/bam file of aligning the Nanopore sequencing to the gene 
database (ie, geneFam.fasta in one of the said databases), the program
can read from this file (note, this is not real-time analysis):
::

   jsa.np.rtStrainTyping -geneDB StrainTyping/Escherichia_coli/ -bamFile ecoli.bam -read 100000 -time 1000 -output output.dat
   
This program can read data from the output stream of an alignment program to
perform analysis in real-time. For example, one can create such a pipeline
to listen on port 3457
::

  jsa.util.streamServer -port 3457 \
  | bwa mem -t 2 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 -a StrainTyping/Escherichia_coli/geneFam.fasta - 2> /dev/null \
  | jsa.np.rtStrainTyping -bam - -geneDB StrainTyping/Escherichia_coli/ -read 0 -time 20 --out EcStrainTyping.dat 2>  kPStrainTyping.log
  
and streams data to this pipeline using npReader:
::

  jsa.np.f5reader -GUI -realtime -folder <DownloadFolder> -fail -output data.fastq -stream serverAddress:3457


*RST*/