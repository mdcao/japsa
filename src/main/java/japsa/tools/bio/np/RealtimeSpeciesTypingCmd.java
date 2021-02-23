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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.bio.np.RealtimeSpeciesTyping;
import japsa.tools.seq.SequenceUtils;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import scala.Int;

/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.np.rtSpeciesTyping", 
	scriptDesc = "Realtime species typing using Nanopore Sequencing data",
	seeAlso = "jsa.np.npreader, jsa.np.rtStrainTyping, jsa.np.rtResistGenes, jsa.util.streamServer, jsa.util.streamClient"
	)
public class RealtimeSpeciesTypingCmd extends CommandLine {

	public RealtimeSpeciesTypingCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("output", "output.dat",  "Output file, - for standard output");		
		addString("bamFile", null,  "The bam file",false);	
		addString("fastqFile", null, "Fastq file", false);
		addString("dbPath",null, "path to databases",true);
		addString("dbs",null, "databases to use in path",true);
		addString("species",null, "species to restrict search",false);
		addBoolean("realtimeAnalysis", false, "whether to run analysis in realtime");
	//	addString("reference", null, "Reference db if fastq is presented", false);
	//	addString("indexFile", null,  "indexFile ",true);
		addString("mm2Preset", null,  "mm2Preset ",false);
		addString("mm2_path", "/sw/minimap2/current/minimap2",  "minimap2 path", false);
		addString("readList", null,  "file with reads to include", false);
		addInt("maxReads",Integer.MAX_VALUE, "max reads to process", false );
		
		//addString("mm2Preset", "splice",  "preset for minimap2", false);
	//	addBoolean("writeBed", false, "whether to write bed",false);
		//addString("mm2Preset", "map-ont",  "preset for minimap2", false);
//		addString("mm2_memory", "4g",  "minimap2 memory", false);
		long mem = (Runtime.getRuntime().maxMemory()-1000000000);
		addString("mm2_memory", mem+"",  "minimap2 memory", false);

		addInt("mm2_threads", 4, "threads for mm2", false);
		
		addDouble("qual", 1,  "Minimum alignment quality");
		addBoolean("twodonly", false,  "Use only two dimentional reads");
		
		addDouble("alpha", 0.05, "Paramater alpha from multinomialCI");
		addInt("minCount", 0, "Mininum number of mapped reads for a species to be considered");
		addString("filter", "", "List of species (separated by semicolon) to excluded from typing");

		addInt("read", 50,  "Minimum number of reads between analyses");		
		addInt("time", 30,   "Minimum number of seconds between analyses");

		addBoolean("web", false, "Whether to use Web visualization.");
		addBoolean("log", false, "Whether to write mapping details to species2reads.map.");
		addStdHelp();		
	} 
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		
		CommandLine cmdLine = new RealtimeSpeciesTypingCmd();		
		args = cmdLine.stdParseLine(args);		

		/**********************************************************************/

		String output    = cmdLine.getStringVal("output");
		String bamFile   = cmdLine.getStringVal("bamFile");		
		String[] fastqFile = cmdLine.getStringVal("fastqFile").split(":");
		if(bamFile==null && fastqFile==null) throw new RuntimeException("must define fastqFile or bam file");
		
		SequenceUtils.mm2_threads= cmdLine.getIntVal("mm2_threads");
		SequenceUtils.mm2_mem = cmdLine.getStringVal("mm2_mem");
		SequenceUtils.mm2_path = cmdLine.getStringVal("mm2_path");
		SequenceUtils.mm2Preset = cmdLine.getStringVal("mm2Preset");
		SequenceUtils.mm2_splicing = null;//
		SequenceUtils.secondary = false;
		//SequenceUtils.mm2Preset = cmdLine.getStringVal("mm2Preset");
		String dbPath = cmdLine.getStringVal("dbPath");
		String[] dbs = cmdLine.getStringVal("dbs").split(":");
	
	//	String reference = cmdLine.getStringVal("reference");
		//String indexFile = cmdLine.getStringVal("indexFile");
		
		
		String filter = cmdLine.getStringVal("filter");
		int maxReads = cmdLine.getIntVal("maxReads");
		int number       = cmdLine.getIntVal("read");
		int time       = cmdLine.getIntVal("time");		
		double qual      = cmdLine.getDoubleVal("qual");				
		boolean twoOnly      = cmdLine.getBooleanVal("twodonly");
		RealtimeSpeciesTyping.JSON = cmdLine.getBooleanVal("web");
		RealtimeSpeciesTyping.OUTSEQ = cmdLine.getBooleanVal("log");
		RealtimeSpeciesTyping.ALPHA = cmdLine.getDoubleVal("alpha");
		RealtimeSpeciesTyping.MIN_READS_COUNT = cmdLine.getIntVal("minCount");
		String[] sample_name = new String[fastqFile.length];
		for(int i=0; i<sample_name.length; i++){
			sample_name[i] = 	(fastqFile!=null  ? fastqFile[i] : bamFile).replace(".gz", "");
			String[] str = sample_name[i].split("\\.");
			sample_name[i] = str[str.length-2];
		}
		for(int k=0; k<sample_name.length; k++){ // do multiple samples sequentially , could consider doing in parallel later
			Collection<String> readList=getReadList(cmdLine.getStringVal("readList"), true);
			String speciesFile=cmdLine.getStringVal("species");
			File modDB = new File("./db");modDB.mkdir();
			Collection<String> species = getReadList(speciesFile,false);
					
			String outdir_new = "./";
			for(int i=0; i<dbs.length; i++){
				
				outdir_new = outdir_new+"_"+dbs[i];
				File outdir = new File(outdir_new+"/"+sample_name[k]);
				outdir.mkdirs();
				String refFile = dbPath+"/"+dbs[i]+"/genomeDB.fna.gz";
				String indexFile=dbPath+"/"+dbs[i]+"/speciesIndex";
				String treef = "commontree.txt.css.mod";
				treef = dbPath+"/"+dbs[i]+"/"+treef;
				String speciesIndex = dbPath+"/"+dbs[i]+"/speciesIndex";
				if(species != null){
					File speciesF = new File(speciesFile);
					long last_m = speciesF.lastModified();
					File refFileOut = new File(modDB,dbs[i]+"_"+speciesFile+"."+last_m+".fna.gz");
					File indexFileOut = new File(modDB,dbs[i]+"_"+speciesFile+"_speciesIndex."+last_m+".txt");
					if(!refFileOut.exists()){
						SequenceUtils.mkdb(refFile, treef, speciesIndex, species, refFileOut.getAbsolutePath(), indexFileOut.getAbsolutePath());
					}
					refFile = refFileOut.getAbsolutePath();
					indexFile  = indexFileOut.getAbsolutePath();
				}
				
				RealtimeSpeciesTyping paTyping =
						new RealtimeSpeciesTyping(indexFile,
												treef, 
													output,outdir,  refFile);
				paTyping.setMinQual(qual);
				paTyping.setTwoOnly(twoOnly);	
				paTyping.setFilter(filter);
				Iterator<SAMRecord> samIter = null;
				SamReader samReader= null;
			
				if(bamFile!=null){
				InputStream bamInputStream;
		
				if ("-".equals(bamFile))
					bamInputStream = System.in;
				else
					bamInputStream = new FileInputStream(bamFile);
					SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
					samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamInputStream));
					samIter= samReader.iterator();
				}else{
					boolean saveSeqs=false;
					String mm2_index = SequenceUtils.minimapIndex(new File(refFile),  false,saveSeqs);
					//String mm2_index = refFile;
					samIter = SequenceUtils.getSAMIteratorFromFastq(fastqFile[k], mm2_index, maxReads, readList);
				}
				paTyping.typing(samIter, number, time);
				if(samReader!=null) samReader.close();
				readList = null;
				fastqFile[k] = paTyping.unmapped_reads;  // unmapped reads taken forward to next database
			}
		}
		
		//paTyping.typing(bamFile, number, time);		
	}
	private static Collection<String> getReadList(String readList, boolean split) {
		if(readList==null || readList=="null") return null;
		List<String> reads = new ArrayList<String>();
		try{
		
		InputStream is = new FileInputStream(new File(readList));
		 if(readList.endsWith(".gz")) is = new GZIPInputStream(is);
		 BufferedReader br = new BufferedReader(new InputStreamReader(is));
		 String st = "";
		 while((st = br.readLine())!=null){
			 
			 reads.add(split ? st.split("\\s+")[0] : st);
		 }
		}catch(IOException exc){
			exc.printStackTrace();
		}
		return reads;
	}
}

/*RST*
----------------------------------------------------------------------------------
*jsa.np.rtSpeciesTyping*: Bacterial species typing with Oxford Nanopore sequencing
----------------------------------------------------------------------------------

*jsa.np.rtSpeciesTyping* identify proportions of species from a DNA sample 
using Oxford Nanopore sequencing in real-time. It reads data in SAM/BAM format
of the alignments of sequence reads to a collection of species genomes.

We provide a genome collection of nearly 1500 bacterial species
on  http://data.genomicsresearch.org/Projects/npAnalysis/.
Refer to the documentation at https://github.com/mdcao/npAnalysis/ for more 
details.
 
<usage> 

~~~~~~~~~~~~~~
Usage examples
~~~~~~~~~~~~~~

If there is a sam/bam file of aligning the Nanopore sequencing to the genome 
collection, the program can read from this
::

   jsa.np.rtSpeciesTyping -bam alignment.sam -index SpeciesTyping/Bacteria/speciesIndex --read 50 -time 60 -out speciesTypingResults.out
   
   
This program can read data from the output stream of an alignment program to
perform analysis in real-time. For example, one can create such a pipeline
to listen on port 3456
::

  jsa.util.streamServer -port 3456 \
  | bwa mem -t 10 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 SpeciesTyping/Bacteria/genomeDB.fasta - 2> /dev/null \
  | jsa.np.rtSpeciesTyping -bam - -index SpeciesTyping/Bacteria/speciesIndex --read 50 -time 60 -out speciesTypingResults.out 2>  speciesTypingResults.log &
  
  
and streams data to this pipeline using npReader:
::

  jsa.np.npreader -GUI -realtime -folder <DownloadFolder> -fail -output data.fastq -stream serverAddress:3456


*RST*/
