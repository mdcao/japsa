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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.bio.np.RealtimeResistanceGene;
import japsa.bio.np.RealtimeSpeciesTyping;
import japsa.bio.phylo.CSSProcessCommand;
import japsa.bio.phylo.NCBITree;
import japsa.bio.phylo.Trie;
import japsa.tools.seq.CachedOutput;
import japsa.tools.seq.SequenceUtils;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

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

	static double q_thresh=7; 
	static double qual=1;
	
	static String filter;
	
	static boolean twoOnly=false;
	static int maxReads=Integer.MAX_VALUE;
	static int number ;
	static int time=30;
	static File resdir = new File("japsa_species_typing");
	public RealtimeSpeciesTypingCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("resdir", "japsa_species_typing", "Results directory");
		addString("output", "output.dat",  "Output file, - for standard output");		
		addString("bamFile", null,  "The bam file",false);	
		addString("fastqFile", null, "Fastq file", false);
		addString("dbPath",null, "path to databases",false);
		addString("resdb",null, "Resistance database",false);
		addString("dbs",null, "databases to use in path",false);
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
		addDouble("fail_thresh", 7.0,  "median phred quality of read", false);
		addInt("mm2_threads", 4, "threads for mm2", false);
		addDouble("qual", 1,  "Minimum alignment quality");
		addBoolean("twodonly", false,  "Use only two dimentional reads");
		addDouble("alpha", 0.05, "Paramater alpha from multinomialCI");
		addInt("minCount", 5, "Mininum number of mapped reads for a species to be considered");
		addString("filter", "", "List of species (separated by semicolon) to excluded from typing");

		addInt("read", 50,  "Minimum number of reads between analyses");		
		addInt("time", 30,   "Minimum number of seconds between analyses");

		addBoolean("web", false, "Whether to use Web visualization.");
		addBoolean("log", false, "Whether to write mapping details to species2reads.map.");
		addBoolean("merge", false, "whether to merge results from multiple bam into single output file");
		
		addBoolean("writeSep" , false, "whether to output fastq for each species detected");
		addBoolean("plasmidOnly", false,  "If writeSep is true, then whether to only output plasmid mapping reads");

		addBoolean("writeUnmapped", false, "whether to output fastq of unmapped reads");
		addString("speciesToIgnore","GRCh38:GRCh38,"," species not to extract sequence for",false);

		addStdHelp();		
	} 
	static class MatchFilter implements FilenameFilter{
		final Pattern p;
		MatchFilter(String st){
			this.p = Pattern.compile(st);
		}
		@Override
		public boolean accept(File dir, String name) {
			return p.matcher(name).find();
		}
	}
	
	
	public static void  getSamIterators(String[] bamFile, String[] fastqFile, String readListSt, int maxReads,double q_thresh,
			List<String> sample_name, List<Iterator<SAMRecord>> iterators,List<SamReader> samReaders,File refFile
			) throws IOException, FileNotFoundException{
		boolean bam = fastqFile==null;
		String[] files = bam ? bamFile : fastqFile;
		
		Collection<String> readList=SequenceUtils.getReadList(readListSt, true);
		
		if(bamFile!=null && bamFile.length==1 && !(new File(bamFile[0])).exists()) {
			files = (new File(".")).list(new MatchFilter(bamFile[0]));	
		}
		if(fastqFile!=null && fastqFile.length==1 &&  !(new File(fastqFile[0])).exists()) {
			files = (new File(".")).list(new MatchFilter(fastqFile[0]));	
		}
		if(files.length==0) {
			System.err.println(bam ? bamFile[0]  : fastqFile[0]);
			throw new RuntimeException("no files match input request");
		}
	//	String[] sample_name = new String[files.length];
		for(int i=0; i<files.length; i++){
			File filek = new File(files[i]);//.replace(".gz","");//.substring(0, files[i].lastIndexOf('.'));
			sample_name.add(filek.getName()	);
			Iterator<SAMRecord> samIter = null;
			SamReader samReader= null;
		if(bam){
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			InputStream bamInputStream;
			if("-".equals(filek)){
				bamInputStream = System.in;
			}else{
				bamInputStream =	new FileInputStream(filek);
			}
				samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamInputStream));
			//	sam_it1 = 
				samReaders.add(samReader);
				samIter = SequenceUtils.getFilteredIterator(samReader.iterator(), readList, maxReads, q_thresh);
		}else{
			boolean saveSeqs=true;
			try{
			String mm2_index = SequenceUtils.minimapIndex(refFile,  false,saveSeqs);
			samIter = SequenceUtils.getSAMIteratorFromFastq(filek.getAbsolutePath(), mm2_index, maxReads, readList, q_thresh);
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		iterators.add(samIter);
		}
	}
	
	static class ReferenceDB{
		// this class encapsulates everything required for a referenceDB
		File refFile, modDB,speciesIndex;
		
		String dbs;
		NCBITree tree;
		public ReferenceDB(String dbPath, String dbs, String speciesFile)  throws IOException{
			this.dbs  = dbs;
			File dbdir = new File(dbPath+"/"+dbs);
			 refFile = new File(dbdir, "genomeDB.fna.gz");
			 File taxaDir = new File(dbPath+"/taxdump");
			speciesIndex=Trie.getIndexFile(taxaDir, refFile);
			String treef = CSSProcessCommand.getTree(taxaDir,dbdir, speciesIndex, false).getAbsolutePath();
			boolean useTaxaAsSlug=false;
			String treef_mod = treef+".mod";
			tree = new NCBITree(new File(treef), useTaxaAsSlug);
			tree.addSpeciesIndex(speciesIndex);
			//tree.print(new File(treef_mod)); // this prints out to tree
	//		if(true)System.exit(0);
			// treef = dbPath+"/"+dbs+"/"+ "commontree.txt.css.mod";
			 modDB = new File("./db");
			 if(speciesFile!=null){
				 this.update(speciesFile, treef);
			 }
			
		}
		public void update(String speciesFile, String treef) {
			
			Collection<String> species = SequenceUtils.getReadList(speciesFile,false);
			if(species!=null){
				 modDB.mkdir();
				File speciesF = new File(speciesFile);
				long last_m = speciesF.lastModified();
				File refFileOut = new File(modDB,dbs+"_"+speciesFile+"."+last_m+".fna.gz");
				File indexFileOut = new File(modDB,dbs+"_"+speciesFile+"_speciesIndex."+last_m+".txt");
				if(!refFileOut.exists()){
					SequenceUtils.mkdb(refFile, treef, speciesIndex.getAbsolutePath(), species, refFileOut.getAbsolutePath(), indexFileOut.getAbsolutePath());
				}
				refFile = refFileOut;//.getAbsolutePath();
				speciesIndex  = indexFileOut;
			}
		}
		
	}
	
	
	
	static void setParams(CommandLine cmdLine){
		
		SequenceUtils.mm2_threads= cmdLine.getIntVal("mm2_threads");
		SequenceUtils.mm2_mem = cmdLine.getStringVal("mm2_mem");
		SequenceUtils.mm2_path = cmdLine.getStringVal("mm2_path");
		SequenceUtils.mm2Preset = cmdLine.getStringVal("mm2Preset");
		SequenceUtils.mm2_splicing = null;//
		SequenceUtils.secondary = false;
		RealtimeSpeciesTypingCmd.q_thresh = cmdLine.getDoubleVal("fail_thresh");
		RealtimeSpeciesTypingCmd.filter = cmdLine.getStringVal("filter");
		RealtimeSpeciesTypingCmd.maxReads = cmdLine.getIntVal("maxReads");
		RealtimeSpeciesTypingCmd.number       = cmdLine.getIntVal("read");
		RealtimeSpeciesTypingCmd.time       = cmdLine.getIntVal("time");		
		RealtimeSpeciesTypingCmd.qual      = cmdLine.getDoubleVal("qual");				
		RealtimeSpeciesTypingCmd.twoOnly      = cmdLine.getBooleanVal("twodonly");
		RealtimeSpeciesTyping.JSON = cmdLine.getBooleanVal("web");
		RealtimeSpeciesTyping.OUTSEQ = cmdLine.getBooleanVal("log");
		RealtimeSpeciesTyping.ALPHA = cmdLine.getDoubleVal("alpha");
		RealtimeSpeciesTyping.MIN_READS_COUNT = cmdLine.getIntVal("minCount");
		RealtimeSpeciesTyping.writeSep = cmdLine.getBooleanVal( "writeSep");
		RealtimeSpeciesTyping.plasmidOnly = cmdLine.getBooleanVal("plasmidOnly");
		RealtimeSpeciesTyping.writeUnmapped = cmdLine.getBooleanVal("writeUnmapped");
		RealtimeSpeciesTyping.speciesToIgnore = Arrays.asList(cmdLine.getStringVal("speciesToIgnore").split(":"));
		CachedOutput.MIN_READ_COUNT = RealtimeSpeciesTyping.MIN_READS_COUNT;
		RealtimeSpeciesTyping.realtimeAnalysis = cmdLine.getBooleanVal("realtimeAnalysis");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLine cmdLine = new RealtimeSpeciesTypingCmd();		
		args = cmdLine.stdParseLine(args);		
		setParams(cmdLine);
		resdir = new File(cmdLine.getStringVal("resdir"));
		resdir.mkdirs();
	    String output    = cmdLine.getStringVal("output");
		String bamFile   = cmdLine.getStringVal("bamFile");		
		String fastqFile = cmdLine.getStringVal("fastqFile");
		String resDB = cmdLine.getStringVal("resdb");
		if(resDB!=null){
			RealtimeSpeciesTyping.plasmidOnly = false;
			RealtimeSpeciesTyping.writeSep = true;
		}
		if(bamFile==null && fastqFile==null) throw new RuntimeException("must define fastqFile or bam file");
		String dbPath = cmdLine.getStringVal("dbPath");
		String dbs = cmdLine.getStringVal("dbs");//.split(":");
		String readList = cmdLine.getStringVal("readList");
		String speciesFile=cmdLine.getStringVal("species");
		ReferenceDB refDB = new ReferenceDB(dbPath, dbs, speciesFile);
		List<String> out_fastq = new ArrayList<String>();
		speciesTyping(refDB, resdir, readList, bamFile==null ? null : bamFile.split(":"), 
										fastqFile==null ? null : fastqFile.split(":"), output, out_fastq);
	
		if(resDB!=null && out_fastq.size()>0){
			SequenceUtils.secondary = true;
			CachedOutput.MIN_READ_COUNT=2;
			List<String> outfiles = new ArrayList<String>();
			File outdir = new File(".");
			RealtimeResistanceGene.writeSep=false; // no need to write fastq files
			RealtimeResistanceGeneCmd.resistanceTyping(new File(resDB),resdir,  null,
					out_fastq.toArray(new String[0]), readList, outdir, output, outfiles);
		}
	}
	public static void speciesTyping(ReferenceDB refDB, File resdir, String readList,
		 String [] bamFile, String[] fastqFile, String output,	List<String> out_fastq 
			) throws IOException{
			List<String> sample_names = new ArrayList<String>();	
			List<Iterator<SAMRecord>> iterators =  new ArrayList<Iterator<SAMRecord>>();
			List<SamReader> readers =  new ArrayList<SamReader>();
			RealtimeSpeciesTypingCmd.getSamIterators(bamFile, fastqFile, readList, maxReads, q_thresh, sample_names,iterators, readers,  refDB.refFile);
			if(false){
				// this merges multiple iterators into one
				//Iterator<SAMRecord> it = SequenceUtils.getCombined(iterators.toArray(new Iterator[0]), false, true);
			}
			File outdir_new = 
			 new File(resdir,refDB.dbs);
			outdir_new.mkdir();
			for(int k=0; k<iterators.size(); k++){ // do multiple samples sequentially , could consider doing in parallel later
				File outdir = new File(outdir_new+"/"+sample_names.get(k));
				outdir.mkdirs();
				Iterator<SAMRecord> samIter = iterators.get(k);
				SamReader samReader = readers.size()>0 ? readers.get(k) : null;
				RealtimeSpeciesTyping paTyping =
						new RealtimeSpeciesTyping(refDB.speciesIndex,	refDB.tree, output,outdir,  refDB.refFile);
				paTyping.setMinQual(qual);
				paTyping.setTwoOnly(twoOnly);	
				paTyping.setFilter(filter);
				try{
				paTyping.typing(samIter, number, time);
				}catch(InterruptedException exc){
					exc.printStackTrace();
				}
				if(samReader!=null) samReader.close();
				readList = null;
				paTyping.getOutfiles(out_fastq);
				//files[k] = paTyping.unmapped_reads;  // unmapped reads taken forward to next database
			}
	//	}dbs
		
		//paTyping.typing(bamFile, number, time);		
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
