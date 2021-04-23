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
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.bio.np.AllRecords;
import japsa.bio.np.RealtimeResistanceGene;
import japsa.bio.np.RealtimeSpeciesTyping;
import japsa.bio.phylo.KrakenTree;
import japsa.bio.phylo.NCBITree;
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
static boolean reduceToSpecies = false;// whether to re-run after reducing db to identified species
static boolean buildConsensus = false;// this re-runs analysis and builds consensus;
	static double q_thresh=7; 
	static double qual=1;
	static boolean deleteUnmappedIntermediates = true;
	static String filter;
	
	public static boolean mkTree = true;
	public static boolean flipName = false;
	public static boolean blast = false;
	
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
		addInt("max_EM_iterations",5000, "How many max iterations of E-M");
		addString("output", "output.dat",  "Output file, - for standard output");		
		addString("bamFile", null,  "The bam file",false);
		addString("consensusFile", null,  "consensus file",false);
		addBoolean("deleteUnmapped",true, "whehter to delete the unmapped reads",false);
		addString("fastqFile", null, "Fastq file", false);
		addString("dbPath",null, "path to databases",false);
		addString("resdb",null, "Resistance database",false);
		addString("dbs",null, "databases to use in path",false);
		addDouble("removeLikelihoodThresh", 0.0, "what is the relative likelihood for removing species");
		addString("speciesFile",null, "species to restrict search",false);
		addBoolean("realtimeAnalysis", false, "whether to run analysis in realtime");
		addBoolean("alignedOnly", false, "whether to output only the aligned portion of a read in fasta file");
		//addDouble("removeLikelihoodThresh", 0.0, "likelihood proportion to remove");
		addBoolean("reduceToSpecies",false, "whether to re-run alignment on reduced set of species, after E-M training and trimming");
	//	addString("reference", null, "Reference db if fastq is presented", false);
	//	addString("indexFile", null,  "indexFile ",true);
		addString("mm2Preset", "map-ont",  "mm2Preset ",false);
		addString("mm2_path", "/sw/minimap2/current/minimap2",  "minimap2 path", false);
		addString("abpoa_path", "/sw/abpoa/current/abpoa",  "abpoa path", false);

		addString("readList", null,  "file with reads to include", false);
		addInt("maxReads",Integer.MAX_VALUE, "max reads to process", false );
		addInt("minCoverage",2, "minimum coverage from consensus file", false );
		addInt("minLength",500, "minimum length from consensus file", false );
		//addString("mm2Preset", "splice",  "preset for minimap2", false);
	//	addBoolean("writeBed", false, "whether to write bed",false);
		//addString("mm2Preset", "map-ont",  "preset for minimap2", false);
//		addString("mm2_memory", "4g",  "minimap2 memory", false);
		long mem = (Runtime.getRuntime().maxMemory()-1000000000);
		addString("mm2_memory", mem+"",  "minimap2 memory", false);
		addString("excludeFile",null,  "file of regions to exclude", false);
		addBoolean("buildConsensus",false,  "builds the consensus by running twice", false);
		addDouble("fail_thresh", 7.0,  "median phred quality of read", false);
		addInt("mm2_threads", 4, "threads for mm2", false);
		addDouble("qual", 0,  "Minimum alignment quality");
		addBoolean("twodonly", false,  "Use only two dimentional reads");
		addDouble("alpha", 0.05, "Paramater alpha from multinomialCI");
		addInt("minCount", 5, "Mininum number of mapped reads for a species to be considered");
		addString("filter", "", "List of species (separated by semicolon) to excluded from typing");

		addInt("read", 50,  "Minimum number of reads between analyses");		
		addInt("time", 30,   "Minimum number of seconds between analyses");

		addBoolean("web", false, "Whether to use Web visualization.");
		addBoolean("log", false, "Whether to write mapping details to species2reads.map.");
		//addBoolean("merge", false, "whether to merge results from multiple bam into single output file");
		
		addString("writeSep" , null, "strings to match for what to write fastq file out, which can be colon separated, e.g. plasmid:phage or all");

		//addBoolean("writeUnmapped", false, "whether to output fastq of unmapped reads");
		addString("speciesToIgnore","GRCh38:GRCh38,"," species not to extract sequence for",false);
		addBoolean("flipName",false, "whether to exchange reference and read name in the SAM record (necessary for blast");
		addBoolean("mkTree", true, "whether to generate a tree for analysis");
		addBoolean("blast", false,"whether this is blast analysis");
		addString("tag",null, "which tag to use for calculating the likleihood of a read in the  EM calculation.  Default is using mapQ scores");
		addStdHelp();		
	} 
	static class MatchFilter implements FilenameFilter{
		final Pattern p,p1;
		
		MatchFilter(String st, String suffix){
			this.p = Pattern.compile(st);
			this.p1 = Pattern.compile(suffix);
		}
		@Override
		public boolean accept(File dir, String name) {
			boolean res =  p.matcher(name).find() && p1.matcher(name).find();
			return res;
		}
	}
	
	public static Iterator<SAMRecord>  getSamIteratorsFQ( File[] fastqFile, String readListSt, int maxReads,double q_thresh,	File refFile,Stack<File> keepBAM
		
			) throws IOException, FileNotFoundException, InterruptedException{
		File[] files =  fastqFile;
		Collection<String> readList=SequenceUtils.getReadList(readListSt, true);
		if(fastqFile.length==1 &&  !(fastqFile[0].exists())) {
			File parent = fastqFile[0].getParentFile();
			if(parent==null) parent = new File(".");
			files =parent.listFiles(new MatchFilter(fastqFile[0].getName(), "fq|fastq"));	
		}
	//	String[] fastqFile_new = fastqFile;
		
		if(files.length==0) {
			throw new RuntimeException("no files match input request");
		}
		boolean saveSeqs=true;
		
		String mm2_index = SequenceUtils.minimapIndex(refFile,  false,saveSeqs);
		//sample_name.add (new File(files[0]));
		return  SequenceUtils.getSAMIteratorFromFastq(getStringVal(files), mm2_index, maxReads, readList, q_thresh, keepBAM);
	}
	
	private static String[] getStringVal(File[] files) {
		String[] res = new String[files.length];
		for(int i=0; i<res.length; i++)res[i] = files[i].getAbsolutePath();
		return res;
	}

	//public static boolean flipRefAndReadName = false;
	public static Iterator<SAMRecord>  getSamIteratorsBam(File[] bamFile ,String readListSt, int maxReads,double q_thresh,
			List<SamReader> samReaders,File refFile, boolean flipRefAndReadName
			) throws IOException, FileNotFoundException{
		 System.err.println(Arrays.asList(bamFile));
		Collection<String> readList=SequenceUtils.getReadList(readListSt, true);
		File[] files = bamFile;
		if(bamFile.length==1 && !(files[0].exists())) {
			File parent = files[0].getParentFile();
			if(parent==null) parent = new File(".");
			files = parent.listFiles(new MatchFilter(bamFile[0].getName(), "bam|sam"));	
		}
		
		if(files.length==0) {
		
		//	else System.err.println(Arrays.asList(fastqFile));
			throw new RuntimeException("no files match input request");
		}
	//	String[] sample_name = new String[files.length];
		//sample_name.add (new File(files[0]));
		List<Iterator> iterators = new ArrayList<Iterator>();
		for(int i=0; i<files.length; i++){
			File filek =files[i];//.replace(".gz","");//.substring(0, files[i].lastIndexOf('.'));
			Iterator<SAMRecord> samIter = null;
			SamReader samReader= null;
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			InputStream bamInputStream;
			if("-".equals(filek)){
				bamInputStream = System.in;
			}else{
				bamInputStream =	new FileInputStream(filek);
			}
				samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamInputStream));
				samReaders.add(samReader);
				samIter = SequenceUtils.getFilteredIterator(samReader.iterator(), readList, maxReads, q_thresh, flipRefAndReadName);
				iterators.add(samIter);
		}
		return SequenceUtils.getCombined(iterators.toArray(new Iterator[0]), false, true);
	}
	
	
	
	static void setParams(CommandLine cmdLine){
		RealtimeSpeciesTypingCmd.deleteUnmappedIntermediates = cmdLine.getBooleanVal("deleteUnmapped");
		SequenceUtils.mm2_threads= cmdLine.getIntVal("mm2_threads");
		SequenceUtils.mm2_mem = cmdLine.getStringVal("mm2_mem");
		SequenceUtils.mm2_path = cmdLine.getStringVal("mm2_path");
		SequenceUtils.mm2Preset = cmdLine.getStringVal("mm2Preset");
		SequenceUtils.mm2_splicing = null;//
		SequenceUtils.apboa_path = cmdLine.getStringVal("abpoa_path");
		SequenceUtils.secondary = true;
		RealtimeSpeciesTyping.removeLikelihoodThresh= cmdLine.getDoubleVal("removeLikelihoodThresh");
		RealtimeSpeciesTyping.EMiterations = cmdLine.getIntVal("max_EM_iterations");
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
		RealtimeSpeciesTyping.writeSep =cmdLine.getStringVal( "writeSep")==null ? null :  Pattern.compile(cmdLine.getStringVal( "writeSep"));
		RealtimeSpeciesTyping.alignedOnly     = cmdLine.getBooleanVal("alignedOnly");//"qual");				
		RealtimeSpeciesTyping.minCoverage = cmdLine.getIntVal("minCoverage");
		RealtimeSpeciesTyping.minlength = cmdLine.getIntVal("minLength");

		//boolean match = RealtimeSpeciesTyping.writeSep.matcher("abcde").find();
	//System.err.println(match);
	//if(true) System.exit(0);;
		//RealtimeSpeciesTyping.writeUnmapped = cmdLine.getBooleanVal("writeUnmapped");
		RealtimeSpeciesTyping.speciesToIgnore = Arrays.asList(cmdLine.getStringVal("speciesToIgnore").split(":"));
		CachedOutput.MIN_READ_COUNT = RealtimeSpeciesTyping.MIN_READS_COUNT;
		RealtimeSpeciesTyping.realtimeAnalysis = cmdLine.getBooleanVal("realtimeAnalysis");
		
		AllRecords.tag = cmdLine.getStringVal("tag");
		blast = cmdLine.getBooleanVal("blast");
		mkTree = cmdLine.getBooleanVal("mkTree");
		flipName = cmdLine.getBooleanVal("flipName");
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
		RealtimeSpeciesTyping.removeLikelihoodThresh = cmdLine.getDoubleVal("removeLikelihoodThresh");
	    String output    = cmdLine.getStringVal("output");
		String bamFile   = cmdLine.getStringVal("bamFile");		
		String fastqFile = cmdLine.getStringVal("fastqFile");
		reduceToSpecies = cmdLine.getBooleanVal("reduceToSpecies");
		if(RealtimeSpeciesTyping.removeLikelihoodThresh<=0.00001) reduceToSpecies=false;
		String resDB = cmdLine.getStringVal("resdb");
		if(bamFile==null && fastqFile==null) throw new RuntimeException("must define fastqFile or bam file");
		String dbPath = cmdLine.getStringVal("dbPath");
		String[] dbs =cmdLine.getStringVal("dbs")==null ? new String[] {null}:  cmdLine.getStringVal("dbs").split(":");
		String readList = cmdLine.getStringVal("readList");
		final String speciesFile=cmdLine.getStringVal("speciesFile");
		List<File> out_fastq = new ArrayList<File>();
		String exclfile = cmdLine.getStringVal("excludeFile");
		buildConsensus = cmdLine.getBooleanVal("buildConsensus");
		//File excl = exclfile==null ? null : new File(exclfile);
		//File consensus = consensusFile==null ? null : new File(consensusFile);
		File[] fastqFiles = getFiles(null,fastqFile,":");
		File[] bamFiles = getFiles(null, bamFile,":");
		List<File> unmapped_reads = dbs.length>1 ? new ArrayList<File>(): null;
		File outdirTop = null;
		
		
		inner: for(int i=0; i<dbs.length; i++){
		//	System.err.println(dbs[i]);
			ReferenceDB refDB = null;
			 File taxdmp = new File(dbPath+"/taxdump");
			 if(speciesFile!=null){
				 File modDB = new File("./db");
				 File newDB = new File(modDB, dbs[i]+"_"+speciesFile);//+"."+last_m);
				 if(newDB.exists()){
					 refDB = new ReferenceDB(newDB, taxdmp, mkTree);
				 }
			}
			if(refDB==null && dbs[i]!=null){
				File db_dir = new File(dbPath+"/"+dbs[i]);
				refDB = new ReferenceDB(db_dir, new File(db_dir, "genomeDB.fna.gz"),mkTree);
		
				if(speciesFile!=null){
					refDB = refDB.update(new File(speciesFile));
				}
			}
			if(refDB==null && dbs[i]==null){
				File db_dir = bamFiles[i].getParentFile();
				if(db_dir==null) db_dir = new File("./");
				refDB = new ReferenceDB(db_dir, taxdmp,new File(bamFiles[i].getAbsolutePath()+".fa"), mkTree, ".index",3, true);
			}
			
			
			if(blast){
				File[] outDs = speciesTyping(refDB, i==0 ? resdir : null, readList,  bamFiles, fastqFiles, output,
						out_fastq, i==dbs.length-1 ? null : unmapped_reads, exclfile, null, null, true, null, null);
				return;
			}
			
			
			String consensusFile = cmdLine.getStringVal("consensusFile");
			List<String> species = new ArrayList<String>();
			Stack<File> bamOut = buildConsensus ?new Stack<File>() : null;
		//	ReferenceDB refDB, File resdir, String readList,
		//	 File [] bamFile, File[] fastqFile, String output,	List<String> out_fastq , 
		//	 List<File> unmapped_reads, String exclude, String consensus,
		//	 List<String> species, boolean runAnalysis, Stack<File> keepBAM, File outdir
			File[] outDs = speciesTyping(refDB, i==0 ? resdir : null, readList,  bamFiles, fastqFiles, output,
							out_fastq, i==dbs.length-1 ? null : unmapped_reads, exclfile, consensusFile, species, true, bamOut, null);
			if(buildConsensus && consensusFile==null && !dbs[i].equals("Human")  ){
				String consensusFile1 = consensusFile==null ? outDs[1].getAbsolutePath(): consensusFile;
				System.err.println("rerunning to build consensus");
				File[] bamFiles1 = bamFiles;
				File bamO = bamOut.size()>0 ? bamOut.pop() : null;
				//System.err.println(bamO);
				if(bamFiles1==null)bamFiles1 = new File[] {bamO};
				outDs = speciesTyping(refDB, i==0 ? resdir : null, readList,  bamFiles1, null, output,
						null, null, exclfile, consensusFile1, null, false, null, outDs[0]);
				if(bamO!=null) bamO.delete();
				File consensus = SequenceUtils.makeConsensus(new File(outDs[0], "fastqs"),4, true);
				speciesTyping(refDB, null, null, null, new File[] {consensus}, "output.dat",
						null, null, null, null, null, true, null, new File(consensus.getAbsolutePath()+".jST"));
			
			}
			if(speciesFile ==null && species.size()>0 &&  reduceToSpecies && !dbs[i].equals("Human") ){
				File specFile1 = new File(resdir, dbs[i]+"."+System.currentTimeMillis()+".txt");
				PrintWriter pw = new PrintWriter(new FileWriter(specFile1));
				for(int k=0; k<species.size(); k++){
					pw.println(species.get(k));
				}
				pw.close();
				List<String> species1 = new ArrayList<String>();
				refDB = refDB.update(specFile1);
			//	consensusFile = cmdLine.getStringVal("consensusFile");
				File[] fastqFiles1 = fastqFiles;
				System.err.println("running on subset");
				outDs = speciesTyping(refDB, i==0 ? resdir : null, readList,  null, fastqFiles1, output,
						out_fastq,  null , exclfile, consensusFile, species1, true, bamOut, null);
				System.err.println(species1.size());
				if(buildConsensus && consensusFile==null){
					String consensusFile1 = consensusFile==null ? outDs[1].getAbsolutePath(): consensusFile;
					System.err.println("rerunning to build consensus");
					File[] bamFiles1 = bamFiles;
					File bamO = bamOut.pop();
					if(bamFiles1==null)bamFiles1 = new File[] {bamO};
					outDs = speciesTyping(refDB, i==0 ? resdir : null, readList,  bamFiles1, null, output,
							null,  null , exclfile, consensusFile1,null, false, null, outDs[0]);
					bamO.delete();
					File consensus = SequenceUtils.makeConsensus(new File(outDs[0], "fastqs"),2, true);
				}
			}
			if(outdirTop==null && !dbs[i].equals("Human")) outdirTop = outDs[0];
			bamFiles = null;
			if(unmapped_reads==null) break inner;
			fastqFiles = unmapped_reads.toArray(new File[0]);
			if(deleteUnmappedIntermediates){
			if(i!=dbs.length-1){
				//dont delete the final unmapped reads
				for(int j=0; j<unmapped_reads.size(); j++){
							((unmapped_reads.get(j))).deleteOnExit();
				}
			}
			}
			unmapped_reads.clear();
			if(fastqFiles.length==0) break inner;
		}
		if(outdirTop!=null){
		KrakenTree overall = new  KrakenTree(outdirTop, "results.krkn");
		//overa.trim(1e-16);
		overall.print(new File(outdirTop,"results_combined.krkn"), new String[]{NCBITree.count_tag,NCBITree.count_tag1}, new String[] {"%d","%d"}, true);
		}
		if(resDB!=null && out_fastq.size()>0){
			SequenceUtils.secondary = true;
			CachedOutput.MIN_READ_COUNT=2;
			List<File> outfiles = new ArrayList<File>();
			File outdir = new File(".");
			RealtimeResistanceGene.writeSep=false; // no need to write fastq files
			//File resDB, File resdir, File[] bamFile, 
//			File[] fastqFile,String readList, File outdir, String output, List<File> outfiles
			RealtimeResistanceGeneCmd.resistanceTyping(new File(resDB),resdir,  null,
					out_fastq.toArray(new File[0]), readList, outdir, output, outfiles);
		}
	}
	
	
public static File[]  getFiles(File base, String str, String sp) {
		if(str==null) return null;
		String[] split = sp==null ? new String[] {str} : str.split(sp);
		File[] f = new File[split.length];
		for(int i=0; i<split.length; i++){
			f[i] = base==null ? new File(split[i]) : new File(base,split[i]);
		}
		return f;
	}

	public static File[] speciesTyping(ReferenceDB refDB, File resdir, String readList,
		 File [] bamFile, File[] fastqFile, String output,	List<File> out_fastq , 
		 List<File> unmapped_reads, String exclude, String consensus,
		 List<String> species, boolean runAnalysis, Stack<File> keepBAM, File outdir
			) throws IOException, InterruptedException{
		
		if(fastqFile!=null) {
		
			System.err.println("running on "+Arrays.asList(fastqFile));
		}
			List<SamReader> readers =  new ArrayList<SamReader>();
			Iterator<SAMRecord> samIter= 
					bamFile!=null ? 	RealtimeSpeciesTypingCmd.getSamIteratorsBam(bamFile,  readList, maxReads, q_thresh, readers,  refDB.refFile, flipName) : 
						RealtimeSpeciesTypingCmd.getSamIteratorsFQ( fastqFile, readList, maxReads, q_thresh, refDB.refFile, keepBAM);
			File outdir_new  = null;
			if(resdir!=null){
				outdir_new =  new File(resdir,refDB.dbs);
				outdir_new.mkdir();
			}
			File sample_namesk = bamFile!=null ?  bamFile[0]: fastqFile[0];
				if(outdir==null)  outdir = outdir_new!=null ?  new File(outdir_new+"/"+sample_namesk.getName()) : new File(sample_namesk.getAbsolutePath()+"."+refDB.dbs+".jST");
				outdir.mkdirs();
			
			//	SamReader samReader = readers.size()>0 ? readers.get(k) : null;
				RealtimeSpeciesTyping paTyping =
						new RealtimeSpeciesTyping(refDB,	exclude,consensus,
								 output,outdir,  refDB.refFile, unmapped_reads!=null);
				
				paTyping.setMinQual(qual);
				paTyping.setTwoOnly(twoOnly);	
				paTyping.setFilter(filter);
				try{
				paTyping.typing(samIter, number, time, species, runAnalysis);
				}catch(InterruptedException exc){
					exc.printStackTrace();
				}
				for(int i=0; i<readers.size(); i++) readers.get(i).close();
				readList = null;
				if(out_fastq!=null && runAnalysis) paTyping.getOutfiles(out_fastq);
				if(paTyping.fqw_unmapped!=null){
					paTyping.fqw_unmapped.getOutFile(unmapped_reads);
					System.err.println("unmapped reads "+Arrays.asList(unmapped_reads));
				}
				//files[k] = paTyping.unmapped_reads;  // unmapped reads taken forward to next database
	//	}dbs
		
		return new File[] {outdir,paTyping.consensus_file_out} ;
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
