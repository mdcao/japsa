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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
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
//static boolean reduceToSpecies = false;// whether to re-run after reducing db to identified species
static boolean buildConsensus = false;// this re-runs analysis and builds consensus;
static boolean saveSAM= true;
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
		addString("todo",null, "List of input files", false);
		addString("dbs",null, "databases to use in path",false);
		addDouble("removeLikelihoodThresh", 0.0, "what is the relative likelihood for removing species");
		addString("speciesFile",null, "species to restrict search",false);
		addBoolean("realtimeAnalysis", false, "whether to run analysis in realtime");
		addBoolean("alignedOnly", false, "whether to output only the aligned portion of a read in fasta file");
		//addDouble("removeLikelihoodThresh", 0.0, "likelihood proportion to remove");
		//addBoolean("reduceToSpecies",false, "whether to re-run alignment on reduced set of species, after E-M training and trimming");
	//	addString("reference", null, "Reference db if fastq is presented", false);
	//	addString("indexFile", null,  "indexFile ",true);
		addString("mm2Preset", "map-ont",  "mm2Preset ",false);
		addString("mm2_path", "minimap2",  "minimap2 path", false);
		addString("abpoa_path", "abpoa",  "abpoa path", false);

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
	
	public static Iterator<SAMRecord> [] getSamIteratorsFQ( File[] fastqFile, String readListSt, int maxReads,double q_thresh,	File refFile,Stack<File> keepBAM
		
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
	//
		
		String mm2_index = SequenceUtils.minimapIndex(refFile,  false,saveSAM);
		//sample_name.add (new File(files[0]));
		return  SequenceUtils.getSAMIteratorFromFastq(getStringVal(files), mm2_index, maxReads, readList, q_thresh, keepBAM);
	}
	
	private static String[] getStringVal(File[] files) {
		String[] res = new String[files.length];
		for(int i=0; i<res.length; i++)res[i] = files[i].getAbsolutePath();
		return res;
	}

	//public static boolean flipRefAndReadName = false;
	public static Iterator<SAMRecord>[]  getSamIteratorsBam(File[] bamFile ,String readListSt, int maxReads,double q_thresh,
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
		return iterators.toArray(new Iterator[0]);//SequenceUtils.getCombined(iterators.toArray(new Iterator[0]), false, true);
	}
	
	
	
	static void setParams(CommandLine cmdLine){
		RealtimeSpeciesTypingCmd.deleteUnmappedIntermediates = cmdLine.getBooleanVal("deleteUnmapped");
		SequenceUtils.mm2_threads= cmdLine.getIntVal("mm2_threads");
		SequenceUtils.mm2_mem = cmdLine.getStringVal("mm2_mem");
		SequenceUtils.mm2_path = cmdLine.getStringVal("mm2_path");
		SequenceUtils.mm2Preset = cmdLine.getStringVal("mm2Preset");
		SequenceUtils.mm2_splicing = null;//
		SequenceUtils.abpoa_path = cmdLine.getStringVal("abpoa_path");
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
		String todo = cmdLine.getStringVal("todo");
		
		//String [] inputFiles_;
		if(todo!=null){
			BufferedReader br = new BufferedReader(new FileReader(todo));
			String st = "";
			StringBuffer ifs = new StringBuffer();
			boolean first = true;
			boolean bam = false;
			while((st=br.readLine())!=null){
				if(st.length()>0 && !st.trim().startsWith("#")){
					if(!first)ifs.append(":");
					ifs.append(st);
					first = false;
					if(st.endsWith(".bam") || st.endsWith(".sam") || st.endsWith("sam.out")){
						bam=true;
					}else{
						if(bam) throw new RuntimeException("inconsistent");
					}
				}
			}
			String inp = ifs.toString();
			if(bam) bamFile = inp; else fastqFile=inp;
			br.close();
		}
		
	//	reduceToSpecies = cmdLine.getBooleanVal("reduceToSpecies");
		//if(RealtimeSpeciesTyping.removeLikelihoodThresh<=0.00001) reduceToSpecies=false;
		String resDB = cmdLine.getStringVal("resdb");
		if(bamFile==null && fastqFile==null) throw new RuntimeException("must define fastqFile or bam file");
		String dbPath = cmdLine.getStringVal("dbPath");
		String dbs =cmdLine.getStringVal("dbs");
		String readList = cmdLine.getStringVal("readList");
		final String speciesFile=cmdLine.getStringVal("speciesFile");
		String exclfile_ = cmdLine.getStringVal("excludeFile");

		buildConsensus = cmdLine.getBooleanVal("buildConsensus");
		//File excl = exclfile==null ? null : new File(exclfile);
		//File consensus = consensusFile==null ? null : new File(consensusFile);
		File[] fastqFiles = getFiles(null,fastqFile,":");
		File[] bamFiles = getFiles(null, bamFile,":");
		int num_sources = fastqFiles==null ? bamFiles.length : fastqFiles.length;
		List<File>[] unmapped_reads = new List[num_sources];
		List<File>[] out_fastq = new List[num_sources];

		File[] outdir = new File[num_sources];
		boolean all_equal  = outdir.length>1;
		for(int i=0; i<outdir.length; i++){
			File sample_namesk = bamFile!=null ?  bamFiles[i]: fastqFiles[i];
			 outdir[i] =  new File(resdir+"/"+sample_namesk.getName()+".jST");//.getParentFile().getName());// : new File(sample_namesk.getAbsolutePath()+"."+refDB.dbs+".jST");
			 if(!outdir[i].getName().equals(outdir[0].getName())) all_equal = false;
		}
		if(all_equal){
			for(int i=0; i<outdir.length; i++){
				File sample_namesk = bamFile!=null ?  bamFiles[i]: fastqFiles[i];
				 outdir[i] =  new File(resdir+"/"+sample_namesk.getParentFile().getName()+".jST");// : new File(sample_namesk.getAbsolutePath()+"."+refDB.dbs+".jST");
			}
		}
		for(int i=0; i<outdir.length; i++){outdir[i].mkdirs();}
		
		for(int i=0; i<num_sources; i++)  {
			unmapped_reads[i] = new ArrayList<File>();
			out_fastq[i] = new ArrayList<File>();
		}

		File outdirTop = null;
			ReferenceDB refDB = null;
			 File taxdmp = new File(dbPath+"/taxdump");
			 if(speciesFile!=null){
				 File modDB = new File("./db");
				 File newDB = new File(modDB, dbs+"_"+speciesFile);//+"."+last_m);
				 if(newDB.exists()){
					 refDB = new ReferenceDB(newDB, taxdmp, mkTree);
				 }
			}
			if(refDB==null && dbs!=null){
				File db_dir = new File(dbPath+"/"+dbs);
				refDB = new ReferenceDB(db_dir, new File(db_dir, "genomeDB.fna.gz"),mkTree);
		
				if(speciesFile!=null){
					refDB = refDB.update(new File(speciesFile));
				}
			}
			if(refDB==null && dbs==null){
				// blast
				File db_dir =  new File("./blastDB"); //bamFiles[i].getParentFile();
				db_dir.mkdir();
					File inp = new File(db_dir,"combined_index.fa.index");
					File inp1 = new File(db_dir,"combined_index.fa");
					ReferenceDB.combineIndex(bamFiles, ".fa.index", inp);
				refDB = new ReferenceDB(db_dir, taxdmp,inp1, mkTree, ".index",3, true);
				refDB.tree.modAll(0, num_sources);
			}
			
			boolean runAnalysis = true;
			boolean keepNames = true;
			File[] exclfile = null;
			File[] consensus_regions = null;
			File[] consensus_regions1 = null;
			File[] coverage_output = null;
			File[] reads_output =  null;

			Stack<File> bamOut = null;
			if(blast){
				reads_output =  RealtimeSpeciesTyping.newF(outdir,"reads.blast.txt.gz", false);
				RealtimeSpeciesTyping.start_thresh = 0.00;
				//refDB  = new ReferenceDB(refDB,0, "Viruses".split(":"), new int[] {2});
				fastqFiles = new File[bamFiles.length];
				for(int i=0; i<bamFiles.length; i++){
					fastqFiles[i] = new File(bamFiles[i].getParentFile(), "consensus_output.fa");
				}
				coverage_output =  RealtimeSpeciesTyping.newF(outdir,"coverage_blast.zip", false);

				File[][] outDs = speciesTyping(refDB,  outdir , readList,  bamFiles, fastqFiles, output,
						out_fastq,  unmapped_reads,  exclfile, consensus_regions, coverage_output, reads_output, runAnalysis, bamOut, keepNames);
				System.err.println(Arrays.asList(outDs));
				return;
			}
			
		
			
			//String consensusFile = cmdLine.getStringVal("consensusFile");
			//List<String> species = new ArrayList<String>();
			bamOut = buildConsensus ?new Stack<File>() : null; //bamOut saves the bam files if created
			File[][] outDs;
			exclfile = 			RealtimeSpeciesTyping.newF(outdir, exclfile_, false);
			consensus_regions = RealtimeSpeciesTyping.newF(outdir, "consensus_regions.txt", false);
			consensus_regions1 = RealtimeSpeciesTyping.newF(outdir, "consensus_regions1.txt", false);
			File[] consensus_output =  RealtimeSpeciesTyping.newF(outdir,"consensus_output.fa", false);
			coverage_output =  RealtimeSpeciesTyping.newF(outdir,"coverage.zip", false);

/*if(consensus_output[0].exists()){
	// third step
	speciesTyping(refDB, outdir, null, null, consensus_output, "output1.dat",
			null, null, null, null, null,  true, null, false);
}else if (consensus_regions[0].exists()){
	//this is second setp and it uses coverage.txt to generate consensus
	
outDs = speciesTyping(refDB, outdir , readList,  bamFiles, (File[]) null, output,
		(List[]) null, (List[]) null, exclfile, consensus_regions, coverage_output,  false, null,  false);
  
  System.err.println("consensus "+consensus);
}else{
				// this is the first step, and it generates coverage.txt
				outDs = speciesTyping(refDB, outdir , readList,  bamFiles, fastqFiles, output,
							out_fastq, unmapped_reads, exclfile, consensus_regions, coverage_output, true, bamOut, false);
				
		}*/
boolean makeConsensus=false;
if(consensus_output[0].exists()){
	//third step
	fastqFiles = consensus_output; 
	bamFiles = null;
	
	out_fastq = null;
	unmapped_reads = null;
	exclfile=null;
	coverage_output =  RealtimeSpeciesTyping.newF(outdir,"coverage1.zip", false);
	consensus_regions = RealtimeSpeciesTyping.newF(outdir, "consensus_regions1.txt", false);
	bamOut = null;
	output = "output1.dat";
	keepNames = false;
}else if(consensus_regions[0].exists()){
	//this is second step and it uses coverage.txt to generate consensus	
	fastqFiles = null;
	out_fastq = null; 
	unmapped_reads = null; 
   runAnalysis = false; 
   keepNames = false;
   makeConsensus=true;
   reads_output = RealtimeSpeciesTyping.newF(outdir,"reads.txt.gz", false);
}else{
	System.err.println("first step");
	// this if the first step
}
outDs = speciesTyping(refDB, outdir , readList,  bamFiles, fastqFiles, output,
		out_fastq, unmapped_reads, exclfile, consensus_regions, coverage_output, reads_output,  runAnalysis, bamOut, keepNames);	
if(makeConsensus){
	File consensus = SequenceUtils.makeConsensus(new File(outDs[0][0], "fastqs"),4, true);
}
System.err.println(outDs[1][0]);
System.err.println(outDs[0][0]);
System.err.println("h");			
			
			//}
		
		/*if(outdirTop!=null){
		KrakenTree overall = new  KrakenTree(outdirTop, "results.krkn");
		//overa.trim(1e-16);
		OutputStreamWriter osw = new OutputStreamWriter(new FileOutputStream(new File(outdirTop,"results_combined.krkn")));
		overall.print(osw, new String[]{NCBITree.count_tag,NCBITree.count_tag1}, new String[] {"%5.3g","%5.3g"}, true, false);
		osw.close();
		}*/
		
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

	public static File[][] speciesTyping(
			ReferenceDB refDB, 
			File[] outdir, String readList,
		 File [] bamFile, File[] fastqFile, String output,	
		 List<File>[] out_fastq , 
		 List<File>[] unmapped_reads, File[] exclude, File[] consensus_regions,File[] coverage_output,File[] reads_output,
		 boolean runAnalysis, Stack<File> keepBAM, boolean keepNames
			) throws IOException, InterruptedException{
		if(exclude==null) exclude = new File[outdir.length];
		if(consensus_regions==null) consensus_regions = new File[outdir.length];
		int num_sources = bamFile!=null ? bamFile.length : fastqFile.length;
		File[][] result = new File[num_sources][];
		if(fastqFile!=null) {
			System.err.println("running on "+Arrays.asList(fastqFile));
		}
			List<SamReader> readers =  new ArrayList<SamReader>();
			Iterator<SAMRecord>[] samIter =
					bamFile!=null ? 	RealtimeSpeciesTypingCmd.getSamIteratorsBam(bamFile,  readList, maxReads, q_thresh, readers,  refDB.refFile, flipName) : 
						RealtimeSpeciesTypingCmd.getSamIteratorsFQ( fastqFile, readList, maxReads, q_thresh, refDB.refFile, keepBAM);
		
			String[] src_names =getSourceNames(bamFile!=null ? bamFile: fastqFile);
				
			//if(resdir==null) throw new RuntimeException("resdir is null");
		//	File outdir_new  = resdir;
			//resdir.mkdir();
			
			/*if(resdir!=null){
				outdir_new =  new File(resdir,refDB.dbs);
				outdir_new.mkdir();
			}*/
			
			
			
			//	SamReader samReader = readers.size()>0 ? readers.get(k) : null;
				RealtimeSpeciesTyping paTyping =
						new RealtimeSpeciesTyping(refDB,	exclude,consensus_regions,coverage_output, reads_output, 
								 output,outdir,  refDB.refFile, unmapped_reads!=null, keepNames, src_names);
				
				paTyping.setMinQual(qual);
				paTyping.setTwoOnly(twoOnly);	
				paTyping.setFilter(filter);
				try{
				paTyping.typing(samIter, number, time, runAnalysis);
				}catch(InterruptedException exc){
					exc.printStackTrace();
				}
				for(int i=0; i<readers.size(); i++) readers.get(i).close();
				readList = null;
				if(out_fastq!=null && runAnalysis) paTyping.getOutfiles(out_fastq);
				for(int src_index=0; src_index<num_sources; src_index++){
				if(paTyping.fqw_unmapped!=null){
					paTyping.fqw_unmapped[src_index].getOutFile(unmapped_reads[src_index]);
					System.err.println("unmapped reads "+Arrays.asList(unmapped_reads[src_index]));
				}
				}
				//files[k] = paTyping.unmapped_reads;  // unmapped reads taken forward to next database
	//	}dbs
				for(int i=0; i<num_sources; i++){
					result[i] = new File[] {outdir[i], paTyping.consensus_file_out[i]};
				}
				refDB.printKrakenFormat(resdir, paTyping.samples(), keepNames ? fastqFile : null);
				return new File[][] {outdir, paTyping.consensus_file_out};
		//paTyping.typing(bamFile, number, time);		
	}

	private static String[] getSourceNames(File[]  bamFile) {
		int num_sources = bamFile.length;
		String[] src_names = new String[num_sources];
		boolean all_equal = true;
		for(int i=0; i<src_names.length; i++){
			src_names[i] = bamFile[i].getName();
			if(i>0 && src_names[i].equals(src_names[0])){
				all_equal = false;
			}
		}
		if(!all_equal && num_sources>1){
			for(int i=0; i<src_names.length; i++){
				src_names[i] = bamFile[i].getParentFile().getName();
			}
		}
		return src_names;
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
