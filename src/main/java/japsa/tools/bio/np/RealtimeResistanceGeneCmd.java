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
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
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
	scriptName = "jsa.np.rtResistGenes", 
	scriptDesc = "Realtime identification of antibiotic resistance genes from Nanopore sequencing",
	seeAlso = "jsa.np.npreader, jsa.np.rtSpeciesTyping, jsa.np.rtStrainTyping, jsa.util.streamServer, jsa.util.streamClient"
	)
public class RealtimeResistanceGeneCmd extends CommandLine{	
	public RealtimeResistanceGeneCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("writeABX" , null, "strings to match for what to write fastq file out, which can be colon separated, e.g. fosfomycin|vancomycin");
		addInt("minCountResistance", 1, "Mininum number of mapped reads for a species to be considered (for species typing step)");
		addInt("minCountSpecies", 1, "Mininum number of mapped reads for a species to be considered (for species typing step)");
		addString("abpoa_path", "/sw/abpoa/current/abpoa",  "abpoa path", false);

		addString("output", "output.dat",  "Output file");
		addString("bamFile", null,  "The bam file");
		addString("fastqFile", null,  "fastq input");
		addBoolean("runKAlign", false, "whether to run msa to get high confidence calls");
		addDouble("score", 0.0001,  "The alignment score threshold");
		addString("msa", "kalign",
			"Name of the msa method, support poa, kalign, muscle and clustalo");
		addString("tmp", "_tmpt",  "Temporary folder");				
		addString("resDB", null,  "Path to resistance database", true);
		addString("dbs", null,  "For subsequent species typing", false);
		addString("dbPath",null, "path to databases",false);
		addString("resdir", "japsa_resistance_typing", "Results directory");
		addDouble("qual", 0,  "Minimum alignment quality");
		addBoolean("twodonly", false,  "Use only two dimentional reads");				
		addInt("read", 50,  "Minimum number of reads between analyses");		
		addInt("time", 1800,   "Minimum number of seconds between analyses");
		addBoolean("log", false, "Whether to write mapping details to genes2reads.map.");
		addBoolean("buildConsensus", true, "whether to try to build species consensus during species typing ");

		addInt("thread", 1,   "Number of threads to run");

		addString("mm2Preset", null,  "mm2Preset ",false);
		addString("mm2_path", "/sw/minimap2/current/minimap2",  "minimap2 path", false);
		addString("readList", null,  "file with reads to include", false);
		addInt("maxReads",Integer.MAX_VALUE, "max reads to process", false );
		long mem = (Runtime.getRuntime().maxMemory()-1000000000);
		addString("mm2_memory", mem+"",  "minimap2 memory", false);
		addDouble("fail_thresh", 7.0,  "median phred quality of read", false);
		addInt("mm2_threads", 4, "threads for mm2", false);
		addDouble("qual", 1,  "Minimum alignment quality");
		
		
		addStdHelp();
	} 
	
	static double scoreThreshold = 0.0010;
	static int readPeriod = 50;
	static int time = 100;
	static int thread = 1;
	static boolean twodonly = false;
   static  int maxReads = Integer.MAX_VALUE;
   static String msa = "kalign";
   static double q_thresh=7;
   static String tmp="_tmpt";
	static File resdir = new File("japsa_resistance_typing");
public static Pattern writeABX = null;
	public static void main(String[] args) throws IOException, InterruptedException{
		CommandLine cmdLine = new RealtimeResistanceGeneCmd();
		args = cmdLine.stdParseLine(args);		
		scoreThreshold = cmdLine.getDoubleVal("score");		
		readPeriod = cmdLine.getIntVal("read");
		writeABX = cmdLine.getStringVal("writeABX") == null ? null : Pattern.compile(cmdLine.getStringVal("writeABX"));
		time = cmdLine.getIntVal("time");
		thread = cmdLine.getIntVal("thread");
		twodonly = cmdLine.getBooleanVal("twodonly");
		msa = cmdLine.getStringVal("msa");
		maxReads = cmdLine.getIntVal("maxReads");
		q_thresh = cmdLine.getDoubleVal("fail_thresh");
		resdir = new File(cmdLine.getStringVal("resdir"));
		resdir.mkdirs();
		SequenceUtils.mm2_threads= cmdLine.getIntVal("mm2_threads");
		SequenceUtils.mm2_mem = cmdLine.getStringVal("mm2_mem");
		SequenceUtils.mm2_path = cmdLine.getStringVal("mm2_path");
		SequenceUtils.mm2Preset = cmdLine.getStringVal("mm2Preset");
		SequenceUtils.mm2_splicing = null;//
		SequenceUtils.secondary = true;
		SequenceUtils.apboa_path  = cmdLine.getStringVal("abpoa_path");
		CachedOutput.MIN_READ_COUNT=cmdLine.getIntVal("minCountResistance"); // this for detection of abx
		
		RealtimeSpeciesTyping.MIN_READS_COUNT = cmdLine.getIntVal("minCountSpecies");
		RealtimeSpeciesTyping.mincount=2;
		RealtimeSpeciesTyping.minCoverage=2;
		RealtimeSpeciesTyping.minlength=300;
		RealtimeResistanceGene.OUTSEQ = cmdLine.getBooleanVal("log");
		RealtimeResistanceGene.runKAlign = cmdLine.getBooleanVal("runKAlign");
		tmp = cmdLine.getStringVal("tmp");		

		
		String output = cmdLine.getStringVal("output");
		String bamFile = cmdLine.getStringVal("bam");		
		File resDB = new File(cmdLine.getStringVal("resDB"));
		String fastqFile = cmdLine.getStringVal("fastqFile");
		String readList = cmdLine.getStringVal("readList");
		
		if(bamFile==null && fastqFile==null) throw new RuntimeException("must define fastqFile or bam file");

		File outdir = new File("./");
		List<File> outfiles = new ArrayList<File>();
		resistanceTyping(resDB,resdir,  RealtimeSpeciesTypingCmd.getFiles(null, bamFile, ":"), RealtimeSpeciesTypingCmd.getFiles(null, fastqFile, ":"),
				readList, outdir, output, outfiles);
		//now do species typing on the resistance genes;
		String dbPath =  cmdLine.getStringVal("dbPath");
		String[] dbs = cmdLine.getStringVal("dbs") == null ? null : cmdLine.getStringVal("dbs").split(":");
		File[] fastqFiles = outfiles.toArray(new File[0]);
		String excl = null;// can add in excl file here
		File currDir = new File(".");
		boolean buildConsensus =cmdLine.getBooleanVal("buildConsensus");
		String consensusFile= null;
		if(dbPath!=null && dbs!=null && outfiles.size()>0){
			
			CachedOutput.MIN_READ_COUNT=RealtimeSpeciesTyping.MIN_READS_COUNT;
			RealtimeSpeciesTyping.writeSep = Pattern.compile("[a-z]");
			SequenceUtils.secondary = false;
			if(fastqFiles.length>0){
			outer: for(int k=fastqFiles.length-1; k>=0; k--){
				
				List<File> unmapped_reads = dbs.length>1 ? new ArrayList<File>(): null;
				File[] fqFiles = new  File[] {fastqFiles[k]};
				File outdirTop = null;
				inner: for(int i=0; i<dbs.length; i++){
					System.err.println("species typing for "+fastqFiles[k]+" "+dbs[i]+"...");
					ReferenceDB refDB = new ReferenceDB(new File(dbPath+"/"+dbs[i]));
					List<File> species_output_files = new ArrayList<File>();
				//	if(fastqFiles.length==0) break;
					List<String> species = new ArrayList<String>();
					if(fqFiles.length==0) break inner;
					Stack<File> bamOut = buildConsensus ?new Stack<File>() : null;

					File[] outD = RealtimeSpeciesTypingCmd.speciesTyping(refDB, null, null, null, fqFiles,  "output.dat", species_output_files,
							i==dbs.length-1 ? null : unmapped_reads, excl, null, species, true,  bamOut, null);
					System.err.println("..done");
					if(buildConsensus && consensusFile==null && !dbs[i].equals("Human") ){
						String consensusFile1 = consensusFile==null ? outD[1].getAbsolutePath(): consensusFile;
						File bamO = bamOut.pop();
						File[] bamFiles1 = new File[] {bamO};
						System.err.println("re-typing for consensus "+bamO.getName()+" "+dbs[i]);

						outD = RealtimeSpeciesTypingCmd.speciesTyping(refDB, i==0 ? resdir : null, readList,  bamFiles1, null, "output.dat",
								null, null, excl, consensusFile1, null, false, null, outD[0]);
						System.err.println("..done");

						bamO.delete();
						System.err.println("using abpoa to make consensus "+SequenceUtils.apboa_path);
						File consensus = SequenceUtils.makeConsensus(new File(outD[0], "fastqs"),4, true);
						System.err.println("done  "+consensus.getName());
					//	System.err.println("re typing for consensus "+consensus.getAbsolutePath()+" "+dbs[i]);
//						RealtimeSpeciesTypingCmd.speciesTyping(refDB, null, null, null, new File[] {consensus}, "output.dat",
	//						null, null, null, null, null, true, null, new File(consensus.getAbsolutePath()+".jST"));
					
					}
					
					if(outdirTop==null && !dbs[i].equals("Human"))  outdirTop = outD[0];
					if(unmapped_reads==null) break inner;
					fqFiles = unmapped_reads.toArray(new File[0]);
					for(int j=0; j<unmapped_reads.size(); j++){
						unmapped_reads.get(j).deleteOnExit();
					}
					unmapped_reads.clear();
				}
				if(outdirTop!=null){
					KrakenTree overall = new  KrakenTree(outdirTop, "results.krkn");
					//overa.trim(1e-16);
					overall.print(new File(outdirTop,"results_combined.krkn"), new String[]{NCBITree.count_tag,NCBITree.count_tag1}, new String[] {"%d","%d"}, true);
					}
			}
			}
		}
	}

	public static void resistanceTyping(File resDB, File resdir, File[] bamFile, 
			File[] fastqFile,String readList, File outdir, String output, List<File> outfiles)  throws IOException, InterruptedException{
	
		List<SamReader> readers =  new ArrayList<SamReader>();
		Iterator<SAMRecord> samIter= 
				bamFile!=null ? 	RealtimeSpeciesTypingCmd.getSamIteratorsBam( bamFile,  readList, maxReads, q_thresh, readers, new File(resDB,"DB.fasta")) : 
					RealtimeSpeciesTypingCmd.getSamIteratorsFQ( fastqFile, readList, maxReads, q_thresh, new File(resDB,"DB.fasta"), null);
				File sample_namek = bamFile!=null ? bamFile[0]: fastqFile[0];
	//	RealtimeSpeciesTypingCmd.getSamIterators(bamFile==null ? null : bamFile, 
	//			fastqFile==null ? null : fastqFile, readList, maxReads, q_thresh, sample_names,iterators, readers, 
	//			new File(resDB,"DB.fasta"));
		
		//for(int k=0; k<iterators.size(); k++){
			File outdir1 =new File(resdir, sample_namek+".resistance");
			
			outdir1.mkdirs();
			RealtimeResistanceGene paTyping = new RealtimeResistanceGene(readPeriod, time, outdir1,outdir1.getAbsolutePath()+"/"+output, resDB.getAbsolutePath(), tmp);		
			paTyping.msa = msa;		
			paTyping.scoreThreshold = scoreThreshold;
			paTyping.twoDOnly = twodonly;
			paTyping.numThead = thread;
			paTyping.typing(samIter);	
			paTyping.close();
			
			paTyping.getOutfiles(outfiles);
		//}
		for(int i=0; i<readers.size(); i++){
			readers.get(i).close();
		}
	
	}
}

/*RST*
-------------------------------------------------------------------------------------------------------
*jsa.np.rtResistGenes*: Antibiotic resistance gene identification in real-time with Nanopore sequencing 
-------------------------------------------------------------------------------------------------------

*jsa.np.rtResistGenes* identifies antibiotic resistance genes from real-time sequencing
with Nanopore MinION. 

<usage> 

~~~~~~~~~~
Setting up
~~~~~~~~~~

Refer to the documentation at https://github.com/mdcao/npAnalysis/ for more 
details.


*RST*/
