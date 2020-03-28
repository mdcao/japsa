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

/*                                                    
 * 01/03/2020 - Lachlan Coin                                       
 ****************************************************************************/

package japsa.tools.bio.hts;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.TranscriptUtils1;
import japsa.util.TranscriptUtils1.Annotation;
import japsa.util.TranscriptUtils1.IdentityProfile1;
import japsa.util.deploy.Deployable;

/**
 * @author Lachlan Coin
 *
 */

@Deployable(scriptName = "jsa.hts.viralanalysis", scriptDesc = "Error analysis of sequencing data")
public class ViralTranscriptAnalysisCmd2 extends CommandLine {
//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

	public ViralTranscriptAnalysisCmd2() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("bamFile", null, "Name of bam file", true);
		addString("breaks", null, "Position File, for looking for specific breaks");
		addString("reference", null, "Name of reference genome", true);
		addString("annotation", null, "ORF annotation file", true);
		addString("resdir", "results"+System.currentTimeMillis(), "results directory");

		addInt("maxReads", Integer.MAX_VALUE, "ORF annotation file");

		addString("pattern", null, "Pattern of read name, used for filtering");
		addInt("qual", 0, "Minimum quality required");
		addInt("bin", 1, "Bin size for coverage");
		addInt("startThresh", 80, "Threshold for having 5'");
		addInt("endThresh", 80, "Threshold for having 3'");
		addDouble("overlapThresh", 0.95, "Threshold for overlapping clusters");
		addBoolean("coexpression", false, "whether to calc coexperssion matrices (large memory for small bin size)");
		addStdHelp();
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLine cmdLine = new ViralTranscriptAnalysisCmd2();
		args = cmdLine.stdParseLine(args);

		String reference = cmdLine.getStringVal("reference");
		int qual = cmdLine.getIntVal("qual");
		int bin = cmdLine.getIntVal("bin");
		String pattern = cmdLine.getStringVal("pattern");
		String bamFile = cmdLine.getStringVal("bamFile");
		String annotFile = cmdLine.getStringVal("annotation");
		String positionsFile = cmdLine.getStringVal("breaks");
		String resdir = cmdLine.getStringVal("resdir");
		boolean coexp = cmdLine.getBooleanVal("coexpression");
		if(coexp && bin <10) {
			throw new Error(" this not good idea");
		}
		double overlapThresh = cmdLine.getDoubleVal("overlapThresh");
		int startThresh = cmdLine.getIntVal("startThresh");
		int endThresh = cmdLine.getIntVal("endThresh");
		int maxReads = cmdLine.getIntVal("maxReads");
		errorAnalysis(bamFile, reference, annotFile, positionsFile, resdir,pattern, qual, bin, coexp, overlapThresh, startThresh, endThresh,maxReads);

		// paramEst(bamFile, reference, qual);
	}

	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 */
	static void errorAnalysis(String bamFiles, String refFile, String annot_file, String positionsFile,  String resdir, String pattern, int qual, int round, final boolean coexp, double overlapThresh, int startThresh, int endThresh, int max_reads) throws IOException {
		boolean cluster_reads = true;
		boolean calcTree = false;
		String[] bamFiles_ = bamFiles.split(":");
		
		
		
		TranscriptUtils1.CigarHash.round = (double) round;
		int len = bamFiles_.length;
		// len = 1;
		
		ArrayList<Sequence> genomes = SequenceReader.readAll(refFile, Alphabet.DNA());

		// get the first chrom
		File resDir =  new File(resdir);
		if(!resDir.exists()) resDir.mkdir();
		else if(!resDir.isDirectory()) throw new RuntimeException(resDir+"should be a directory");
		ArrayList<IdentityProfile1> profiles = new ArrayList<IdentityProfile1>();
	//	List<List<String>> genes_all = new ArrayList<List<String>>();
		PrintWriter genes_all_pw = new PrintWriter(new FileWriter(new File(resdir+"/genes.txt")));
		for (int jj = 0; jj < genomes.size(); jj++) {
			Sequence ref = genomes.get(jj);
			List<Integer[]> pos = new ArrayList<Integer[]>();
			if(positionsFile!=null){
			BufferedReader br = new BufferedReader(new FileReader(new File(positionsFile)));
			List<String> header =Arrays.asList( br.readLine().split(",")); //assuming chrom_index, start, end, gene_name
			String st = "";
			int[] inds = new int[] {header.indexOf("genes"), header.indexOf("start"), header.indexOf("end")};
			List<String> genes = new ArrayList<String>();
			while((st=br.readLine())!=null){
				String[] str = st.split(",");
				
				String gene = str[inds[0]];
				int gene_index = genes.indexOf(gene);
				if(gene_index<0){
					gene_index = genes.size();
					genes_all_pw.println(gene+","+jj+","+gene_index);
					genes.add(gene);
				}
				if(Integer.parseInt(str[0])==jj){
					pos.add(new Integer[] {Integer.parseInt(str[inds[1]])-1, Integer.parseInt(str[inds[2]])-1, gene_index});
				}
			}
			br.close();
			}
			//genes_all.add(genes);
			profiles.add(new IdentityProfile1(ref, resDir, pos, len,jj,coexp, overlapThresh, startThresh, endThresh));

		}
		genes_all_pw.close();
		for (int ii = 0; ii < len; ii++) {
			int source_index = ii;
			int currentIndex = 0;
			Sequence chr = genomes.get(currentIndex);
			
			String bamFile = bamFiles_[ii];
			File bam = new File( bamFile);
			for (int jj = 0; jj < genomes.size(); jj++) {
				profiles.get(jj).updateSourceIndex(ii);
			}
			
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			SamReader samReader = null;// SamReaderFactory.makeDefault().open(new File(bamFile));

			if ("-".equals(bamFile))
				samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				samReader = SamReaderFactory.makeDefault().open(bam);

			SAMRecordIterator samIter = samReader.iterator();
			// Read the reference genome
			

			long totReadBase = 0, totRefBase = 0;
			int numReads = 0;

			int numNotAligned = 0;
		//	int max_reads = Integer.MAX_VALUE;
			// String log =
			// "###Read_name\tRead_length\tReference_length\tInsertions\tDeletions\tMismatches\n";
			for (int cntr = 0; samIter.hasNext() && cntr < max_reads; cntr++) {
				SAMRecord sam = samIter.next();

				if (pattern != null && (!sam.getReadName().contains(pattern)))
					continue;

				// make the read seq
				Sequence readSeq = new Sequence(Alphabet.DNA(), sam.getReadString(), sam.getReadName());
				if (readSeq.length() <= 1) {
					// LOG.warn(sam.getReadName() +" ignored");
					// TODO: This might be secondary alignment, need to do something about it
					continue;
				}

				numReads++;

				if (sam.getReadUnmappedFlag()) {
					numNotAligned++;
					continue;
				}

				int flag = sam.getFlags();

				if (sam.getMappingQuality() < qual) {
					numNotAligned++;
					continue;
				}

				// int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index
				int refIndex = sam.getReferenceIndex();

				// if move to another chrom, get that chrom
				if (refIndex != currentIndex) {
					currentIndex = refIndex;
					chr = genomes.get(currentIndex);
				}

				TranscriptUtils1.identity1(chr, readSeq, sam, profiles.get(currentIndex), source_index, cluster_reads);

				// log+=sam.getReadName() + "\t" + profile.readBase + "\t" + profile.refBase +
				// "\t" + profile.baseIns + "\t" + profile.baseDel + "\t" + profile.mismatch+
				// "\n";

				// numReadsConsidered ++;

			}
			samReader.close();
			
		}
		//genes_all_pw.close();
		
		Annotation annot = new Annotation(new File(annot_file));
		for(int i=0; i<profiles.size(); i++) {
			profiles.get(i).finalise();
			profiles.get(i).getConsensus(annot);
			if(calcTree) profiles.get(i).printTree();
		}
		
		// Done

		// System.out.println(log);
	}

}

/*
 * RST* ----------------------------------------------------------
 * jsa.hts.errorAnalysis*: Error analysis of sequencing data
 * ----------------------------------------------------------
 * 
 * jsa.hts.errorAnalysis* assesses the error profile of sequencing data by
 * getting the numbers of errors (mismatches, indels etc) from a bam file.
 * Obviously, it does not distinguish sequencing errors from mutations, and
 * hence consider mutations as errors. It is best to use with the bam file from
 * aligning sequencing reads to a reliable assembly of the sample.
 * 
 * <usage>
 * 
 * RST
 */