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
 * 10 Dec 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsadev.tools;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.tools.seq.ExtractGeneSequenceCmd;
import japsa.util.DoubleArray;
import japsa.util.Logging;
import japsa.util.ProcessManagement;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import com.google.common.collect.LinkedHashMultimap;

/**
 * @author minhduc
 * A class implementing functions in the PanCoreGenesCmd
 */
public class PanCoreGenes {

	/**********************************************************************/	
	static final String annotationPath   = "Annotations/";
	static final String geneSequencePath = "GeneSequences/";	
	static final String variationPath  = "Variation/" ;
	static final String variationLogPath  = "Variation/logs/" ;
	static final String variationVarPath  = "Variation/var/" ;
	static final String assemblyPath   = "Assemblies/";
	static final String panGenome  = assemblyPath + "blast/pangenome";	
	static final String panGenesFile   = geneSequencePath + "bwaIndex/pangenes.fasta";


	static ArrayList<SampleRecord> sampleList = new ArrayList<SampleRecord>();
	static HashMap<String, SampleRecord> sampleMap = new HashMap<String, SampleRecord>();
	static final int FLANK = 300;
	static final int SMALL_FLANK = 60;
	static int threads = 8;
	/**********************************************************************/

	static void readConfig(String cfgFile) throws IOException{
		BufferedReader bf = SequenceReader.openFile(cfgFile);
		SampleRecord samRecord = null;

		String line;			
		int lineNo = 0;

		while ((line = bf.readLine())!=null){
			lineNo ++;
			line = line.trim();
			if (line.startsWith("#"))
				continue;

			if (line.length() <=1)
				continue;

			if (line.indexOf('=') < 0){
				Logging.error("I expect to see char = at line " + lineNo);
				bf.close();						
				System.exit(1);
			}

			if (line.startsWith("[")){//end of a sample
				if (samRecord != null){
					if (samRecord.complete()){
						Logging.info(samRecord.toString());
						sampleList.add(samRecord);
						samRecord = null;
					}else{
						Logging.error("Sample " + samRecord.sampleID + " not complete at line " + lineNo);
						bf.close();						
						System.exit(1);
					}
				}				
			}//if line

			if (line.startsWith("[SAMPLE")){
				String ID = line.substring(7).trim();
				if (!ID.startsWith("ID")){
					Logging.error("I dont understand line " + lineNo);
					bf.close();						
					System.exit(1);					
				}				
				ID = ID.split("=", 2)[1].trim().replaceAll(" ", "_");

				int endPos = ID.indexOf(']');
				if (endPos > 0)
					ID = ID.substring(0, endPos);

				samRecord = new SampleRecord(ID);
				continue;
			}

			if (line.startsWith("GENUS")){
				if (samRecord == null){
					Logging.error("Wrong specs at line " + lineNo);
					bf.close();						
					System.exit(1);					
				}

				samRecord.genus = line.split("=", 2)[1].trim().replaceAll(" ", "_");;
				continue;
			}
			if (line.startsWith("SPECIES")){
				if (samRecord == null){
					Logging.error("Wrong specs at line " + lineNo);
					bf.close();						
					System.exit(1);					
				}
				samRecord.species = line.split("=", 2)[1].trim().replaceAll(" ", "_");;
				continue;
			}
			if (line.startsWith("STRAIN")){
				if (samRecord == null){
					Logging.error("Wrong specs at line " + lineNo);
					bf.close();						
					System.exit(1);					
				}
				samRecord.strain = line.split("=", 2)[1].trim().replaceAll(" ", "_");
				continue;
			}

			if (line.startsWith("READ_1")){
				if (samRecord == null){
					Logging.error("Wrong specs at line " + lineNo);
					bf.close();						
					System.exit(1);					
				}
				samRecord.readFile1 = line.split("=", 2)[1].trim();
				continue;
			}

			if (line.startsWith("READ_2")){
				if (samRecord == null){
					Logging.error("Wrong specs at line " + lineNo);
					bf.close();						
					System.exit(1);					
				}
				samRecord.readFile2 = line.split("=", 2)[1].trim();
				continue;
			}			

			if (line.startsWith("GENOME")){
				if (samRecord == null){
					Logging.error("Wrong specs at line " + lineNo);
					bf.close();						
					System.exit(1);					
				}
				samRecord.genome = line.split("=", 2)[1].trim();
				continue;
			}

			if (line.startsWith("ANNOTATION")){
				if (samRecord == null){
					Logging.error("Wrong specs at line " + lineNo);
					bf.close();						
					System.exit(1);					
				}
				samRecord.annotation = line.split("=", 2)[1].trim();
				continue;
			}

			if (line.startsWith("TYPE")){
				if (samRecord == null){
					Logging.error("Wrong specs at line " + lineNo);
					bf.close();						
					System.exit(1);					
				}
				samRecord.isReference = line.split("=", 2)[1].trim().toUpperCase().startsWith("REF");
				continue;
			}
		}//while

		if (samRecord != null){
			if (samRecord.complete()){
				sampleList.add(samRecord);		
				Logging.info(samRecord.toString());
				samRecord = null;
			}else{
				Logging.error("Sample " + samRecord.sampleID + " not complete");
				bf.close();						
				System.exit(1);
			}
		}
		bf.close();

		String[] paths = {assemblyPath,annotationPath,geneSequencePath,variationPath,variationLogPath,variationVarPath,assemblyPath + "blast/"};
		for (String path:paths){
			//System.out.println(path);
			File dir = new File(path);			
			if (!dir.isDirectory()){
				Logging.info("Creating path " + dir.getPath());
				if (!dir.mkdirs()){
					Logging.exit("Cannot create path " +dir.getPath(), 1);
				}
			}
		}
	}

	/**
	 * The most important methods, read in the config
	 * @param cfgFile
	 * @throws IOException
	 */
	static void readSampleList(String cfgFile) throws IOException{				
		BufferedReader bf = SequenceReader.openFile(cfgFile);
		String line;		
		int lineNo = 0;
		while ((line = bf.readLine())!=null){
			lineNo ++;
			if (line.startsWith("#"))
				continue;

			String [] toks = line.split("\t");

			if (toks.length < 4){
				Logging.exit("Line " + lineNo + " unexpected", 1);
			}

			SampleRecord sample = new SampleRecord(toks[0],toks[1],toks[2],toks[3]);
			sample.setReadFile1(toks[4]);
			sample.setReadFile2(toks[5]);

			sampleList.add(sample);
		}		
		bf.close();
		Logging.info("Read " + sampleList.size() + " samples");


		//Make sure the folders are there or make them
		String[] paths = {assemblyPath,annotationPath,geneSequencePath,variationPath,variationLogPath,variationVarPath,assemblyPath + "blast/"};
		for (String path:paths){
			//System.out.println(path);
			File dir = new File(path);			
			if (!dir.isDirectory()){
				Logging.info("Creating path " + dir.getPath());
				if (!dir.mkdirs()){
					Logging.exit("Cannot create path " +dir.getPath(), 1);
				}
			}
		}
	}


	/**
	 * Check status of the analysis
	 */
	static void checkStatus(){
		System.out.println("=== Check data : ");		

		System.out.println("=== Check assemblies ======================");	

		System.out.println("Possible commands :");
	}


	/**
	 * Assemble the draft genomes of all samples
	 * 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	static void assemble() throws IOException, InterruptedException{

		boolean changed = false;
		for  (SampleRecord sampleRecord:PanCoreGenes.sampleList){
			String sampleID  = sampleRecord.sampleID;
			String sampleAssembly = assemblyPath + sampleID + ".fasta";
			File file = new File(sampleAssembly);

			if (!file.exists()){
				changed = true;
				if (sampleRecord.genome == null){					
					ProcessBuilder pb = new ProcessBuilder("spades.py", 
						"-k", "21,33,55,77,99,127", 
						"--careful",
						"--pe1-1", String.valueOf(sampleRecord.readFile1),
						"--pe1-2", String.valueOf(sampleRecord.readFile2),
						"-m", "60",
						"-t", String.valueOf(threads),
						"-o", assemblyPath + "spades_output_" + sampleID
						);
					//System.out.println("qsub -b y -cwd -pe smp 12 -N a" +sampleID + " /sw/spades/SPAdes-3.6.0-Linux/bin/spades.py -k 21,33,55,77,99,127 --careful --pe1-1 " +  String.valueOf(sampleRecord.readFile1) + " --pe1-2 " + String.valueOf(sampleRecord.readFile2) + " -m 60 -t 16 -o " + assemblyPath + "spades_output_" + sampleID);
					int status = ProcessManagement.runProcess(pb);
					if (status ==0 ){			
						Logging.info("Successfully assembled " + sampleID);
						//continue;
					}else{
						Logging.exit("Assemble " + sampleID + " FAIL",1);
					}					


					//b. Sort, rename contigs IDs
					String contigFile =	assemblyPath + "spades_output_" + sampleID + "/contigs.fasta";							
					ArrayList<Sequence> seqList = SequenceReader.readAll(contigFile, Alphabet.DNA());

					Collections.sort(seqList, 
						new Comparator<Sequence>(){
						public int compare(Sequence seq1, Sequence seq2) {
							//sort to reverse order
							return seq2.length() - seq1.length();	      
						}
					});

					//c. Rename contig and write to file
					int ind = 0;
					//how many chars needed 
					int fieldSize = String.valueOf(seqList.size()).length();
					SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(sampleAssembly);

					for (Sequence seq:seqList){
						ind ++;
						String name = String.valueOf(ind);
						//pad in 0
						while (name.length() < fieldSize)
							name = "0" + name;

						seq.setDesc(seq.getName() +" " + seq.getDesc());
						seq.setName(sampleID + "|C" + name);
						seq.writeFasta(sos);										
					}
					Logging.info("Assembled sample " + sampleID);
					sos.close();
				}else{//if sample.genome

					ArrayList<Sequence> seqList = SequenceReader.readAll(sampleRecord.genome, Alphabet.DNA());

					//Dont sort it unless the annotation can be sorted as well
					//Collections.sort(seqList, 
					//	new Comparator<Sequence>(){
					//	public int compare(Sequence seq1, Sequence seq2) {
					//		//sort to reverse order
					//		return seq2.length() - seq1.length();	      
					//	}
					//});

					int ind = 0;
					//how many chars needed 
					int fieldSize = String.valueOf(seqList.size()).length();

					HashMap<String, String> mapNames = new HashMap<String, String>();

					SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(sampleAssembly);
					for (Sequence seq:seqList){
						ind ++;
						String name = String.valueOf(ind);
						//pad in 0
						while (name.length() < fieldSize)
							name = "0" + name;

						mapNames.put(seq.getName(), sampleID + "|C" + name); 
						seq.setDesc(seq.getName() +" " + seq.getDesc());
						seq.setName(sampleID + "|C" + name);
						seq.writeFasta(sos);
						//for (String st)
					}
					sos.close();
					Logging.info("Copy assembly for sample " + sampleID);

					if (sampleRecord.annotation != null){
						SequenceOutputStream gffOS = SequenceOutputStream.makeOutputStream(assemblyPath + sampleID + ".gff");
						BufferedReader bf = SequenceReader.openFile(sampleRecord.annotation);
						String line;						
						while ((line = bf.readLine())!=null){
							for (String key:mapNames.keySet()){
								line = line.replaceAll(key, mapNames.get(key));
							}
							gffOS.print(line + "\n");							
						}
						bf.close();				
						gffOS.close();
						Logging.info("Copy annotation for sample " + sampleID);
					}
				}
			}
		}

		//d. build blast index only if there is change

		if (changed){
			SequenceOutputStream panOs = SequenceOutputStream.makeOutputStream(panGenome + ".fasta");
			for  (SampleRecord sampleRecord:PanCoreGenes.sampleList){
				String sampleID  = sampleRecord.sampleID;
				String sampleAssembly = assemblyPath + sampleID + ".fasta";
				ArrayList<Sequence> seqList = SequenceReader.readAll(sampleAssembly, Alphabet.DNA());
				for (Sequence seq:seqList){
					seq.writeFasta(panOs);										
				}					
				Logging.info("Assembly for " + sampleID + " added to blast");				
			}

			panOs.close();

			Logging.info("Building blast index");			
			ProcessBuilder pb = new ProcessBuilder("makeblastdb", 
				"-dbtype", "nucl", 
				"-in",  panGenome + ".fasta", 
				"-title","Pan Genome",
				"-out", panGenome);

			int status = ProcessManagement.runProcess(pb);
			if (status ==0 )
				Logging.info("Building blast index done");
			else{
				Logging.warn("Building blast index FAIL");
			}
		}	
	}

	static void annotate() throws IOException, InterruptedException{
		//2. annotate
		//2aa TODO:submit assemblies to rast to annotate
		//2ab TODO:download annotations back
		//2ac TODO:fix gff format (using jsa.dev.fixRastGff)

		//"prokka
		for  (SampleRecord sampleRecord:PanCoreGenes.sampleList){
			String sampleID = sampleRecord.sampleID;
			String sampleAnnotation = annotationPath  +  sampleID + "/" + sampleID + ".gff";
			File annoFile = new File(sampleAnnotation);
			if (annoFile.exists()){
				System.out.println(sampleID  + " exist ");
				continue;
			}			
			String cmd =   
				"prokka --cpus " + PanCoreGenes.threads 
				+ " --mincontiglen 200 --addgenes --force --debug --outdir " + annotationPath  +  sampleID 
				+ " --genus "   + sampleRecord.genus 
				+ " --species " + sampleRecord.species
				+ " --strain "  + sampleRecord.strain
				+ " --prefix "  + sampleID
				+ " --locustag "+ sampleID;

			//if rast annotation exists
			File file;

			file = new File(assemblyPath + sampleID + ".gff");
			if (file.exists())
				cmd = cmd + " --rast " + file.getPath();
			else{
				file = new File(annotationPath + "RAST/" + sampleID + "/" + sampleID + ".gff");
				if (file.exists())
					cmd = cmd + " --rast " + file.getPath();
			}

			cmd = cmd + " " + assemblyPath + sampleID + ".fasta";
			Logging.info("Run " + cmd);

			ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
			int status = ProcessManagement.runProcess(pb);
			if (status ==0 )
				Logging.info("Annotation of " + sampleID + " done");
			else{
				Logging.warn("Annotation of " + sampleID + " FAIL");
			}
		}
	}


	/**
	 * Extract gene sequences from the annotations, group them
	 * @throws IOException
	 * @throws InterruptedException
	 */
	static void groupGenes() throws IOException, InterruptedException{
		//3. Extract gene sequences
		String cdHitOut = PanCoreGenes.geneSequencePath + "/cd-hit.out";
		String allGeneFile = geneSequencePath + "allgenes.fasta";

		boolean redo = true;
		if (redo){
			SequenceOutputStream out = SequenceOutputStream.makeOutputStream(allGeneFile);
			for  (SampleRecord sampleRecord:sampleList){			
				String sampleID   = sampleRecord.sampleID;			
				//String sequence = sampleRecord.genome;
				//if (sequence == null) 
				String sequence = annotationPath  + sampleID + "/" + sampleID + ".fna";
				String gff = annotationPath  + sampleID + "/" + sampleID + ".gff";

				//String type     = "CDS,tRNA,rRNA,tmRNA";
				String type     = "CDS";
				ExtractGeneSequenceCmd.extractGenes(sequence, gff, type, SMALL_FLANK , out);				
			}
			out.close();

			//cmds.add("/sw/cd-hit/v2015/bin/cd-hit-est");
			String cmd = "cd-hit-est -M 30000 -d 0 -s 0.65 -aS 0.8 -uS 0.3 -uL 0.3 -p 1 -g 1 -T " + PanCoreGenes.threads + " -i " + PanCoreGenes.geneSequencePath + "allgenes.fasta -o " + cdHitOut;			
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s+"));
			Logging.info("Running " + cmd);

			int status = ProcessManagement.runProcess(pb);			
			if (status == 0 )
				Logging.info("Successfully group with cd-hit");
			else{
				Logging.exit("Unsuccessfully group with cd-hit", 1);
			}
		}
		///////////////////////////////////////////////////////////////
		/****************************************************************/
		{
			HashMap<String, GeneRecord> panGenes = new HashMap<String, GeneRecord>(); 
			BufferedReader bf = new BufferedReader(new FileReader(cdHitOut+".clstr"));
			String line = "";
			//ArrayList<GeneRecord> clusters = new ArrayList<GeneRecord>();
			GeneRecord mainGene = null;			
			int clusterIndex = 0;
			while ((line = bf.readLine())!=null){			
				if (line.startsWith(">")){					
					if (mainGene != null){						
						panGenes.put(mainGene.geneID, mainGene);
						System.out.println(">Cluster " + clusterIndex);
						mainGene = null;
					}
					clusterIndex ++;
					continue;
				}			
				String [] toks = line.trim().split(" ");


				if (toks.length < 3){
					bf.close();
					Logging.exit("Not expected 1", -1);
				}
				String name = toks[1].substring(1, toks[1].length()-3);
				GeneRecord gene = new GeneRecord(name);

				if (mainGene == null ||gene.compareTo(mainGene) < 0){
					mainGene = gene;
				}				
			}//while
			bf.close();
			if (mainGene != null){						
				panGenes.put(mainGene.geneID, mainGene);
				System.out.println(">Cluster " + clusterIndex);
			}

			//Get description of the gene
			Logging.info("Getting gene description");
			bf = new BufferedReader(new FileReader(allGeneFile));
			while ((line = bf.readLine())!=null){		
				if (!line.startsWith(">"))
					continue;

				String [] toks = line.trim().substring(1).split("\\s+",2);
				GeneRecord gene = panGenes.get(toks[0]);
				if (gene != null)
					gene.geneDesc = toks[1];
			}			


			//get gene sequences			
			ArrayList<Sequence> panGenones = SequenceReader.readAll(panGenome + ".fasta", Alphabet.DNA16());
			HashMap<String, Sequence> panGenomeMap = new HashMap<String, Sequence> ();

			for (Sequence seq:panGenones)
				panGenomeMap.put(seq.getName(), seq);

			Logging.info("Getting gene sequences for " + panGenes.size() + " from " + panGenomeMap.size());
			int index = 0;
			int fieldLength = 6;

			SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(panGenesFile);
			SequenceOutputStream gffOS = SequenceOutputStream.makeOutputStream(geneSequencePath + "pangenes.gff");
			gffOS.print("##gff-version 3\n");
			int geneIndex = 0;

			for (GeneRecord gene:panGenes.values()){
				String contig = gene.contig;
				Sequence seq = panGenomeMap.get(contig);

				Sequence geneWithFlank = seq.subsequenceWithFlank(gene.start - 1, gene.end, FLANK);
				if (gene.strand == '-')
					geneWithFlank = Alphabet.DNA.complement(geneWithFlank);
				//Got the gene

				index ++;
				String name = String.valueOf(index);
				while (name.length() < fieldLength)
					name = "0" + name;

				name = "G" + name;

				geneWithFlank.setDesc("Pos=" + gene.geneID + ";" + gene.geneDesc);
				geneWithFlank.setName(name);
				geneWithFlank.writeFasta(sos);

				gffOS.print("##sequence-region " + name + " 1 " + geneWithFlank.length() + '\n');
				gffOS.print(name + "\tpcgene\tgene\t" + (FLANK+1) + "\t" + (geneWithFlank.length() - FLANK)+"\t.\t+\t0\t" + gene.geneDesc + '\n');
				gffOS.print(name + "\tpcgene\tCDS\t" + (FLANK+1) + "\t" + (geneWithFlank.length() - FLANK)+"\t.\t+\t0\t" + gene.geneDesc + '\n');

			}
			gffOS.close();
			sos.close();
		}
		/*****************************************************************/

		ProcessBuilder pb = new ProcessBuilder("bwa", "index", panGenesFile);
		int status = ProcessManagement.runProcess(pb);			
		if (status == 0 )
			Logging.info("bwa index successful");
		else{
			Logging.exit("bwa index fail", 1);
		}

		pb = new ProcessBuilder("samtools", "faidx", panGenesFile);
		status = ProcessManagement.runProcess(pb);			
		if (status == 0 )
			Logging.info("samtools index successful");
		else{
			Logging.exit("samtools index fail", 1);
		}		
	}

	static void alignment() throws IOException, InterruptedException{
		int jobNo = 0;
		for  (SampleRecord sampleRecord:PanCoreGenes.sampleList){
			jobNo ++;


			System.out.println("#=#" + jobNo);
			String sampleID = sampleRecord.sampleID;
			String file1 = sampleRecord.readFile1;
			String file2 = sampleRecord.readFile2;
			if (file1 == null)
				file1 = assemblyPath + sampleID + "_R1.fq.gz";

			if (file2 == null)
				file2 = assemblyPath + sampleID + "_R2.fq.gz";

			String cmd = "bwa mem -t " + threads + " -R \"@RG\\tID:"+sampleID+"\\tSM:" + sampleID+"\\tPL:ILLUMINA\" "  + panGenesFile + " " + file1+ " " + file2 + " > " + variationPath + sampleID + ".sam";
			System.out.println(cmd);
			cmd = "samtools view -bSu " + variationPath + sampleID + ".sam | samtools sort - " + variationPath + sampleID;
			System.out.println(cmd);
			System.out.println();
		}		
	}


	/**
	 * TODO: check this function
	 * @param species
	 * @throws IOException
	 */
	static void getSNPSites() throws IOException{

		//I want to preserve the order of keys as well as value hence LinkedHashMultimap
		LinkedHashMultimap<String, SampleRecord> speciesMap = LinkedHashMultimap.create();		
		for (SampleRecord sampleRecord:sampleList){
			String speciesName = (sampleRecord.genus + " " + sampleRecord.species).replaceAll(" ", "_");
			speciesMap.put(speciesName, sampleRecord);			
		}

		for (String key:speciesMap.keySet()){			
			Set<SampleRecord> strainSet = speciesMap.get(key);
			System.out.println(key + " " + strainSet.size());
			for (SampleRecord sampleRecord:strainSet){
				System.out.println("    " + sampleRecord.sampleID);				
			}

			//if (strainSet.size() < 3) continue;


			HashMap<String, SequenceBuilder> seqs = new HashMap<String, SequenceBuilder>();			
			for  (SampleRecord sampleRecord:strainSet){
				SequenceBuilder seq = new SequenceBuilder(Alphabet.DNA(),8012);
				seq.setName(sampleRecord.sampleID);

				seqs.put(sampleRecord.sampleID, seq);
			}

			SequenceReader panGenesReader = SequenceReader.getReader(panGenesFile);
			SequenceOutputStream pos = SequenceOutputStream.makeOutputStream(key + "_snps.pos");

			Sequence genSeq;
			while ((genSeq = panGenesReader.nextSequence(Alphabet.DNA()))!=null){			
				collectSNPFromVCF(genSeq,strainSet,seqs, pos,true);
			}
			panGenesReader.close();
			pos.close();


			SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(key + "_snps.fasta");		
			for (SequenceBuilder seq:seqs.values()){			
				System.out.println(seq.getName() + " " + seq.length());
				seq.writeFasta(sos);
			}
			sos.close();
		}
		//Get set of samples
	}


	private static void collectSNPFromVCF(Sequence gene, Set<SampleRecord> sampleSet, HashMap<String, SequenceBuilder> seqs, SequenceOutputStream pos, boolean onlyCoding) throws IOException{
		VCFFileReader vcf = new  VCFFileReader(new File(variationPath +  "var/" + gene.getName() + ".vcf.gz"));
		CloseableIterator<VariantContext> iter = vcf.iterator();

		char [] SNPs = new char[sampleSet.size()];

		//Go through the vcf
		while (iter.hasNext()){
			VariantContext var = iter.next();

			int snpSite = var.getStart();			
			if (onlyCoding){
				if ((snpSite <= FLANK))
					continue;

				if ((snpSite > gene.length() - FLANK))
					break;
			}							

			GenotypesContext gTypes = var.getGenotypes();
			Arrays.fill(SNPs, 'N');
			int index = 0;

			boolean good = true;

			for (SampleRecord sampleRecord:sampleSet){
				String sampleID = sampleRecord.sampleID;				
				Genotype gType = gTypes.get(sampleID);
				if (gType == null){
					System.err.println("null " + sampleID);
					good = false;
					break;//for
				}	

				if (gType.getAlleles().size()  < 1){
					System.err.println("<1 " + sampleID);
					good = false;
					break;//for
				}	

				Allele allele = gType.getAllele(0);

				String bases = allele.getDisplayString();
				if (bases.length() != 1){
					good = false;
					break;//for
				}

				SNPs[index] = bases.charAt(0);
				index ++;					

			}//for			
			if (!good)
				continue;//while -- move to the next var

			if (index != SNPs.length){
				Logging.exit("Unexpected  " + index + " vs " + SNPs.length,1);				
			}
			pos.print(gene.getName() + "\t" + snpSite + "\n");

			index = 0;
			for (SampleRecord sampleRecord:sampleSet){
				SequenceBuilder seq = seqs.get(sampleRecord.sampleID);
				switch (SNPs[index]){
				case 'A':
				case 'a':
					seq.append((byte)Alphabet.DNA.A);
					break;

				case 'C':
				case 'c':
					seq.append((byte)Alphabet.DNA.C);
					break;

				case 'G':
				case 'g':
					seq.append((byte)Alphabet.DNA.G);
					break;

				case 'T':
				case 't':
					seq.append((byte)Alphabet.DNA.T);
					break;

				default:
					seq.append((byte)Alphabet.DNA.N);
				}
				index ++;
			}			

		}
		vcf.close();
	}


	static void associate(String phenoFile, String vafFile) throws IOException{

		//1. Read phenotypes
		BufferedReader bf = SequenceReader.openFile(phenoFile);

		String line = bf.readLine();
		String [] toks = line.trim().split("\t");

		ArrayList<String> pSample = new ArrayList<String>();
		//ArrayList<double[]> pPhenotypes = new ArrayList<double[]>();		
		HashMap<String, double[]> pMap = new HashMap<String, double[]>(); 

		ArrayList<String> pType = new ArrayList<String>();


		for (int i = 1; i < toks.length;i++){
			pType.add(toks[i]);
		}
		while ((line = bf.readLine())!=null){
			if (line.startsWith("#"))
				continue;
			toks = line.trim().split("\t");
			double [] t = new double[pType.size()];

			for (int i = 0; i < t.length;i++){
				try{
					t[i] = Double.valueOf(toks[i+1]);
				}catch(Exception e){
					t[i] = Double.NaN;
				}
			}//for
			pMap.put(toks[0],t);
			pSample.add(toks[0]);
		}//while		
		bf.close();

		//for (String sample:pSample){
		//	System.out.print(sample);
		//	double [] t = pMap.get(sample);
		//	for (int i = 0; i < t.length;i++){
		//		if (t[i] == Double.NaN)
		//			System.out.print("  -");
		//		else					
		//			System.out.print("  " + t[i]);
		//	}
		//	System.out.println();			
		//}


		PearsonsCorrelation pearsonTest = new PearsonsCorrelation();

		//2. Read genotypes
		bf = SequenceReader.openFile(vafFile);		
		line = bf.readLine();		
		ArrayList<String> gSample = new ArrayList<String>();
		toks = line.trim().split("\t");
		for (int i = 1; i < toks.length;i++){
			gSample.add(toks[i]);
		}

		HashMap<String, Integer>  gMap = new HashMap<String, Integer>();

		DoubleArray pValues = new DoubleArray();
		DoubleArray gValues = new  DoubleArray();

		while ((line = bf.readLine())!=null){
			if (line.startsWith("#"))
				continue;


			toks = line.trim().split("\t");
			String var = toks[0];

			gMap.clear();
			int typeIndex = 0;

			int [] gTypes = new int[gSample.size()];

			for (int x = 0; x< gTypes.length;x++){
				String  myGTypeStr = toks[x+1];
				Integer myGType = gMap.get(myGTypeStr);
				if (myGType == null){
					myGType = typeIndex;
					gMap.put(myGTypeStr, myGType);					
					typeIndex ++;
				}
				gTypes[x] = myGType;				
			}
			if (typeIndex <=1)
				continue;//while


			System.out.print("##" + var);
			for (String key:gMap.keySet()){
				System.out.print(" " + gMap.get(key)+":" + key);
			}
			System.out.println();
			
			
			//System.out.print(var

			for (int i = 0; i < pType.size();i++){
				String type = pType.get(i);
				pValues.clear();
				gValues.clear();

				for (int x = 0; x < gSample.size();x++){
					double[] pThisSample = pMap.get(gSample.get(x));
					if (pThisSample != null){
						if (!Double.isNaN(pThisSample[i])){
							pValues.add(pThisSample[i]);
							gValues.add(gTypes[x]);
						}
					}//if
				}//for x
				double pearson = -100;
				try{
					//double corr = new PearsonsCorrelation().correlation(y, x);
					pearson = new PearsonsCorrelation().correlation(pValues.toArray(), gValues.toArray());
				}catch(Exception e){
					pearson = -100;
				}

				if (pearson > -100){				
					System.out.print(var + " " + typeIndex + " " + type + " " + pearson + " " + pValues.size());
					//System.out.print("[");
					//for (int y=0; y< pValues.size();y++)
					//	System.out.print(pValues.get(y)+ ",");
					//System.out.print("][");
					//for (int y=0; y< gValues.size();y++)
					//	System.out.print(gValues.get(y)+ ",");
					//System.out.print("]");				
					System.out.println();
				}
			}//for i
		}//while
		bf.close();
	}

	//Assume FreeBayes has run
	static void analyseVar() throws IOException{		
		ArrayList<VariantContext> varList = new ArrayList<VariantContext>();
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream("-");

		sos.print("##LOCI");
		for (SampleRecord sampleRecord:sampleList){			
			sos.print("\t" + sampleRecord.sampleID);
		}
		sos.println();

		SequenceReader reader = SequenceReader.getReader(PanCoreGenes.panGenesFile);
		Sequence seq;
		//int varCount = 0;
		while ( (seq = reader.nextSequence(Alphabet.DNA16())) != null){
			String vcfFile = variationPath +  "var/" + seq.getName() + ".vcf.gz";			
			collectVar(vcfFile, varList,sos);
		}		
		reader.close();		
		sos.print("##There are " + varList.size());
		sos.println();
		sos.close();
	}

	static void analyseRefVar(String vcfFile) throws IOException{
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream("-");		
		sos.print("##LOCI");
		for (SampleRecord sampleRecord:sampleList){
			sos.print("\t" + sampleRecord.sampleID);
		}
		sos.println();

		ArrayList<VariantContext> varList = new ArrayList<VariantContext>();
		collectVar(vcfFile, varList,sos);
		sos.close();
	}	

	private static void collectVar(String vcfFile, ArrayList<VariantContext> varList, SequenceOutputStream sos) throws IOException{		
		VCFFileReader vcf = new  VCFFileReader(new File(vcfFile));//variationPath +  "var/" + gene.getName() + ".vcf.gz"));
		CloseableIterator<VariantContext> iter = vcf.iterator();

		HashMap<String, Integer> typeMap = new HashMap<String, Integer>();

		//Go through the vcf
		while (iter.hasNext()){
			VariantContext var = iter.next();
			varList.add(var);

			//		int type = 0;

			GenotypesContext gTypes = var.getGenotypes();
			sos.print(var.getChr()+":" + var.getStart());			

			//for (SampleRecord sampleRecord:sampleSet){
			for (SampleRecord sampleRecord:sampleList){				
				String sampleID = sampleRecord.sampleID;				
				Genotype gType = gTypes.get(sampleID);
				sos.print('\t');


				//				String myGType = "";

				if (gType == null){
					sos.print('.');
					//myGType = ".";					
				}else{
					for (int x = 0; x < gType.getAlleles().size();x++){
						Allele allele = gType.getAllele(x);
						sos.print(allele.getBaseString());
						sos.print('|');
					}//for
				}

				//Integer myType = typeMap.get(myGType);
				//if (myType == null){
				//	myType = new Integer(type);
				//	typeMap.put(myGType, myType);					
				//	type ++;
				//}

				//sos.print("\t" + " " + myGType);
			}//for
			sos.println();			
		}//while
		vcf.close();
	}


	/**
	 * Run freebayes to call variations from bam files against the core genes.
	 * Paralellisation obtained by running on each core-gene separately
	 * 
	 * @param sampleList
	 * @throws IOException
	 * @throws InterruptedException
	 */
	static void runFreeBayes() throws IOException, InterruptedException{
		String freeBayesBin = "/sw/freebayes/current/bin/freebayes";

		ArrayList<String> fbCmd = new ArrayList<String>();		
		fbCmd.add(freeBayesBin);

		fbCmd.add("--fasta-reference");
		fbCmd.add(PanCoreGenes.panGenesFile);

		for (SampleRecord sampleRecord:sampleList){
			String sampleID = sampleRecord.sampleID;
			fbCmd.add("-b");
			fbCmd.add(PanCoreGenes.variationPath + sampleID + ".bam");				
		}

		SequenceReader reader = SequenceReader.getReader(PanCoreGenes.panGenesFile);
		Sequence seq;

		while ( (seq = reader.nextSequence(Alphabet.DNA16())) != null){
			ArrayList<String> qsubCmd = new ArrayList<String>();
			qsubCmd.add("qsub");
			qsubCmd.add("-b");
			qsubCmd.add("y");

			qsubCmd.add("-cwd");

			//qsubCmd.add("-pe");
			//qsubCmd.add("1");

			qsubCmd.add("-N");
			qsubCmd.add("fb"+seq.getName());

			qsubCmd.add("-o");
			qsubCmd.add(variationLogPath + seq.getName() + ".o");

			qsubCmd.add("-e");
			qsubCmd.add(variationLogPath +seq.getName() + ".e");


			qsubCmd.addAll(fbCmd);			

			qsubCmd.add("--vcf");
			qsubCmd.add(variationVarPath + seq.getName() + ".vcf");

			qsubCmd.add("--region");
			qsubCmd.add(seq.getName());	

			ProcessBuilder pb = new ProcessBuilder(qsubCmd);
			pb.redirectErrorStream(true);

			Process process = pb.start();
			BufferedReader pbOut = new BufferedReader(new InputStreamReader(process.getInputStream()));
			String outStr = "", outLine = "";			 

			while ((outLine = pbOut.readLine())!=null){
				outStr += outLine.trim() + "==";
			}			
			pbOut.close();

			//getCmd(pb);

			int status = process.waitFor();

			if (status != 0){
				Logging.warn(getCmd(pb) + " #for " + seq.getName() + " not succcessful!!! : " + outStr);
			}else{
				Logging.info(getCmd(pb) + " #for " + seq.getName() + " succcessful");
			}

			int jobsInQueue = 2000;

			while (jobsInQueue > 1000){
				Thread.sleep(200);
				//check status
				ProcessBuilder qstatPb = new ProcessBuilder("qstat");
				Process qstatProc = qstatPb.start();
				jobsInQueue = 0;
				BufferedReader br = new BufferedReader(new InputStreamReader(qstatProc.getInputStream()));
				String line;

				while ((line = br.readLine()) != null) {
					jobsInQueue ++;
				}

				if (jobsInQueue > 1000){
					Logging.info("Queue after " + seq.getName() + " (" + jobsInQueue +")");
				}				
				//System.out.println(getCmd(pb));
			}
		}
	}

	static void genePresence() throws IOException, InterruptedException{
		String blastn = "blastn";

		ProcessBuilder pb = new ProcessBuilder(blastn, 
			"-db", panGenome,
			"-query", panGenesFile,
			"-num_threads",String.valueOf(threads),			
			"-outfmt", "7 qseqid qlen qstart qend sseqid slen sstart send length frames pident nident gaps mismatch score bitscore");		                

		Process process = pb.start();

		double minCov = 0.6;//Minimum coverage of the gene
		double minID =  0.85;//minim sequence identity

		BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));
		String line;

		while ((line = br.readLine()) != null) {
			if (line.startsWith("# Query: ")){
				System.out.println(line);
				continue;
			}
			if (line.startsWith("#"))
				continue;

			String [] toks = line.trim().split("\t");
			int length = Integer.parseInt(toks[8]);//length of the alignment
			int qlen   = Integer.parseInt(toks[1]);//length of the query
			if (minCov * qlen > length)
				continue;

			//sequence identity = toks[10]
			if (Double.parseDouble(toks[10]) < minID * 100)
				continue;
			//pass		
			String geneID = toks[0];
			String contig = toks[4];			

			System.out.println(geneID + "\t" + contig + "\t" + toks[6] + "\t" + toks[7] + "\t" + toks[2] + "\t" + toks[3] + "\t" + toks[8] + "\t" + toks[9]);

			//res.add(toks[0]);			
		}
		br.close();
		int status = process.waitFor();//Do i need this???
		//	return status;
	}

	private static String getCmd(ProcessBuilder pb){
		String cmd = "";
		List<String> cmds = pb.command();

		for (String s:cmds)
			cmd = cmd + s + " ";		

		return cmd;
	}

	static private void buildString(StringBuilder builder, String str, int field){
		if (str.length() < field){
			builder.append(str);
			for (int x = str.length();x < field;x++)
				builder.append(' ');	
		}else
			builder.append(str.substring(0, field));
	}

	static void runBlast() throws IOException, InterruptedException{

		SequenceOutputStream outputFull = SequenceOutputStream.makeOutputStream("geneVarFull.dat");
		SequenceOutputStream outputPart = SequenceOutputStream.makeOutputStream("geneVar.dat");

		SequenceReader panReader = SequenceReader.getReader(panGenesFile);		
		String blastn = "blastn";		 
		ProcessBuilder pb = new ProcessBuilder(blastn, 
			"-db", panGenome,
			"-query", panGenesFile,
			"-num_threads",String.valueOf(threads),			
			"-outfmt", "5",
			"-out", "-");

		Process process = pb.start();
		BufferedReader pbOut = new BufferedReader(new InputStreamReader(process.getInputStream()));
		String outLine = "";			 

		String queryDefPrefix = "<Iteration_query-def>";
		String queryLengthPrefix = "<Iteration_query-len>";		

		String hitDefPrefix = "<Hit_def>";
		String hitLengthPrefix = "<Hit_len>";

		String hspBitScorePrefix = "<Hsp_bit-score>";
		String hspScorePrefix  = "<Hsp_score>";
		String hspEValuePrefix ="<Hsp_evalue>";
		String hspQueryFromPrefix = "<Hsp_query-from>";
		String hspQueryToPrefix = "<Hsp_query-to>";
		String hspHitFromPrefix = "<Hsp_hit-from>";
		String hspHitToPrefix = "<Hsp_hit-to>";
		String hspQueryFramePrefix = "<Hsp_query-frame>";
		String hspHitFramePrefix = "<Hsp_hit-frame>";
		String hspIdenPrefix = "<Hsp_identity>";
		String hspPositvePrefix = "<Hsp_positive>";
		String hspGapPrefix = "<Hsp_gaps>";
		String hspAlignLenPrefix = "<Hsp_align-len>";	
		String hspQuerySeqPrefix = "<Hsp_qseq>";
		String hspHitSeqPrefix = "<Hsp_hseq>";
		String hspMidlinePrefix = "<Hsp_midline>";
		String interEndPrefix = "</Iteration>";


		String queryName = "";
		String queryInfo = "";		
		String hitName = "";
		String hitInfo = "";		
		int queryLength = 0, hitLength = 0;
		double bitScore = 0, hspScore = 0, eValue;
		int hitFrom, hitTo, queryFrom, queryTo;
		int queryFrame = 0, hitFrame = 0;

		//ArrayList<HSPAlignment> alignments = new ArrayList<HSPAlignment>();
		//ArrayList<HSPAlignment> finalAlignments = new ArrayList<HSPAlignment>();
		ArrayList<String> hspNames = new ArrayList<String>();
		HashSet<String>   hspStrainSet = new HashSet<String>(); 


		int arrayIndex = -1;
		Sequence geneSequence = null;
		NodeAlignment head = null;//new NodeAlignment(aGene.symbolAt(0));
		NodeAlignment tail = null;// head;

		int lineNo = 0;		
		while ((outLine = pbOut.readLine())!=null){
			lineNo ++;
			outLine = outLine.trim();
			if (outLine.startsWith(interEndPrefix)){
				//End of an query iternation, print out results
				if (arrayIndex >= 0){
					int countVar = 0;
					StringBuilder queryFullStr = new StringBuilder();
					StringBuilder queryPartStr = new StringBuilder();
					buildString(queryFullStr, queryName, 29);queryFullStr.append(' ');
					buildString(queryPartStr, queryName, 29);queryPartStr.append(' ');


					ArrayList<StringBuilder> fullStr = new ArrayList<StringBuilder>(arrayIndex + 1);
					ArrayList<StringBuilder> partStr = new ArrayList<StringBuilder>(arrayIndex + 1);
					for (int i = 0; i <= arrayIndex;i++){
						StringBuilder fBuilder = new StringBuilder();						
						StringBuilder pBuilder = new StringBuilder();

						buildString(fBuilder, hspNames.get(i),29);fBuilder.append(' ');
						buildString(pBuilder, hspNames.get(i),29);pBuilder.append(' ');

						fullStr.add(fBuilder);
						partStr.add(pBuilder);												
					}
					NodeAlignment node = head;
					int genePos = -1;//-1 = promoter, 0= cds, 1 = post
					while (node != null){
						//check if end of CDS
						if (genePos == 0 && node.refPos == queryLength - FLANK){
							genePos = 1;

							queryFullStr.append('|');
							for (StringBuilder builder:fullStr)
								builder.append('|');

							queryPartStr.append('|');
							for (StringBuilder builder:partStr)
								builder.append('|');
						}

						boolean same = true;
						queryFullStr.append(node.refNucl);
						for (int i = 0; i<=arrayIndex;i++){
							fullStr.get(i).append(node.site[i]);
							if (node.site[i] != node.refNucl)
								same = false;
						}

						if (!same){
							queryPartStr.append(node.refNucl);
							for (int i = 0; i<=arrayIndex;i++){
								partStr.get(i).append(node.site[i]);
							}	
							countVar ++;
						}
						//check if start of CDS
						if (genePos < 0 && node.refPos == FLANK - 1){
							genePos = 0;

							queryFullStr.append('|');
							for (StringBuilder builder:fullStr)
								builder.append('|');

							queryPartStr.append('|');
							for (StringBuilder builder:partStr)
								builder.append('|');							
						}	
						node = node.next;

					}//while


					///Print full version
					outputFull.print("============================================================================================\n");
					outputFull.print(queryName + " : " + queryInfo);
					outputFull.println();

					outputFull.print(queryFullStr.toString());
					outputFull.println();

					for (StringBuilder builder:fullStr){
						outputFull.print(builder.toString());
						outputFull.println();
					}
					outputFull.print("###Count var = " + countVar + " for " + queryName + " : " + queryInfo);
					outputFull.println();

					//print part version										
					outputPart.print("============================================================================================\n");
					outputPart.print(queryName + " : " + queryInfo);
					outputPart.println();

					outputPart.print(queryPartStr.toString());
					outputPart.println();

					for (StringBuilder builder:partStr){
						outputPart.print(builder.toString());
						outputPart.println();
					}
					outputPart.print("###Count var = " + countVar + " for " + queryName + " : " + queryInfo);					
					outputPart.println();
				}//if
				continue;
			}

			if (outLine.startsWith(queryDefPrefix)){
				//starting new iteration
				outLine = outLine.substring(queryDefPrefix.length(), outLine.length() - queryDefPrefix.length() - 1);
				String [] toks = outLine.split("\\s",2);
				queryName = toks[0];
				queryInfo = (toks.length > 1)?toks[1]:"";

				geneSequence = panReader.nextSequence(Alphabet.DNA());
				head = tail = new NodeAlignment(geneSequence.charAt(0));
				head.refPos = 0;
				for (int i = 1; i< geneSequence.length();i++){
					NodeAlignment node = new NodeAlignment(geneSequence.charAt(i));
					node.refPos = i;
					node.prev = tail;
					tail.next = node;
					tail = node;		
				}

				hspNames.clear();
				hspStrainSet.clear();
				arrayIndex = -1;//Reset

				//assert geneSequence.name == queryName

			}else if (outLine.startsWith(queryLengthPrefix)){
				outLine = outLine.substring(queryLengthPrefix.length(), outLine.length() - queryLengthPrefix.length() - 1);
				queryLength = Integer.parseInt(outLine);
			}else if (outLine.startsWith(hitDefPrefix)){
				outLine = outLine.substring(hitDefPrefix.length(), outLine.length() - hitDefPrefix.length() - 1);
				String [] toks = outLine.split("\\s",2);
				hitName = toks[0];
				hitInfo = (toks.length > 1)?toks[1]:"";				
			}else if (outLine.startsWith(hitLengthPrefix)){
				outLine = outLine.substring(hitLengthPrefix.length(), outLine.length() - hitLengthPrefix.length() - 1);
				hitLength = Integer.parseInt(outLine);
			}
			else if ("<Hsp>".equals(outLine)){
				pbOut.readLine();//<Hsp_num>				

				outLine = pbOut.readLine().trim();//<Hsp_bit-score>2433.16</Hsp_bit-score>
				outLine = outLine.substring(hspBitScorePrefix.length(), outLine.length() - hspBitScorePrefix.length() - 1);
				bitScore = Double.parseDouble(outLine);


				outLine = pbOut.readLine().trim();//<Hsp_score>1317</Hsp_score>
				outLine = outLine.substring(hspScorePrefix.length(), outLine.length() - hspScorePrefix.length() - 1);
				hspScore = Double.parseDouble(outLine);

				outLine = pbOut.readLine().trim();//<Hsp_evalue>0</Hsp_evalue>
				outLine = outLine.substring(hspEValuePrefix.length(), outLine.length() - hspEValuePrefix.length() - 1);
				eValue = Double.parseDouble(outLine);

				outLine = pbOut.readLine().trim();//<Hsp_query-from>1</Hsp_query-from>
				outLine = outLine.substring(hspQueryFromPrefix.length(), outLine.length() - hspQueryFromPrefix.length() - 1);
				queryFrom = Integer.parseInt(outLine);


				outLine = pbOut.readLine().trim();//<Hsp_query-to>1317</Hsp_query-to>
				outLine = outLine.substring(hspQueryToPrefix.length(), outLine.length() - hspQueryToPrefix.length() - 1);
				queryTo = Integer.parseInt(outLine);				

				outLine = pbOut.readLine().trim();//<Hsp_hit-from>46233</Hsp_hit-from>
				outLine = outLine.substring(hspHitFromPrefix.length(), outLine.length() - hspHitFromPrefix.length() - 1);
				hitFrom = Integer.parseInt(outLine);

				outLine = pbOut.readLine().trim();//<Hsp_hit-to>47549</Hsp_hit-to>
				outLine = outLine.substring(hspHitToPrefix.length(), outLine.length() - hspHitToPrefix.length() - 1);
				hitTo = Integer.parseInt(outLine);

				outLine = pbOut.readLine().trim();//<Hsp_query-frame>1</Hsp_query-frame>
				outLine = outLine.substring(hspQueryFramePrefix.length(), outLine.length() - hspQueryFramePrefix.length() - 1);
				queryFrame = Integer.parseInt(outLine);
				//queryFrane == 1

				outLine = pbOut.readLine().trim();//<Hsp_hit-frame>1</Hsp_hit-frame>
				outLine = outLine.substring(hspHitFramePrefix.length(), outLine.length() - hspHitFramePrefix.length() - 1);
				hitFrame = Integer.parseInt(outLine);				

				outLine = pbOut.readLine().trim();//<Hsp_identity>1317</Hsp_identity>

				outLine = pbOut.readLine().trim();//<Hsp_positive>1317</Hsp_positive>

				outLine = pbOut.readLine().trim();//<Hsp_gaps>0</Hsp_gaps>

				outLine = pbOut.readLine().trim();//<Hsp_align-len>1317</Hsp_align-len>

				outLine = pbOut.readLine().trim();//<Hsp_qseq>
				String qSeq = outLine.substring(hspQuerySeqPrefix.length(), outLine.length() - hspQuerySeqPrefix.length() - 1);

				outLine = pbOut.readLine().trim();//<Hsp_hseq>
				String hSeq = outLine.substring(hspHitSeqPrefix.length(), outLine.length() - hspHitSeqPrefix.length() - 1);

				outLine = pbOut.readLine().trim();//<Hsp_midline>
				String midLineStr = outLine.substring(hspMidlinePrefix.length(), outLine.length() - hspMidlinePrefix.length() - 1);

				//Checking
				int left = FLANK;
				if (queryFrom > left)
					left = queryFrom;

				int right = queryLength - FLANK - 1;
				if (queryTo < right)
					right = queryTo;

				if ((right - left) < ((queryLength -  FLANK - FLANK) * 0.5)){
					//cover < half of CDS
					continue;
				}				

				//System.out.println(queryName + " (" + queryLength +") " + queryInfo + " vs " + hitName + ":" + queryFrame + " vs " + hitFrame);
				if (arrayIndex >= NodeAlignment.LENGTH - 1){
					//TODO: extend this
					continue;					
				}


				String strain = hitName.split("\\|")[0];

				if (hspStrainSet.contains(strain)){ 
					continue;
				}

				hspStrainSet.add(strain);				

				arrayIndex ++;
				//while (arrayIndex >= alignments.size())
				//	alignments.add(new HSPAlignment());


				hspNames.add(hitName + ":" + hitFrom + "-" + hitTo);	


				//HSPAlignment alignment = alignments.get(arrayIndex);				
				//alignment.queryID = queryName;
				//alignment.hitID   = hitName;
				//alignment.queryFrom = queryFrom;
				//alignment.queryTo   = queryTo;
				//alignment.hitFrom = hitFrom;
				//alignment.hitTo   = hitTo;
				//alignment.qSeq = qSeq;
				//alignment.hSeq = hSeq;
				//alignment.sameFrame = (hitFrame == 1);
				//if (hitFrame == queryFrame){
				//	alignment.hitFrom = hitFrom;
				//	alignment.hitTo   = hitTo;
				//}else{
				//	alignment.hitTo = hitFrom;
				//	alignment.hitFrom  = hitTo;
				//}
				NodeAlignment node = head;

				//search for the first alignment
				int pos = node.refPos;//should be 0
				while (pos < queryFrom - 1){
					node.site[arrayIndex] = '*';//DNA.GAP;

					node = node.next;
					if (node == null){
						throw new RuntimeException("This node shouldnt be null");
					}
					pos = node.refPos;					
				}

				//node point to the next point to fill
				for (int i = 0; i< qSeq.length();i++){
					char nuclQuery = qSeq.charAt(i);
					char nuclHit   = hSeq.charAt(i);

					if (nuclQuery == '-'){//1.insertion
						if (node.prev == null){
							throw new RuntimeException("This node shouldnt be head");
						}else if (node.refNucl == '-'){ 
							//1.1 This site is already an insertion						
							node.site[arrayIndex] = nuclHit;
							node = node.next;
						}else{//not yet an insertion, need to make a node before current node
							NodeAlignment newNode = new NodeAlignment('-');
							newNode.refPos = node.refPos - 1;
							for (int j = 0; j < arrayIndex; j++){
								if (node.site[j] == '*' || node.prev.site[j] =='*')
									newNode.site[j] = '*';
								else
									newNode.site[j] = '-';

							}

							newNode.site[arrayIndex] = nuclHit;

							newNode.prev = node.prev;
							newNode.next = node;
							node.prev.next = newNode;
							node.prev = newNode;							
						}												

					}else{//3.alignment and deletion
						while (node.refNucl == '-'){
							node.site[arrayIndex] = '-';
							node = node.next;
							if (node == null){
								throw new RuntimeException("This node shouldnt be null 3");
							}
						}//while
						//assert: node.refNuc == nucQuery
						node.site[arrayIndex] = nuclHit;
						node = node.next;
					}//else
				}//for				

				//The end
				while (node != null){
					node.site[arrayIndex] = '*';
					node = node.next;
				}
			}else{

				//dont understand	
			}
			//if (lineNo > 100000)
			//	System.exit(1);
		}		
		pbOut.close();

		if (arrayIndex >= 0){
			//FIXME
		}
		int status = process.waitFor();		
		outputFull.close();
		outputPart.close();
	}


	static void runBredSeq() throws IOException{
		for (SampleRecord sampleRecord:sampleList){
			String sampleID = sampleRecord.sampleID;
			System.out.println("qsub -b y -cwd -pe smp 8 -N bre"+ sampleID + " -V breseq -n " + sampleID + " -o BreSeq/"+ sampleID + " -j 8 -p -r GeneSequences/pangenes.gbk RAW_READs/" + sampleID + "_1.fastq RAW_READs/" + sampleID + "_2.fastq");
		}
	}

	/**
	 * Represent an alignment from blastn
	 * @author minhduc
	 *
	 */
	static class HSPAlignment{
		String queryID = null, hitID = null;
		int queryFrom, queryTo, hitFrom, hitTo;
		String qSeq = null, hSeq = null;
		boolean sameFrame = true;
		HSPAlignment(){

		}

	}

	static class NodeAlignment {
		static final int LENGTH = 1000;//max these many accepted
		NodeAlignment next = null, prev = null;
		char[] site = new char[LENGTH];//
		int refPos = 0;
		char refNucl = '*';

		public NodeAlignment(char base){
			refNucl = base;
		}	
	}

	static class SampleRecord{
		String  sampleID, genus, species, strain;
		String  readFile1 = null, readFile2 = null;
		boolean isReference = false;
		String  genome = null;
		String  annotation = null;

		SampleRecord(String id){
			sampleID = id;			
		}

		SampleRecord(String id, String gen, String sp, String st){
			sampleID = id;
			genus = gen;
			species = sp;
			strain = st;			
		}

		void setReadFile1(String fileName){
			readFile1 = fileName;			
		}

		void setReadFile2(String fileName){
			readFile2 = fileName;			
		}

		/**
		 * Only complete if either read file(s) or genome file is specified
		 * @return
		 */
		boolean complete(){
			if (readFile1 == null && genome == null)
				return false;

			return true;
		}

		public String toString(){
			StringBuilder ret = new StringBuilder("SAMPLE " + sampleID + "\n");
			if (genus != null)
				ret.append(" GENUS = " + genus + '\n');

			if (species != null)
				ret.append(" GENUS = " + species + '\n');

			if (strain != null)
				ret.append(" STRAIN = " + strain + '\n');

			if (readFile1 != null)
				ret.append(" READ 1 = " + readFile1 + '\n');
			if (readFile2 != null)
				ret.append(" READ 2 = " + readFile2 + '\n');

			if (genome != null)
				ret.append(" GENOME = " + genome + '\n');
			if (annotation != null)
				ret.append(" ANNOTATION = " + annotation + '\n');


			ret.append(" COMPLETE = " + complete() + '\n');

			return ret.toString();
		}

	}


	static class GeneRecord implements Comparable<GeneRecord>{
		String geneID;
		int length;
		int flanks;		
		int isReference = 0;
		String geneDesc;
		String contig = "";
		int start;
		int end;
		char  strand;

		GeneRecord(String gene){
			//GN042|C038:10071-19686:-:19724
			geneID = gene;
			String [] ts = geneID.split(":");
			contig = ts[0];
			strand = ts[2].charAt(0);			
			int gLength = Integer.parseInt(ts[3]);

			ts = ts[1].split("-");

			start = Integer.parseInt(ts[0]);
			end   = Integer.parseInt(ts[1]);

			length = end + 1 -start;


			int toEnd = gLength - end;//convert to distance to end of contig					
			flanks = ((start>PanCoreGenes.FLANK)?PanCoreGenes.FLANK:start) + ((toEnd>PanCoreGenes.FLANK)?PanCoreGenes.FLANK:toEnd);

			SampleRecord sample =  sampleMap.get(contig.split("\\|")[0]);
			if (sample!=null && sample.isReference){
				isReference = 1;
			}


		}
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(GeneRecord o) {
			int comp = Integer.compare(o.isReference, isReference); 

			if (comp==0)
				comp = Integer.compare(o.flanks,  flanks);

			if (comp == 0){
				comp = Integer.compare(o.length, length);
			}
			return comp;
		}

	}	




}
