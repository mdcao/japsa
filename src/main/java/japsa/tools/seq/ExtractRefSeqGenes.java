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

/**************************     REVISION HISTORY    **************************
 * 07/09/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.tools.seq;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import japsa.bio.BuildGeneDatabase;
import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;


/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.seq.refseq2genesDb", 
	scriptDesc = "Extract gene sequences from refseq anotation and group them"
	)
public class ExtractRefSeqGenes  extends CommandLine{
	//CommandLine cmdLine;
	public ExtractRefSeqGenes(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());						

		addString("db", "-", "name of db file");		
		addString("gene", "genes.ffn", "name of gene file");
		addString("annoType", "RS", "Annotation type, refseq or prokka");

		addString("family", null, "Output family file");
		addString("allele", null, "Output all genes in the database");
		addString("out", null, "out",true);
		addString("prefix", null, "prefix");		
		addDouble("threshold", 0.9, "Threshold identity");

		addStdHelp();
	}

	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws Exception 
	 * @throws OutOfMemoryError 
	 */
	public static void main(String[] args) throws IOException, InterruptedException{		 		
		CommandLine cmdLine = new ExtractRefSeqGenes();
		args = cmdLine.stdParseLine(args);		
		String db   = cmdLine.getStringVal("db");
		String gene   = cmdLine.getStringVal("gene");
		String annoType   = cmdLine.getStringVal("annoType");


		String familyOption = cmdLine.getStringVal("family");
		String alleleOption = cmdLine.getStringVal("allele");		
		String outOption = cmdLine.getStringVal("out");		

		String prefix = cmdLine.getStringVal("prefix");
		if (prefix == null)
			prefix = System.currentTimeMillis() + "";

		double thresholdOption = cmdLine.getDoubleVal("threshold");
		BuildGeneDatabase.ratio = thresholdOption;

		geneOS = SequenceOutputStream.makeOutputStream(gene);
		processDB(db, annoType, familyOption, alleleOption, outOption, prefix);
		geneOS.close();
	}


	static SequenceOutputStream geneOS;//, proOS;	 

	private static void processDB(String dbFile, String annoType, String family, String allele, String outOption, String prefix) throws IOException, InterruptedException{		
		BufferedReader bf = SequenceReader.openFile(dbFile);		
		HashSet<String> organismSet = new HashSet<String>();
		HashSet<String> stSet = new HashSet<String>();

		BuildGeneDatabase geneDB = new BuildGeneDatabase(prefix);
		SequenceOutputStream sos =  SequenceOutputStream.makeOutputStream(outOption);

		String line = "";
		while ( (line = bf.readLine())!=null){
			if (line.startsWith("#"))
				continue;

			String [] toks = line.trim().split("\t");
			String strainID = toks[4];

			//if (toks[2].equals(toks[1]) && toks[3].equals("")){
			//	Logging.info(strainID + " Ignored because of no strain information");
			//	continue;				
			//}

			String fnaFile = toks[5];
			String gffFile = fnaFile.replace(".fna.gz",".gff.gz");										

			String species = toks[1];
			String strain = toks[3];

			double n50 = Double.parseDouble(toks[7]);

			if (n50 < 100000){
				Logging.info(strainID + " Ignored because of low n50 " + n50);
				continue;//while
			}


			String organismName = toks[2];			
			if (!organismName.equals(toks[1] + " " + toks[3])){
				organismName += " " + toks[3];
			}

			String strainName = organismName.replaceAll(" ", "_");
			strainName = strainName.replaceAll("/", "_");
			strainName = strainName.replaceAll("'", "_");
			strainName = strainName.replaceAll("\"", "_");
			strainName = strainName.replaceAll(";", "_");	
			strainName = strainName.replaceAll(":", "_");
			strainName = strainName.replaceAll("__*", "_");//Make sure no double hyphen

			if (!organismSet.add(organismName)){
				Logging.info(strainID + " Ignored because of dupNAME\t" + organismName + "\t" + n50);
				continue;
			}					


			if (toks.length >= 14){
				if (!toks[13].equals("0")){
					Logging.info(strainID + " Ignored because of not good ST " + toks[13]);
					continue;//while
				}
				String ST = toks[12];
				if (!stSet.add(toks[11])){
					Logging.info(strainID + " Ignored because of dupST\t" + organismName + "\t" + n50 + " ST_" + ST + " (" + toks[11] +")");
					continue;
				}				
				
				strainName += "_ST" + ST;
			}
			
			/*****************************************************/
			HashMap<String, Sequence> myGenes = "RS".equals(annoType)?readStrainRefSeq (strainName, gffFile, fnaFile):readStrainProkka (strainName, toks[6]);
			HashMap <String, String> mapped = geneDB.addGeneMap(myGenes, true);
			for (String key:myGenes.keySet()){
				Sequence keySeq = myGenes.get(key);
				String dbID = mapped.get(key);
				if (dbID != null)
					Logging.info("Added " + key + " as "+ dbID +" G");
				else{
					dbID = geneDB.addGene(myGenes.get(key));
					Logging.info("Added " + key + " as "+ dbID +" B");
				}//else					
				keySeq.setDesc("JSA=" + dbID+";"+keySeq.getDesc());
				keySeq.writeFasta(sos);
			}//for key
			Logging.info("BIG END " + strainName + " " + geneDB.geneDatabase.size());			
			/*****************************************************/			
			System.out.println(strainID + "\t" + strainName + "\t" + species + "\t" + strain + "\t" + toks[8] + "\t" + n50);
		}

		bf.close();

		sos.close();
		geneDB.geneDatabase.write2File(family,false);
		geneDB.geneDatabase.write2File(allele,true);
		geneDB.cleanUp();
	}

	private static HashMap<String, Sequence> readStrainProkka(String strainID, String gff) throws IOException{		
		/*******************   Read in protein  ***************************/
		//read in
		InputStream in  = new FileInputStream(gff);
		ArrayList<JapsaAnnotation> annoGFF = JapsaAnnotation.readMGFF(in,0,0, "CDS");
		in.close();

		HashMap<String, Sequence> myGenes = new HashMap<String, Sequence>();

		for (JapsaAnnotation anno:annoGFF){
			Sequence seq = anno.getSequence();
			for (JapsaFeature feature:anno.getFeatureList()){									
				String desc = feature.getDesc();
				String geneID = feature.getID();

				geneID = strainID + "__" + geneID;

				int start = feature.getStart() - 1;
				int end = feature.getEnd();				

				Sequence geneSeq = seq.subSequence(start, end);

				if (feature.getStrand() == '-')
					geneSeq = Alphabet.DNA.complement(geneSeq);

				geneSeq.setName(geneID);
				geneSeq.setDesc(desc);

				geneSeq.writeFasta(geneOS);
				myGenes.put(geneID, geneSeq);
			}//for				
		}//for
		return myGenes;
	}


	private static HashMap<String, Sequence> readStrainRefSeq (String strainID, String gff, String fnaFile) throws IOException{		
		//read in
		InputStream in  = new FileInputStream(gff);
		ArrayList<JapsaAnnotation> annoGFF = JapsaAnnotation.readMGFF(in,0,0, "CDS");
		in.close();

		HashMap<String, Sequence> myGenes = new HashMap<String, Sequence>();

		ArrayList<Sequence> seqs = SequenceReader.readAll(fnaFile, Alphabet.DNA());

		for (JapsaAnnotation anno:annoGFF){
			Sequence seq = anno.getSequence();
			if (seq == null){
				for (Sequence s:seqs){
					if (s.getName().equals(anno.getAnnotationID())){
						seq = s;
						break;
					}
				}
			}

			if (seq == null){
				Logging.info("ERROR: not found sequence for " + anno.getAnnotationID());
				continue;
			}

			for (JapsaFeature feature:anno.getFeatureList()){									
				String desc = feature.getDesc();
				String geneID = feature.getID();

				//Sequence pro = proteinMap.get(geneID);
				geneID = strainID + "__" + geneID;
				//totGenes ++;


				int start = feature.getStart() - 1;
				int end = feature.getEnd();				

				Sequence geneSeq = seq.subSequence(start, end);

				if (feature.getStrand() == '-')
					geneSeq = Alphabet.DNA.complement(geneSeq);

				geneSeq.setName(geneID);
				geneSeq.setDesc(desc);

				geneSeq.writeFasta(geneOS);
				myGenes.put(geneID, geneSeq);
			}				
		}
		return myGenes;
	}


}

