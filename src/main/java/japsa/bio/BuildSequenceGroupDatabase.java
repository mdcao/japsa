/*
 * Copyright (c) 2017  Minh Duc Cao (minhduc.cao@gmail.com).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the names of the institutions nor the names of the contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/**************************     REVISION HISTORY    **************************
 * 07/09/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio;



import japsa.bio.gene.GeneDatabase;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.seq.Alphabet.DNA;
import japsa.util.HTSUtilities;
//import japsa.util.JapsaTimer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

/**
 * Implement a database of genes
 * TODO: to finalise the design of the database, and move this class out here
 * @author minhduc
 *
 */

public class BuildSequenceGroupDatabase {
    private static final Logger LOG = LoggerFactory.getLogger(BuildSequenceGroupDatabase.class);
	public GeneDatabase geneDatabase;

	String prefix;
	String fFile;
	String gFile;
	String bFile;

	public BuildSequenceGroupDatabase(String p) throws IOException{
		//geneFamilies = new ArrayList<GeneDatabase.GeneFamily>();
		geneDatabase = new  GeneDatabase();
		prefix = p;

		fFile =prefix + "geneFam.fasta";
		gFile =prefix + "gene.fasta";
		bFile =prefix + "gene.sam";

		prepreScript();
		prepreScript2();
	}

	String bwaStr = "bwa mem -t 12 -a -k11 -A1 -B1 -O1 -E1 -L0 -Y ";  
	private void prepreScript() throws IOException{
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(prefix + "runBWA.sh");

		sos.print("bwa index " + fFile);
		sos.println();
		sos.print(bwaStr+ fFile + " " + gFile + "  > " + bFile + " 2> " + prefix + "gene.log");
		sos.println();
		sos.close();
	}

	private void prepreScript2() throws IOException{
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(prefix + "runBWA2.sh");

		sos.print("bwa index " + gFile);
		sos.println();
		sos.print(bwaStr + gFile + " " + gFile + "  > " + bFile + " 2> " + prefix + "gene2.log");
		sos.println();
		sos.close();
	}	

	public void cleanUp() throws IOException, InterruptedException{
		LOG.info("rm -f "+ prefix +"*");
		Process process = Runtime.getRuntime().exec("rm -f "+ prefix +"*");
		process.waitFor();
	}

	/**
	 * Add a collection of genes to the database as a batch
	 * 
	 * @param seqs: new list of sequences to be added
	 * @param checkGeneID: insert based on geneID first
	 * @return A map of what sequences were inserted to what
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public HashMap<String, String> addGeneMap(final Map<String, Sequence> seqs, boolean checkGeneID) throws IOException, InterruptedException{
		HashMap<String, String>  strMap =  new HashMap<String, String> ();

		//0: initialise grouping within the new sequences
		HashMap<String, ArrayList<String>> setMap = new HashMap<String, ArrayList<String>>();
		for (Sequence seq:seqs.values()){
		    //tSet: set of sequences belong to a group. Potentially, every sequence makes up one group
			ArrayList<String> tSet =  new ArrayList<String>();
			tSet.add(seq.getName());
			setMap.put(seq.getName(), tSet);
		}
		//note: setMap: maps each set to an ID (of the first sequence)
        LOG.info("Total " + seqs.size() + " sequences");
		//0.5 Merge based on annotation
		if (checkGeneID){
			HashMap<String, ArrayList<String>> annoMap = new HashMap<String, ArrayList<String>>();
			int getIntrisciID = 0;
			for (Sequence seq:seqs.values()){
				//Try to find the annotation ID to a map
				String desc = seq.getDesc();
				String annoID = null;
				String [] toks = desc.split(";");

				//Gene name is not good, need to use protein
				for (int x = 0; x < toks.length;x++){
					if (toks[x].startsWith("protein_id=")){
						annoID = "ProteinID:" + toks[x].substring(11);
						break;//for x
					}						
				}
				//if not known yet -> search for UniProtKB
				if (annoID == null){
					toks = desc.split("[:;]");
					for (int x = 0; x < toks.length - 1;x++){
						if (toks[x].startsWith("UniProtKB")){
							annoID = "UniProtKB_"+toks[x+1];
							break;//for x
						}						
					}	
				}
				//if still not known yet -> search for CLUSTERS
				if (annoID == null){
					toks = desc.split("[:;]");
					for (int x = 0; x < toks.length - 1;x++){
						if (toks[x].startsWith("CLUSTERS")){
							annoID = "CLUSTERS_"+toks[x+1];
							break;//for x
						}						
					}	
				}
				//if still not known yet -> try once more with Pfam
				if (annoID == null){
					toks = desc.split("[:;]");
					for (int x = 0; x < toks.length - 1;x++){
						if (toks[x].startsWith("Pfam")){
							annoID = "Pfam_"+toks[x+1];
							break;//for x
						}//if						
					}//for x	
				}//if


				//got annotation ID
				if (annoID != null){
					getIntrisciID ++;
					annoID = annoID.replaceAll("/", "_");
					ArrayList<String> list = annoMap.get(annoID);
					if (list == null){
						list = new ArrayList<String>();
						annoMap.put(annoID, list);
					}
					list.add(seq.getName());
				}else{
					//LOG.info("NOT found annoID " + seq.getName());
				}
			}//for

            LOG.info(" Step 0.5  " + annoMap.size() + " from " + getIntrisciID);
			for (String annoID:annoMap.keySet()){
				ArrayList<String> list = annoMap.get(annoID);
				String firstName = list.get(0);
				ArrayList<String> setGroup = setMap.get(firstName);

				for (int i = 1;i < list.size();i++){
					String tName = list.get(i);
					setGroup.add(tName);
					setMap.put(tName, setGroup);
				}					
			}//for
		}

		boolean merge = true;
		int iteration = 0;
		while (merge){
			//1. Group these sequences first
            LOG.info(" --Start of " + iteration + "  " + new Date());
			merge = false;
			iteration ++;
			SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(gFile);			
			LOG.info("Iteration " + iteration + " size = " + setMap.size() + "  " + new Date());
			int countI = 0;

			//write the reps of temp groups to file anc cosider them as `reads' as well as refs
			for (String key:setMap.keySet()){
				ArrayList<String> tSet = setMap.get(key);
				if (tSet.get(0) == key){
					countI ++;
					seqs.get(key).writeFasta(sos);
				}
			}
			sos.close();

			LOG.info("Iteration " + iteration + " of " + countI + " sequences " + new Date());
			//JapsaTimer.systemInfo();

			Process process = Runtime.getRuntime().exec("bash " + prefix + "runBWA2.sh");
			process.waitFor();

			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			SamReader samReader = SamReaderFactory.makeDefault().open(new File(bFile));
			SAMRecordIterator samIter = samReader.iterator();

			while (samIter.hasNext()){
				SAMRecord sam = samIter.next();

				if (sam.getReadUnmappedFlag())
					continue;
				//TODO: Remember that need to have softclip (with -Y turn on)

				String readName = sam.getReadName();
				String refName = sam.getReferenceName();

				ArrayList<String> readSet = setMap.get(readName);
				ArrayList<String> refSet = setMap.get(refName);

				if (readSet == refSet){
					//these two have been from group
					continue;
				}

				Sequence readSeq = seqs.get(readName);
				Sequence refSeq = seqs.get(refName);

				//now we find the respresentations from groups are similar, and hence join two groups together
				if (isSimilar(refSeq, readSeq, sam)){
					merge = true;
					refSet.addAll(readSet);
					for (String key:readSet){//merge every one in
						setMap.put(key, refSet);//both points to the same set
					}//basically, no more readSet
					readSet.clear();//help the GC
				}
			}
			samIter.close();
			samReader.close();
			//JapsaTimer.systemInfo();
		}

		LOG.info(" 2 Grouping " + geneDatabase.size());
		//JapsaTimer.systemInfo();

		//2.Try to add to the existing groups
		int G = 0, GG = 0;
		if (geneDatabase.size() > 0){			
			geneDatabase.write2File(fFile, false);

			//creat a map of current reps

            SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(fFile);
            HashMap<String, Sequence> repMap = new HashMap<String, Sequence>();
            for (GeneDatabase.GeneFamily family:geneDatabase){
                Sequence rep = family.represetationSequence();
                repMap.put(rep.getName(),rep);
            }
            sos.close();


			//Run bwa
			//LOG.info("Running bwa for " + seq.getName() + " " + geneFamilies.size() + " family");
			Process process = Runtime.getRuntime().exec("bash " + prefix + "runBWA.sh");
			process.waitFor();

			//Read the sam file
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			SamReader samReader = SamReaderFactory.makeDefault().open(new File(bFile));	

			SAMRecordIterator samIter = samReader.iterator();
			String geneID = null;		
			while (samIter.hasNext()){
				SAMRecord sam = samIter.next();

				if (sam.getReadUnmappedFlag())
					continue;
				//TODO: Remember that need to have softclip (with -Y turn on)

				String readName = sam.getReadName();
				if (strMap.containsKey(readName))
					continue;

				Sequence readSeq = seqs.get(readName);
				if (readSeq == null){
					LOG.error("ERROR 4: sequence " + readName + " not found!");
				}
				String refName = sam.getReferenceName();
				GeneDatabase.GeneFamily family = geneDatabase.getFamily(refName);
				if (family == null){
					LOG.error("ERROR 5: family " + refName + " not found!");
					continue;
				}

				Sequence refSeq = family.represetationSequence();
				if (refSeq == null){
					LOG.error("ERROR 6: rep for family " + refName + " not found!");
				}
				
				if (isSimilar(refSeq, readSeq, sam)){
					ArrayList<String> readSet = setMap.get(readName);

					for (String key:readSet){
						//if (strMap.containsKey(key)){
						//	LOG.error("ERROR 1 : " + key);
						//}
						Sequence keySeq = seqs.get(key);
						geneID = family.addSequence(keySeq);
						strMap.put(key, geneID);
						G++;						
					}
					readSet.clear();
				}
				//LOG.info(" ADD a G " + G + " " + (new Date()));
			}
			samIter.close();
			samReader.close();
		}


		LOG.info(" 3 Grouping " + geneDatabase.size() + " " + seqs.size());
		//JapsaTimer.systemInfo();

		for (Sequence seq:seqs.values()){
			String seqName = seq.getName();
			ArrayList<String> tSet =  setMap.get(seqName);
			if (tSet.isEmpty()){
				if (!strMap.containsKey(seqName)){
                    LOG.error("ERROR 2 : " + seqName);
				}
				continue;//added
			}
			//tSet is not empty
			GeneDatabase.GeneFamily family = null;
			LOG.info(" ADDing GGs " + GG + " of size " + tSet.size() + " " + (new Date()));
			for (String key:tSet){
				Sequence keySeq = seqs.get(key);
				if (family == null){
					String geneID = geneDatabase.addNewFamily(keySeq);
					family = geneDatabase.getFamily(geneID);
					strMap.put(key, geneID);					
				}else{
					String geneID = family.addSequence(keySeq);					
					strMap.put(key, geneID);
				}				
				GG ++;				
			}
			tSet.clear();
		}

        LOG.info("Manage to add " + G + " and " + GG + " " + new Date());
		return strMap;
	}

	public String addGene(Sequence seq) throws IOException, InterruptedException{
		//1.Create first family is not done so already
		if (geneDatabase.size() == 0){
			return geneDatabase.addNewFamily(seq);			
			//GeneDatabase.GeneFamily geneFamily = new GeneDatabase.GeneFamily(0);
			//geneFamilies.add(geneFamily);
			//return geneFamily.addSequence(seq);
		}
		//Prepare index and sequence
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(prefix + "geneFam.fasta");

		for (GeneDatabase.GeneFamily family:geneDatabase){
			family.represetationSequence().writeFasta(sos);
		}
		sos.close();

		sos = SequenceOutputStream.makeOutputStream(prefix + "gene.fasta");
		seq.writeFasta(sos);
		sos.close();

		//Run bwa
		//LOG.info("Running bwa for " + seq.getName() + " " + geneFamilies.size() + " family");
		Process process = Runtime.getRuntime().exec("bash " + prefix + "runBWA.sh");
		process.waitFor();

		//Read the sam file
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(bFile));	

		SAMRecordIterator samIter = samReader.iterator();
		String geneID = null;		
		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();

			if (sam.getReadUnmappedFlag())
				continue;
			//TODO: Remember that need to have softclip (with -Y turn on)
			String refName = sam.getReferenceName();
			GeneDatabase.GeneFamily family = geneDatabase.getFamily(refName);
			if (family == null){
                LOG.error("Check for problem : family " + refName + " not found!");
				continue;
			}
			Sequence refSeq = family.represetationSequence();
			
			if (isSimilar(refSeq, seq, sam)){
				geneID = family.addSequence(seq);
				break;
			}
		}
		samIter.close();
		samReader.close();

		if (geneID != null) return geneID;
		return geneDatabase.addNewFamily(seq);	
	}

	//TODO: no good public
	public static double ratio = 0.9;
	public static double coverageTheshold = 0.9;

	static boolean isSimilar(Sequence refSeq, Sequence readSeq, SAMRecord sam){
		japsa.util.HTSUtilities.IdentityProfile 
		//in the other strand -> need to reverse
		profile = sam.getReadNegativeStrandFlag()?HTSUtilities.identity(refSeq, DNA.complement(readSeq), sam)
			:HTSUtilities.identity(refSeq, readSeq, sam); 

		double m = profile.match / ratio;
		return (m > profile.refBase
                && m > profile.readBase
                && profile.refBase > refSeq.length() * coverageTheshold
                && profile.readBase > readSeq.length() * coverageTheshold
        );
	}

	/**
	 * Set up with a list of strains
	 * @param file
	 * @return
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public HashMap<String, Sequence> setup(String file) throws IOException, InterruptedException{
		BufferedReader bf = SequenceReader.openFile(file);
		String line = null;

		HashMap<String, Sequence> myGenes = new HashMap<String, Sequence>();
		while ( (line = bf.readLine())!=null){
			if (line.startsWith("#"))
				continue;

			String [] toks = line.trim().split("\t");
			if (toks.length < 7)
				continue;					

			String strain = toks[0].replaceAll(" ","");
			String gffFile = toks[1];
			String species = toks[2].replaceAll(" ","");			
			String n50 = toks[4];

			if (Integer.parseInt(n50) < 100000){
                LOG.info("Strain " + strain + ": skipped as N50=" + toks[4]);
				continue;
			}


			FileInputStream gffIn = new FileInputStream(gffFile);
			ArrayList<JapsaAnnotation>
			annoMap = JapsaAnnotation.readMGFF(gffIn,0,0,"CDS");
			gffIn.close();

            LOG.info("Genome " + toks[3] + " " + strain);
			LOG.info("There are " + annoMap.size()+ " annotations here");
			for (JapsaAnnotation anno:annoMap){				
				for (int i = 0; i < anno.numFeatures(); i++){
					//totGenes ++;
					String geneID = null;
					JapsaFeature feature = anno.getFeature(i);
					String desc = feature.getDesc();

					//search for geneID
					toks = desc.split(";");
					for (int x = 0; x < toks.length;x++){
						if (toks[x].startsWith("ID=")){
							geneID = toks[x].substring(3);
							break;//for x
						}						
					}
					if (geneID == null){
						//LOG.error("ERROR = " + desc );
						continue;
					}
					//totGenes ++;
					geneID = species + "_" + strain + "_" + geneID;
					Sequence geneSeq = anno.getSequence().subSequence(feature.getStart() - 1, feature.getEnd());
					geneSeq.setName(geneID);
					geneSeq.setDesc(desc);	
					//geneSeq.writeFasta(sos);

					myGenes.put(geneID, geneSeq);
					//geneID = addGene(geneSeq);			
					//LOG.info("Added " + geneSeq.getName() + " as "+ geneID);
				}//for i
				/*************************************************************
				LOG.info("Trying to add " + anno.getAnnotationID() + " of " + anno.numFeatures() + " genes");
				/*************************************************************/

			}//for anno

		}//while		
		return myGenes;

	}	
}