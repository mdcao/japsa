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

import java.io.IOException;
import java.util.HashMap;

import japsa.seq.Alphabet;
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
	scriptName = "jsa.seq.geneDBBuild", 
	scriptDesc = "Group genes based on their identity and build a database of gene family and their alleles"
	)
public class BuildGeneDatabaseCmd  extends CommandLine{
	//CommandLine cmdLine;
	public BuildGeneDatabaseCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addStdInputFile();
		addString("family", null, "Output family file");
		addString("allele", null, "Output all genes in the database");
		addString("out", null, "out",true);
		addString("prefix", null, "prefix");
		addInt("number", 1000000, "max number of gene per round");
		addBoolean("checkGeneID", true, "Check gene ID");		
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
		CommandLine cmdLine = new BuildGeneDatabaseCmd();
		args = cmdLine.stdParseLine(args);			


		String inputOption = cmdLine.getStringVal("input");
		String familyOption = cmdLine.getStringVal("family");
		String alleleOption = cmdLine.getStringVal("allele");		
		//String listOption = cmdLine.getStringVal("list");
		String outOption = cmdLine.getStringVal("out");
		int number = cmdLine.getIntVal("number");
		boolean checkGeneID = cmdLine.getBooleanVal("checkGeneID");

		String prefix = cmdLine.getStringVal("prefix");
		if (prefix == null)
			prefix = System.currentTimeMillis() + "";

		double thresholdOption = cmdLine.getDoubleVal("threshold");

		BuildGeneDatabase.ratio = thresholdOption;
		BuildGeneDatabase db = new BuildGeneDatabase(prefix);
		SequenceOutputStream sos =  SequenceOutputStream.makeOutputStream(outOption);


		HashMap<String, Sequence> myGenes;
		//if (listOption == null){
		SequenceReader reader = SequenceReader.getReader(inputOption);
		Alphabet.DNA alphabet = Alphabet.DNA();
		Sequence seq;

		myGenes = new HashMap<String, Sequence>();
		while ((seq = reader.nextSequence(alphabet)) != null){									
			if (myGenes.size() >=number){
				Logging.info("BIG TER ");
				HashMap <String, String> mapped = db.addGeneMap(myGenes, checkGeneID);
				for (String key:myGenes.keySet()){
					Sequence keySeq = myGenes.get(key);
					String dbID = mapped.get(key);
					if (dbID != null)
						Logging.info("Added " + key + " as "+ dbID +" G");
					else{
						dbID = db.addGene(myGenes.get(key));
						Logging.info("Added " + key + " as "+ dbID +" B");
					}//else					
					keySeq.setDesc("JSA=" + dbID+";"+keySeq.getDesc());
					keySeq.writeFasta(sos);
				}//for key
				Logging.info("BIG TER END " + db.geneDatabase.size());
				myGenes.clear();				
			}

			myGenes.put(seq.getName(), seq);
		}		
		reader.close();


		HashMap <String, String> mapped = db.addGeneMap(myGenes, checkGeneID);
		for (String key:myGenes.keySet()){
			Sequence keySeq = myGenes.get(key);
			String dbID = mapped.get(key);
			if (dbID != null)
				Logging.info("Added " + key + " as "+ dbID +" G");
			else{
				dbID = db.addGene(myGenes.get(key));
				Logging.info("Added " + key + " as "+ dbID +" B");
			}//else					
			keySeq.setDesc("JSA=" + dbID+";"+keySeq.getDesc());
			keySeq.writeFasta(sos);
		}//for key

		//}else{
		//	myGenes = db.setup(listOption);
		//}

		sos.close();
		db.geneDatabase.write2File(familyOption,false);
		db.geneDatabase.write2File(alleleOption,true);

		db.cleanUp();
	}		
}