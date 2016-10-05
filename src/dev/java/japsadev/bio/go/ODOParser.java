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
 * 20 Jul 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsadev.bio.go;

import japsa.bio.gene.GeneDatabase;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.Logging;
import japsadev.bio.go.GoTerm.GoRelationship;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @author minhduc
 *
 */
public class ODOParser {
	HashMap<String, GoTerm> terms = new HashMap<String,GoTerm>();

	HashMap<String,TypeRelationship> relTypes = 
		new HashMap<String,TypeRelationship>();

	public ODOParser(String fileName) throws IOException{

		//First round to read
		BufferedReader bf = SequenceReader.openFile(fileName);
		String line = "";
		int lineNo = 0;			

		while ((line = bf.readLine())!=null) {
			lineNo ++;
			line = line.trim();
			if (line.equals("[Term]")){
				line = bf.readLine().trim();
				lineNo ++;
				if (!line.startsWith("id: "))
					throw new RuntimeException("Wrong format at line " + lineNo);

				String ID = line.substring(4);
				line = bf.readLine().trim();
				lineNo ++;

				if (!line.startsWith("name: "))
					throw new RuntimeException("Wrong format at line " + lineNo);
				String name = line.substring(5).trim();

				GoTerm currentTerm = new GoTerm(ID);
				currentTerm.name = name;
				terms.put(ID, currentTerm);
				continue;
			}

			if (line.equals("[Typedef]")){
				line = bf.readLine().trim();
				lineNo ++;
				if (!line.startsWith("id: "))
					throw new RuntimeException("Wrong format at line " + lineNo);

				String ID = line.substring(4);
				line = bf.readLine().trim();
				lineNo ++;

				if (!line.startsWith("name: "))
					throw new RuntimeException("Wrong format at line " + lineNo);
				String name = line.substring(5);

				TypeRelationship relationship = new TypeRelationship (ID, name);
				relTypes.put(ID, relationship);				
				continue;				
			}
		}		
		bf.close();		


		//System.out.println(terms.size());
		//System.out.println(relTypes.size());


		//Second round to build relationship		
		bf = SequenceReader.openFile(fileName);
		line = "";
		lineNo = 0;			

		while ((line = bf.readLine())!=null) {
			lineNo ++;
			line = line.trim();


			if (line.equals("[Term]")){
				line = bf.readLine().trim();
				lineNo ++;
				if (!line.startsWith("id: "))
					throw new RuntimeException("Wrong format at line " + lineNo);

				String ID = line.substring(4);

				line = bf.readLine().trim();
				lineNo ++;

				GoTerm currentTerm = terms.get(ID);

				//	System.out.println(currentTerm + " " + lineNo);
				while(true){
					line = bf.readLine();
					lineNo ++;
					if (line == null) 
						break;

					line = line.trim();
					if (line.length() == 0) 
						break;


					if (line.startsWith("is_a: ")){
						String[] toks = line.split(" ",4);
						TypeRelationship relationship =  relTypes.get("is_a");
						GoTerm term = terms.get(toks[1]);
						if (term != null)
							currentTerm.addRelationShip(relationship, term);
						else
							Logging.warn("Term " + toks[1] + " not found at line " + lineNo);

					}
					if (line.startsWith("relationship: ")){
						String[] toks = line.split(" ",4);
						TypeRelationship relationship =  relTypes.get(toks[1]);
						GoTerm term = terms.get(toks[2]);

						if (term != null)
							currentTerm.addRelationShip(relationship, term);
						else
							Logging.warn("Term " + toks[2] + " not found at line " + lineNo);

					}//if								
				}//while
			}//if
		}//while
		bf.close();
	}

	static HashSet<String> conferDrugList = new HashSet<String>();
	static HashSet<String> conferList = new HashSet<String>();
	static StringBuilder sb = new StringBuilder();


	static void confer(GoTerm term){
		for (GoRelationship rel:term.relationship){
			if (rel.relType.ID.equals("confers_resistance_to"))				
				conferList.add(rel.relTerm.name);
			else if (rel.relType.ID.equals("confers_resistance_to_drug"))
				conferDrugList.add(rel.relTerm.name);
			else if (rel.relType.ID.equals("is_a"))
				confer(rel.relTerm);
		}
	}

	static void list(GoTerm term, HashSet<String> relationList){
		for (GoRelationship rel:term.relationship){			
			if (relationList.contains(rel.relType.ID) && (!blackList.contains(rel.relTerm.ID))){
				sb.append(term.ID + "->" + rel.relTerm.ID + "[" + rel.relTerm.name + "];");				
				list(rel.relTerm, relationList);
			}
		}
	}

	static ODOParser odo;		
	//static  HashSet<String> relationList = new HashSet<String>();

	static HashSet<String> blackList = new HashSet<String>();
	static{
		blackList.add("AL513383.1.gene203");
		blackList.add("DQ303918.1.gene1");
		blackList.add("AL513383.1.gene214");

		//CONFUSED GENE
		//blackList.add("KM998962.1.gene1");
		//

		//DQ464881.1.gene2
		//AL513383.1.gene203
		//AY566250.1.gene1
		//AB091338.1.gene1


		//blackList.add("ARO:1000001");//every thing		
		//blackList.add("ARO:3000557");// [antibiotic inactivation enzyme]
		//blackList.add("ARO:3000000");//determinant of antibiotic resistance
		/******************************************************************
		blackList.add("ARO:3000722"); 
		blackList.add("ARO:3000717");
		blackList.add("ARO:3002639");
		blackList.add("ARO:3000744"); //-same as above
		blackList.add("ARO:3000730");
		blackList.add("ARO:3002701");
		blackList.add("ARO:3002703");
		blackList.add("ARO:3000816");
		blackList.add("ARO:3000723");
		blackList.add("ARO:3000498");
		blackList.add("ARO:3000744");
		blackList.add("ARO:3000795");
		blackList.add("ARO:3000796");
		blackList.add("ARO:3002675");
		blackList.add("ARO:3003079");
		blackList.add("ARO:3000730");
		blackList.add("ARO:3003078");
		blackList.add("ARO:3000743");
		blackList.add("ARO:3000822");
		blackList.add("ARO:3002986");
		blackList.add("ARO:3000730");
		blackList.add("ARO:3000717");
		blackList.add("ARO:3000735");
		blackList.add("ARO:3000722");
		blackList.add("ARO:3000734");
		blackList.add("ARO:3000730");
		blackList.add("ARO:3003577");
		blackList.add("ARO:3000299");
		blackList.add("ARO:3000717");
		blackList.add("ARO:3000027");
		blackList.add("ARO:3000730");
		blackList.add("ARO:3000734");
		blackList.add("ARO:3000074");
		blackList.add("ARO:3000838");
		blackList.add("ARO:3000735");
		blackList.add("ARO:3000734");
		blackList.add("ARO:3003301");
		blackList.add("ARO:3000723");
		blackList.add("ARO:3000722");
		blackList.add("ARO:3000735");
		blackList.add("ARO:3000733");
		blackList.add("ARO:3003063");
		blackList.add("ARO:3000781");
		blackList.add("ARO:3000780");
		blackList.add("ARO:3000768");
		blackList.add("ARO:3000782");
		blackList.add("ARO:3000722");
		blackList.add("ARO:3003392");
		blackList.add("ARO:3000730");
		blackList.add("ARO:3001214");
		blackList.add("ARO:3000717");
		blackList.add("ARO:3000737");
		blackList.add("ARO:3003064");
		blackList.add("ARO:3003052");
		blackList.add("ARO:3003053");
		blackList.add("ARO:3003051");
		blackList.add("ARO:3000502");
		blackList.add("ARO:3000502");
		blackList.add("ARO:3000502");
		blackList.add("ARO:3000499");
		blackList.add("ARO:3000499");
		blackList.add("ARO:3000499");
		blackList.add("ARO:3000810");
		blackList.add("ARO:3000811");
		blackList.add("ARO:3003046");
		blackList.add("ARO:3003047");
		blackList.add("ARO:3000804");
		blackList.add("ARO:3000805");
		blackList.add("ARO:3000785");
		blackList.add("ARO:3000254");
		blackList.add("ARO:3000254");
		blackList.add("ARO:3000254");
		blackList.add("ARO:3000254");
		blackList.add("ARO:3000237");
		blackList.add("ARO:3000379");
		blackList.add("ARO:3000379");
		blackList.add("ARO:3000379");
		blackList.add("ARO:3002983");
		blackList.add("ARO:3002983");
		blackList.add("ARO:3003033");
		blackList.add("ARO:3002982");
		blackList.add("ARO:3002982");
		blackList.add("ARO:3003034");
		blackList.add("ARO:3000774");
		blackList.add("ARO:3000774");
		blackList.add("ARO:3000812");
		blackList.add("ARO:3000807");
		blackList.add("ARO:3000806");
		blackList.add("ARO:3000377");
		blackList.add("ARO:3000377");
		blackList.add("ARO:3000794");
		blackList.add("ARO:3000794");
		blackList.add("ARO:3000794");
		blackList.add("ARO:3000792");
		blackList.add("ARO:3000792");
		blackList.add("ARO:3000792");
		blackList.add("ARO:3000793");
		blackList.add("ARO:3000793");
		blackList.add("ARO:3000793");
		blackList.add("ARO:3000800");
		blackList.add("ARO:3000800");
		blackList.add("ARO:3000801");
		blackList.add("ARO:3000801");
		blackList.add("ARO:3000790");
		blackList.add("ARO:3000789");
		blackList.add("ARO:3000791");
		blackList.add("ARO:3000779");
		blackList.add("ARO:3000779");
		blackList.add("ARO:3000779");
		blackList.add("ARO:3000207");
		blackList.add("ARO:3000216");
		blackList.add("ARO:3000809");
		blackList.add("ARO:3000808");
		blackList.add("ARO:3000533");
		blackList.add("ARO:3000378");
		blackList.add("ARO:3000803");
		blackList.add("ARO:3000777");
		blackList.add("ARO:3000778");
		blackList.add("ARO:3000535");
		blackList.add("ARO:3000775");
		blackList.add("ARO:3000775");
		blackList.add("ARO:3003057");
		blackList.add("ARO:3003056");
		blackList.add("ARO:3003055");
		blackList.add("ARO:3003039");
		blackList.add("ARO:3000802");
		blackList.add("ARO:3000776");
		blackList.add("ARO:3000783");
		blackList.add("ARO:3000784");
		blackList.add("ARO:3003010");
		blackList.add("ARO:3003009");
		blackList.add("ARO:3000206");
/******************************************************************/		
		/* 
		 * ARO:3000839
		 * ARO:3000725
		 * ARO:3000815
		 * ARO:3000765
		 * ARO:3000764
		 * ARO:3000508
		 * ARO:3003066
		 * ARO:3003067
		 * ARO:3000725
		 * ARO:3000764
		 * ARO:3000656
		 * ARO:3000502
		 * ARO:3000499
		 * ARO:3000730
		 * 
		 * 
		 * KM998962.1.gene1
		 * 
		 * 
		 */
	}

	public static void resistanceClass(GoTerm term, HashSet<String> allClass, HashSet<String> myClass) throws IOException{
		if (allClass.contains(term.ID))
			myClass.add(term.name);
		else{
			for (GoRelationship rel:term.relationship){			
				if (rel.relType.ID.equals("is_a"))
					resistanceClass(rel.relTerm, allClass, myClass);
			}
		}

	}

	public static HashMap<String, HashSet<String>> checkGenes(ArrayList<Sequence> seqs) throws IOException{
		//Set up known classes		
		HashSet<String> abrGroups = new HashSet<String>();
		abrGroups.add("ARO:3000052");// phenicol resistance gene
		abrGroups.add("ARO:3000102");// fluoroquinolone resistance gene
		abrGroups.add("ARO:3000104");// aminoglycoside resistance gene
		abrGroups.add("ARO:3000129");// beta-lactam resistance gene
		abrGroups.add("ARO:3000240");// streptogramin resistance gene
		abrGroups.add("ARO:3000241");// lincosamide resistance gene
		abrGroups.add("ARO:3000267");// linezolid resistance gene
		abrGroups.add("ARO:3000271");// fosfomycin resistance gene
		abrGroups.add("ARO:3000315");// macrolide resistance gene
		abrGroups.add("ARO:3000362");// mosaic antibiotic resistance gene
		abrGroups.add("ARO:3000383");// rifampin resistance gene
		abrGroups.add("ARO:3000398");// chloramphenicol resistance gene
		abrGroups.add("ARO:3000408");// sulfonamide resistance gene
		abrGroups.add("ARO:3000468");// ethambutol resistance gene
		abrGroups.add("ARO:3000472");// tetracycline resistance gene
		abrGroups.add("ARO:3000477");// aminocoumarin resistance gene
		abrGroups.add("ARO:3000494");// glycopeptide resistance gene
		abrGroups.add("ARO:3000529");// mupirocin resistance gene
		abrGroups.add("ARO:3000751");// peptide antibiotic resistance gene
		abrGroups.add("ARO:3000868");// streptothricin resistance gene
		abrGroups.add("ARO:3001217");// trimethoprim resistance gene
		abrGroups.add("ARO:3001311");// elfamycin resistance gene
		abrGroups.add("ARO:3002984");// polymyxin resistance gene
		abrGroups.add("ARO:3003024");// fusidic acid resistance gene
		abrGroups.add("ARO:3003058");// tunicamycin resistance gene
		abrGroups.add("ARO:3003073");// lipopeptide antibiotic resistance gene
		abrGroups.add("ARO:3003252");// bacitracin resistance gene
		abrGroups.add("ARO:3003432");// isoniazid resistance gene
		abrGroups.add("ARO:3003433");// pyrazinamide resistance gene


		abrGroups.add("ARO:3000004");// class B (metallo-) beta-lactamase
		abrGroups.add("ARO:3000075");// class D beta-lactamase
		abrGroups.add("ARO:3000076");// class C beta-lactamase
		abrGroups.add("ARO:3000078");// class A beta-lactamase	

		HashMap<String, HashSet<String>> ret= new HashMap<String, HashSet<String>>();

		for (Sequence seq : seqs){
			String ID = seq.getName();

			if (blackList.contains(ID))
				continue;

			String desc = seq.getDesc();
			HashSet<String> myGroups = new HashSet<String>();			

			int indexPos = 0;
			while (true){
				indexPos = desc.indexOf("ARO:",indexPos);
				if (indexPos < 0)
					break;

				if (indexPos + 11 < desc.length()){
					String ARO = desc.substring(indexPos, indexPos + 11);
					indexPos += 11;
					if (ARO.equals("ARO:1000001"))
						continue;

					GoTerm term = odo.terms.get(ARO);
					resistanceClass(term,abrGroups, myGroups);					
				}
			}//while
			System.out.print(seq.getName() + " : ");
			for (String group:myGroups)
				System.out.print(group + ";" );

			System.out.println();

			ret.put(seq.getName(), myGroups);
		}
		return ret;
	}

	/**
	 * @param args
	 * @throws IOException 
	 * 
	 */
	public static void main(String[] args) throws IOException {
		odo = new ODOParser("aro.obo");
		//relationList.add("is_a");		
/********************************************************
		ArrayList<Sequence> seqs = SequenceReader.readAll(args[0],Alphabet.DNA());

		HashMap<String, Sequence> allele2Sequence = new HashMap<String, Sequence>();

		for (Sequence seq:seqs){
			allele2Sequence.put(seq.getName(), seq);
		}

		HashMap<String, HashSet<String>> allele2Groups = checkGenes(seqs);

		
		BufferedReader bf = SequenceReader.openFile("togene");
		String line = "";

		HashMap<String, HashSet<String>> gene2Groups = new HashMap<String, HashSet<String>>();
		HashMap<String, String>          gene2dbGeneID = new HashMap<String, String>();

		GeneDatabase geneDB = new GeneDatabase();


		while ( (line = bf.readLine())!=null){
			String [] toks = line.trim().split("\t");			
			String alleleID = toks[0];

			if (!allele2Groups.containsKey(alleleID))
				continue;

			String geneID = ((toks.length > 2)? toks[2]:"None") + "_" + toks[1];

			HashSet<String> agroup = allele2Groups.get(alleleID);
			if (gene2Groups.containsKey(geneID)){
				HashSet<String> ggroup = gene2Groups.get(geneID);	
				if (agroup == null){
					Logging.error(alleleID + " xx  " + geneID );
				}

				if (!ggroup.containsAll(agroup))
					Logging.error( alleleID + " <>  " + geneID );

				if (!agroup.containsAll(ggroup))
					Logging.error(alleleID  + " #  " + geneID );

				GeneDatabase.GeneFamily  dbFamily = geneDB.getFamily(gene2dbGeneID.get(geneID));
				
				if (dbFamily == null)
					Logging.error("Problem finding " + geneID );
				else{
					dbFamily.addSequence(allele2Sequence.get(alleleID));					
				}


			}else{
				gene2Groups.put(geneID, allele2Groups.get(toks[0]));
				Sequence seq = allele2Sequence.get(alleleID);				 
				String famID = geneDB.addNewFamily(seq);
				gene2dbGeneID.put(geneID, famID);		
				String desc = "geneID=" + geneID + ";dg=";

				for (String group:agroup){
					desc += group + ",";
				}
				
				geneDB.getFamily(famID).setDesc(desc);
			}
		}

		geneDB.write2File("F.fasta", false);
		geneDB.write2File("A.fasta", true);
		
		Logging.info(gene2Groups.size() + " group ");


		
		/********************************************************/
		
		for (String relID:odo.relTypes.keySet()){
			System.out.println(relID + "\t" + odo.relTypes.get(relID));
		}

		for (String termID:odo.terms.keySet()){
			GoTerm term = odo.terms.get(termID);			
			for (GoRelationship rel:term.relationship){
				if (rel.relType.ID.equals("is_a"))
					System.out.println(termID + "(" + term.name +")\t" + rel.relTerm.ID + "(" + rel.relTerm.name+")");
			}
			
			conferDrugList.clear();
			conferList.clear();
			sb = new StringBuilder();
			
			odo.confer(term);
			
			System.out.println("   " + term.ID + " :" + term.name + " " + term.desc);
			
			//sodo.list(term, relationList);
			System.out.print("    Confer: ");
			for (String st:conferList){
				System.out.print(st + ";");
			}
			
			System.out.print("\n    Confer drug : ");
			for (String st:conferDrugList){
				System.out.print(st + ";");
			}			
			System.out.print("\n-----------------------------------------------------\n");
			
		}
		/********************************************************/
	}
}
