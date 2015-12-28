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

/*                           Revision History                                
 * 11/01/2012 - Minh Duc Cao: Revised 
 * 01/01/2013 - Minh Duc Cao, revised                                       
 ****************************************************************************/

package japsa.tools.seq;

import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(scriptName = "jsa.seq.gff2fasta",
scriptDesc = "Extract sequences from a gff annotation")
public class ExtractGeneSequenceCmd extends CommandLine{	
	public ExtractGeneSequenceCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("sequence", null, "The sequence (whole chromosome)",true);
		addString("gff", null, "Annotation file in gff format",true);		
		addString("type", "gene", "types of features to be extracted (all, gene, CDS etc)");
		addInt("flank", 0, "Size of flanking regions");
		addStdOutputFile();	

		addStdHelp();		
	} 


	public static void main(String[] args) throws IOException {
		CommandLine cmdLine = new ExtractGeneSequenceCmd();	
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/		
		String sequence = cmdLine.getStringVal("sequence");		
		String gff = cmdLine.getStringVal("gff");
		String type = cmdLine.getStringVal("type");		
		String output = cmdLine.getStringVal("output");
		int flank = cmdLine.getIntVal("flank");	

		/**********************************************************************/
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(output);
		extractGenes(sequence, gff, type, flank, out);
		out.close();

	}

	public static void extractGenes(String sequence, String gff, String type, int flank, SequenceOutputStream out) throws IOException{

		//Read annotation, without upstream and downstream
		FileInputStream aReader = new FileInputStream(gff);		
		ArrayList<JapsaAnnotation> annos = JapsaAnnotation.readMGFF(aReader,0,0,type);
		aReader.close();

		Logging.info("Read " + annos.size());

		SequenceReader reader = SequenceReader.getReader(sequence);
		Sequence seq = reader.nextSequence(Alphabet.DNA());

		for (JapsaAnnotation anno:annos){	
			//Logging.info("Read anno " + anno.numFeatures());
			while (seq != null && !seq.getName().equals(anno.getAnnotationID())){
				seq = reader.nextSequence(Alphabet.DNA());
			}
			if (seq == null){
				Logging.error("Sequence " + anno.getAnnotationID() + " not found");
				reader.close();
				out.close();				
				System.exit(1);
			}

			for (JapsaFeature feature:anno.getFeatureList()){
				int start = feature.getStart() - 1;//note: convert 1-index, inclusive to 0-index inclusive 
				int end = feature.getEnd();//note: 1-index, inclusive == 0-index exclusive
				Sequence //featureSeq = new Sequence(Alphabet.DNA(), end - start + 2 * flank);
				featureSeq = seq.subsequenceWithFlank(start, end, flank);

				if (feature.getStrand() == '-')
					featureSeq = Alphabet.DNA.complement(featureSeq);

				featureSeq.setName(seq.getName() + ":" + (start + 1) + "-" + end +":" + feature.getStrand() + ":" + seq.length());
				//featureSeq.setDesc("Type=" + feature.getType() + ";" + feature.getDesc());
				featureSeq.setDesc(feature.getDesc());


				featureSeq.writeFasta(out);
			}//for feature	
		}//for anno
		reader.close();
	}

	public static HashMap<String, Sequence>  extractGenes(String sequence, String gff, String type, int flank) throws IOException{
		HashMap<String, Sequence> seqList = new HashMap<String, Sequence>();

		//Read annotation, without upstream and downstream
		FileInputStream aReader = new FileInputStream(gff);		
		ArrayList<JapsaAnnotation> annos = JapsaAnnotation.readMGFF(aReader,0,0,type);
		aReader.close();

		Logging.info("Read " + annos.size());

		SequenceReader reader = SequenceReader.getReader(sequence);
		Sequence seq = reader.nextSequence(Alphabet.DNA());

		for (JapsaAnnotation anno:annos){	
			//Logging.info("Read anno " + anno.numFeatures());
			while (seq != null && !seq.getName().equals(anno.getAnnotationID())){
				seq = reader.nextSequence(Alphabet.DNA());
			}

			if (seq == null){
				Logging.error("Sequence " + anno.getAnnotationID() + " not found");
				reader.close();							
				System.exit(1);
			}

			for (JapsaFeature feature:anno.getFeatureList()){
				int start = feature.getStart() - 1;//note: convert 1-index, inclusive to 0-index inclusive 
				int end = feature.getEnd();//note: 1-index, inclusive == 0-index exclusive
				Sequence //featureSeq = new Sequence(Alphabet.DNA(), end - start + 2 * flank);
				featureSeq = seq.subsequenceWithFlank(start, end, flank);

				//for (int i = 0; i < featureSeq.length();i++){
				//	int j = start - flank + i; 
				//	if (j < 0 || j>= seq.length())
				//		featureSeq.setSymbol(i, Alphabet.DNA.N);
				//	else 
				//		featureSeq.setSymbol(i, seq.getBase(j));
				//}//for

				if (feature.getStrand() == '-')
					featureSeq = Alphabet.DNA.complement(featureSeq);

				featureSeq.setName(seq.getName() + ":" + (start + 1) + "-" + end +":" + feature.getStrand() + ":" + seq.length());
				featureSeq.setDesc("Type=" + feature.getType() + ";" + feature.getDesc());
				seqList.put(featureSeq.getName(),featureSeq);
			}//for feature	
		}//for anno
		reader.close();
		return seqList;
	}
}


/*RST*
-------------------------------------------
 *jsa.seq.gff2fasta*: Extract gene sequences 
-------------------------------------------

 *jsa.seq.gff2fasta* extract the functional sequences (genes, CDS, etc) from
a gff file and a sequence file.

<usage> 

 *RST*/
