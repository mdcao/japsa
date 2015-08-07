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
 * 21/06/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.seq.tools;

import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;



/**
 * FIXME: Need testing
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
@Deployable(scriptName = "jsa.seq.bed2jsa",
            scriptDesc = "Convert gene annotation from bed format to Japsa format")
public class Bed2Japsa {
	public static void main(String[] args) throws IOException {
		
		/*********************** Setting up script ****************************/
		Deployable annotation = Bed2Japsa.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options]",
				annotation.scriptDesc());		
		/**********************************************************************/		
		cmdLine.addStdInputFile();		
		cmdLine.addString("output", "-",
				"Name of the output file ( - for standard output)");
		
		cmdLine.addDouble("maxOverlap", 0.3, "Maximum overlap of gene");
		cmdLine.addInt("upstream", 2000, "Length of upstream region");
		cmdLine.addInt("downstream", 2000, "Length of downstream region");
		
		cmdLine.addString("sequence", null,
				"Name of the sequence file");				

		args = cmdLine.stdParseLine_old(args);
		/**********************************************************************/
		
		String inFile = cmdLine.getStringVal("input");		
		String outFile = cmdLine.getStringVal("output");
		String seqFile = cmdLine.getStringVal("sequence");
		
		if ("-".equals(seqFile)){
			System.err.println("ERROR: A sequence file must be suplied (not -)\n" + cmdLine.usageMessage());
			System.exit(-1);
		}		
		
		BufferedReader in = SequenceReader.openFile(inFile);
		ArrayList<JapsaAnnotation> annos = readBED(in, cmdLine.getDoubleVal("maxOverlap"), cmdLine.getIntVal("upstream"),cmdLine.getIntVal("downstream"));
		
		SequenceReader reader = null;
		if (seqFile != null){
			reader = SequenceReader.getReader(seqFile);
		}
		
		SequenceOutputStream outStream = SequenceOutputStream.makeOutputStream(outFile);
		for (int i = 0; i < annos.size();i++){
			JapsaAnnotation anno = annos.get(i); 
			if (reader == null)
				anno.writeAnnotation(outStream);
			else{
				anno.sortFeatures();			
				JapsaAnnotation.write(reader.nextSequence(Alphabet.DNA16()), anno, outStream);
			}
		}
		
		outStream.close();
	}

	

	/**
	 * Read in gene annotations of chromosomes in a bed file.
	 * 
	 * Pre-requisite: The genes in a chromosome are sorted by position
	 * @param in: BufferedReader to read in
	 * @throws IOException 
	 *  
	 */
	static	ArrayList<JapsaAnnotation> readBED(BufferedReader in, double allowOverlap, int upStream, int downStream) throws IOException{

		int noIgnore = 0;
		ArrayList<JapsaAnnotation> annos = new ArrayList<JapsaAnnotation>();
		JapsaAnnotation anno = null;
		JapsaFeature lastFeature = null;
		String line = null;		
		while( (line = in.readLine()) != null){
			String [] toks = line.trim().split("\t");

			if (anno == null || !anno.getAnnotationID().equals(toks[0])){
				anno = new JapsaAnnotation();
				anno.setAnnotationID(toks[0]);
				annos.add(anno);

				lastFeature = null;//the last feature
			}

			int start = Integer.parseInt(toks[1]) + 1;
			int end = Integer.parseInt(toks[2]);

			int cdsStart = Integer.parseInt(toks[6]) + 1;
			int cdsEnd   = Integer.parseInt(toks[7]);//the end point is exclusive hence not plus 1
			
			boolean codingGene = (cdsEnd > cdsStart);//is it a coding gene

			if (lastFeature != null && (lastFeature.getEnd() + 1.0 - start) /(end - start + 1.0) > allowOverlap){
				System.err.println("Gene " + toks[3] + "["+start +","+end+"] ignored because of overlaping with " + 
						lastFeature.getID()+"["+lastFeature.getStart() +"," + lastFeature.getEnd()+"] of " +((lastFeature.getEnd() + 1.0 - start) /(end - start + 1.0)));
				noIgnore ++;
				continue;
			}
			//Admit this 
			lastFeature = new JapsaFeature(start, end, "gene", toks[3], toks[5].charAt(0), toks[0]);
			lastFeature.setScore(Double.parseDouble(toks[4]));
			
			if (!codingGene)
				lastFeature.addDesc("Non-coding gene");

			//now get the various parts of the gene

			
			//int noExons = Integer.parseInt(toks[9]); 
					
			if (lastFeature.getStrand() != '-'){//positive strand
				
				//0. Upstream if needed
				if (upStream > 0)
					anno.add(new JapsaFeature(start - upStream, start - 1, 
							"upstream", lastFeature.getID()+"#U",  lastFeature.getStrand(), lastFeature.getID()));
				
				//1. add feature
				anno.add(lastFeature);
				
				String [] exonLens  = toks[10].split(",");
				String [] exonStart = toks[11].split(",");				
				
				
				//do some checking
				if (Integer.parseInt(exonStart[0]) != 0){
					System.err.println("###GENE" + lastFeature.getID()+"["+lastFeature.getStart() +"," + lastFeature.getEnd()+"] with non-zero first exon : " + exonStart[0]);
				}				
				
				int my_end = Integer.MAX_VALUE;// the end of the last exon
				for (int i = 0; i < exonStart.length; i++){
					//the start of an exon
					int my_start = Integer.parseInt(exonStart[i]) + start;
					if (my_end < my_start - 1){
						anno.add(new JapsaFeature(my_end + 1, my_start - 1, 
								"intron", lastFeature.getID()+"#I",  lastFeature.getStrand(), lastFeature.getID()));
					}					
					my_end = my_start + Integer.parseInt(exonLens[i]) - 1;
					
					if (my_start < cdsStart){
						if (my_end > cdsStart ){
							anno.add(new JapsaFeature(my_start, cdsStart - 1, 
									codingGene?"5UTR":"UTR", lastFeature.getID()+ (codingGene?"#5":"#U"),  lastFeature.getStrand(), lastFeature.getID()));
							my_start = cdsStart;							
						}else{ //cdsStart >= my_end
							anno.add(new JapsaFeature(my_start, my_end, 
									codingGene?"5UTR":"UTR", lastFeature.getID()+ (codingGene?"#5":"#U"),  lastFeature.getStrand(), lastFeature.getID()));
							continue;//for
						}						
					}					
					
					if (my_end <= cdsEnd){
						anno.add(new JapsaFeature(my_start, my_end, 
								"CDS", lastFeature.getID()+"#C",  lastFeature.getStrand(), lastFeature.getID()));
						continue;//for						
					}
					
					//assert my_end > cdsEnd					
					if (my_start <= cdsEnd){
						anno.add(new JapsaFeature(my_start, cdsEnd, 
								"CDS", lastFeature.getID()+"#C",  lastFeature.getStrand(), lastFeature.getID()));
						my_start = cdsEnd + 1;
					}
					
					if (my_start <= my_end)
						anno.add(new JapsaFeature(my_start, my_end, 
							codingGene?"3UTR":"UTR", lastFeature.getID()+ (codingGene?"#3":"#U"), lastFeature.getStrand(), lastFeature.getID()));
				}//for
				if (downStream > 0){			
						anno.add(new JapsaFeature(end + 1, end + downStream, 
								"downstream", lastFeature.getID()+"#D",  lastFeature.getStrand(), lastFeature.getID()));
				}					
				
			}else{//negative strand

				//0. Upstream if needed
				if (downStream > 0)
					anno.add(new JapsaFeature(start - downStream, start - 1, 
							"downstream", lastFeature.getID()+"#D",  lastFeature.getStrand(), lastFeature.getID()));
				
				//1. add feature
				anno.add(lastFeature);
				
				String [] exonLens  = toks[10].split(",");
				String [] exonStart = toks[11].split(",");				
				
				
				//do some checking
				if (Integer.parseInt(exonStart[0]) != 0){
					System.err.println("###GENE" + lastFeature.getID()+"["+lastFeature.getStart() +"," + lastFeature.getEnd()+"] with non-zero first exon : " + exonStart[0]);
				}				
				
				int my_end = Integer.MAX_VALUE;// the end of the last exon
				for (int i = 0; i < exonStart.length; i++){
					//the start of an exon
					int my_start = Integer.parseInt(exonStart[i]) + start;
					if (my_end < my_start - 1){
						anno.add(new JapsaFeature(my_end + 1, my_start - 1, 
								"intron", lastFeature.getID()+"#I",  lastFeature.getStrand(), lastFeature.getID()));
					}					
					my_end = my_start + Integer.parseInt(exonLens[i]) - 1;
					
					if (my_start < cdsStart){
						if (my_end > cdsStart ){
							anno.add(new JapsaFeature(my_start, cdsStart - 1, 
									codingGene?"3UTR":"UTR", lastFeature.getID()+ (codingGene?"#3":"#U"),  lastFeature.getStrand(), lastFeature.getID()));
							my_start = cdsStart;							
						}else{ //cdsStart >= my_end
							anno.add(new JapsaFeature(my_start, my_end, 
									codingGene?"3UTR":"UTR", lastFeature.getID()+ (codingGene?"#3":"#U"),  lastFeature.getStrand(), lastFeature.getID()));
							continue;//for
						}						
					}					
					
					if (my_end <= cdsEnd){
						anno.add(new JapsaFeature(my_start, my_end, 
								"CDS", lastFeature.getID()+"#C",  lastFeature.getStrand(), lastFeature.getID()));
						continue;//for						
					}
					
					//assert my_end > cdsEnd					
					if (my_start <= cdsEnd){
						anno.add(new JapsaFeature(my_start, cdsEnd, 
								"CDS", lastFeature.getID()+"#C",  lastFeature.getStrand(), lastFeature.getID()));
						my_start = cdsEnd + 1;
					}					
					if (my_start <= my_end)
						anno.add(new JapsaFeature(my_start, my_end, 
							codingGene?"5UTR":"UTR", lastFeature.getID()+ (codingGene?"#5":"#U"), lastFeature.getStrand(), lastFeature.getID()));
				}//for
				if (upStream > 0){			
						anno.add(new JapsaFeature(end + 1, end + upStream, 
								"upstream", lastFeature.getID()+"#U",  lastFeature.getStrand(), lastFeature.getID()));
				}	
				///////////////////////////////////////////////////////////////////////////////////////////
			}//else
		}//while	
		
		System.err.println("Ignore " + noIgnore + " genes ");
		
		return annos;
	}
}
