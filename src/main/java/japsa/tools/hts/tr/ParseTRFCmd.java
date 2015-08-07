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
 * 18/03/2012 - Minh Duc Cao: Revised                                        
 * 10/5/2012 - Minh Duc Cao
 *   -Standardise commandline
 *   -Inlude feature for simulation 
 ****************************************************************************/

package japsa.tools.hts.tr;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Iterator;



/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.trv.parseTRF",
            scriptDesc = "Parse trf output to jsa, bed or tr format")
public class ParseTRFCmd extends CommandLine{	
	public ParseTRFCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", null, "Name of the input file (output of TRF), - for standard input", true);
		addString("output", "-", "Name of output file, - for standard output");
		addString("format", "jsa", "Format of the output file. Options: jsa, bed, xaf");
		addString("sequence", null, "Name of the sequence file (has to be the same file to run TRF)");
		addInt("max", 6,"Maximum unit size");
		addInt("min", 2,"Minimum unit size");	
		
		addStdHelp();		
	} 
	
	/**
	 * Parse the result from Tandem Repeat Finder into various annotation format
	 * Sample : 50501 50569 3 23.0 3 100 0 138 33 33 33 0 1.58 AGC
	 * AGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGC
	 * <Start> <End> <Period> <noRun>....
	 * 
	 *  
	 * @param args
	 */
	public static void main(String[] args) throws Exception {		 		
		CommandLine cmdLine = new ParseTRFCmd();
		args = cmdLine.stdParseLine(args);
		
		String inputFile = cmdLine.getStringVal("input");
		String outputFile = cmdLine.getStringVal("output");				
		String format = cmdLine.getStringVal("format");
		
		int min = cmdLine.getIntVal("min");
		int max = cmdLine.getIntVal("max");
		
		int skip = 0, retain = 0, filter=0;
			
		//Read in TRF
		BufferedReader bf = SequenceReader.openFile(inputFile);		
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(outputFile);
		
		
		SequenceReader seqReader = null;		
		if (cmdLine.getStringVal("sequence") != null){
			seqReader = SequenceReader.getReader(cmdLine.getStringVal("sequence"));
		}
		Sequence seq = null;

		String line = "";				
		JapsaAnnotation anno = null;
		
		//JapsaFeature lastFeature = null;
		//double lastScore = 0;
	//	int lastFactor = 1;
		TandemRepeat lastTR = null;
		
		String id = "";
		while ((line = bf.readLine()) != null) {			
			line = line.trim();
			
			String[] tokens = line.split(" +");
			if(tokens.length <= 1)
				continue;
			if (tokens[0].equals("Sequence:")){				
				if (anno != null){
					write(seq,anno,out,format);
				}
				id = tokens[1];
				
				anno = new JapsaAnnotation(null,tokens[1]);
				anno.addDescription("Short Tandem Repeat obtained from trf for " + id);
				lastTR = null;
				if (seqReader != null){
					seq = seqReader.nextSequence(null);
				}
				continue;
			}			
			
			if (tokens[0].equals("Parameters:")){
				if (anno == null){
					anno = new JapsaAnnotation(null,"");
					anno.addDescription("Short Tandem Repeat obtained from trf");
					lastTR = null;
				}
				
				anno.addDescription(line);
			}
//			count++;
			if (tokens.length >=15){
				int unitSize = tokens[13].length();
				if (unitSize > max || unitSize < min){
					filter ++;
					continue;
				}
				
				//each corresponds to a TR
				if (anno == null){
					//no header previously
					anno = new JapsaAnnotation(null,"");
					anno.addDescription("Short Tandem Repeat obtained from trf");
					lastTR = null;
				}//if anno				
				int start = Integer.parseInt(tokens[0]);
				int end = Integer.parseInt(tokens[1]);
				
				int factor = Integer.parseInt(tokens[4]);
//				int rSize = tokens[13].length();				
				double unitNo = Double.parseDouble(tokens[3]);
//				566 601 4 9.0 4 81 0 45 27 25 22 25 2.00 TGAC TGACTGCCTGAATAACTGACTGACTGACTGACTGAC								
				
				TandemRepeat tr = 
						new TandemRepeat(id, start, end, factor, unitNo);
				
				tr.setScore(Double.parseDouble(tokens[7]));
				tr.setUnitNo(unitNo);				
				tr.setUnit(tokens[13]);
				
				
				tr.setID("P"+start);
				tr.addDesc("@F:"+line);
				
				tr.addDesc("@R:" + factor);
				
				//Unit no
				tr.addDesc("@N:" + tokens[3]);				
				// percent Match				
				tr.addDesc("@M:" + tokens[5]);
				// percent Indel
				tr.addDesc("@I:" + tokens[6]);			
				// score
				tr.addDesc("@S:" + tokens[7]);
				// Entropy
				tr.addDesc("@E:" + tokens[12]);
				// Unit
				tr.addDesc("@U:" + tokens[13].substring(0, factor));
				//data sequence
				tr.addDesc("@D:" + tokens[14]);
				
				
				
				//check if overlaping with the last 
				if (lastTR != null && lastTR.getStart() < tr.getEnd() && tr.getStart() < lastTR.getEnd()){
					if (tr.getScore() > lastTR.getScore()){
						anno.remove(lastTR);
						skip ++;retain --;
						Logging.warn("Skip [" + lastTR.getStart() + " " + lastTR.getEnd() + "](" +lastTR.getScore()+") because of [" + tr.getStart() + " " + tr.getEnd() + "](" +tr.getScore() +")");
					}else{
						skip ++;
						Logging.warn("Skip [" + tr.getStart() + " " + tr.getEnd() + "](" +tr.getScore()+") because of [" + lastTR.getStart() + " " + lastTR.getEnd() + "](" +lastTR.getScore() +")");
						continue;
					}
				}
				
				anno.add(tr);
				retain ++;
				lastTR = tr;
			}//if
		}// while
		
		if (anno != null)
			write(seq, anno, out, format);
			
		out.close();			
		Logging.info(" Retain " + retain + " records, skip " + skip + " because of overlapping and filter "+ filter + " of irrelevant period size");
	}
	
	/**
	 * Write the annotation to a stream in the format specified by format 
	 * parameter
	 * 
	 * @param anno
	 * @param out
	 * @param format
	 * @throws IOException
	 */
	
	private static void write(Sequence seq, JapsaAnnotation anno, SequenceOutputStream out, String format)
	throws IOException{		 
		if ("bed".equals(format))
			anno.writeBED(out);
		else if (format.startsWith("xaf")){			
			String[] headers = TandemRepeat.STANDARD_HEADER;			
			
			String desc = anno.getDescription();			
			String [] toks = desc.split("\n");
			for (int i = 0;i < toks.length;i++)
				if (toks[i].length() > 0)
					out.print("##" + toks[i]+"\n");
				
			//write out the header
			out.print("#H:" + headers[0]);
			for (int i = 1; i < headers.length; i++)
				 out.print("\t" + headers[i]);
			 
			out.print('\n');
			Iterator<JapsaFeature> trIter =  anno.iterator();
			
			while (trIter.hasNext()){
				TandemRepeat tr = (TandemRepeat) trIter.next();
				out.print(tr.toString(headers, false) + "\n");
			}			
			out.print('\n');
			
		}else{//default jsa format
			anno.sortFeatures();
			JapsaAnnotation.write(seq, anno, out);
			//anno.write(out);
		}
	}
}



