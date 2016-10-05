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

package japsadev.tools;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.combineProdigal",
	scriptDesc = "Sample script description"
	)
public class CombineProdigalCmd extends CommandLine{	
	public CombineProdigalCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("input","-","Input from prodigal");
		addString("gff",null,"Extra cds",true);

		//addBoolean("reverse",false,"Reverse sort order");	

		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new CombineProdigalCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String input = cmdLine.getStringVal("input");
		String gff = cmdLine.getStringVal("gff");

		FileInputStream aReader = new FileInputStream(gff);		
		ArrayList<JapsaAnnotation> annos = JapsaAnnotation.readMGFF(aReader,0,0,"CDS");
		aReader.close();

		JapsaAnnotation anno = null;
		JapsaFeature feature = null;

		int index = -1;
		int cdsIndex = 0;
		int annoIndex = 0;


		BufferedReader bf = SequenceReader.openFile(input);
		String line;

		while ((line = bf.readLine())!= null){
			line = line.trim();
			if (line.startsWith("#")){
				int pos = line.indexOf("seqhdr=");
				if (pos > 0){
					String ctg = line.substring(pos + 8, line.length() -1);
					index ++;
					if (index > annos.size()){
						Logging.exit("index out of range",1);
					}
					anno = annos.get(index);					
					if (!anno.getAnnotationID().equals(ctg)){
						Logging.exit("Unexpected " + ctg + " vs " + anno.getAnnotationID(),1);
					}
					cdsIndex = 0;
					annoIndex = 0;
					feature = (anno.numFeatures()>annoIndex)?anno.getFeature(annoIndex):null;
				}				
				System.out.println(line);
				continue;
			}

			String [] toks = line.split("_");
			int start = Integer.parseInt(toks[1]);
			int end = Integer.parseInt(toks[2]);
			char strain = toks[3].charAt(0);

			//Output any cds from extra annotation that are before of the current cds
			while (feature != null && feature.getEnd() < start){
				cdsIndex ++;
				System.out.println(">"+cdsIndex + "_" + feature.getStart()+"_"+feature.getEnd()+"_"+feature.getStrand());
				annoIndex ++;
				feature = (anno.numFeatures()>annoIndex)?anno.getFeature(annoIndex):null;				
			}

			//assert  feature==null || feature.getEnd() >= start
			while (feature != null && feature.getStart() <= end){
				annoIndex ++;
				feature = (anno.numFeatures()>annoIndex)?anno.getFeature(annoIndex):null;				
			}

			//assert  feature==null || feature.getStart() > end			
			//if (feature == null || end < feature.getStart()){
			cdsIndex ++;
			System.out.println(">"+cdsIndex+"_"+start+"_"+end+"_"+strain);
		}

	}

}
/*RST*












 *RST*/

