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
 * File: Japsa2TR.java
 * 16/11/2013 - Minh Duc Cao: Created
 *
 ****************************************************************************/

package japsa.tools.hts.tr;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.JapsaFileFormat;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * FIXME: Need to test
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 */
@Deployable(scriptName = "jsa.trv.jsa2tr",
            scriptDesc = "Convert tandem repeat annotation in JAPSA format to TR format")
public class Japsa2TRCmd extends CommandLine{	
	public Japsa2TRCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", null, "Name of input JSA file, - for standard input",true);
		addString("output", "-", "Name of output file, - for standard out");			
		
		addStdHelp();		
	} 

	/**
	 * @param args
	 */		
	public static void main(String[] args) throws IOException{		
		CommandLine cmdLine = new Japsa2TRCmd();		
		args = cmdLine.stdParseLine(args);
		
		
		String inputFile = cmdLine.getStringVal("input");		
		String outputFile = cmdLine.getStringVal("output");
		
		SequenceOutputStream ps = SequenceOutputStream.makeOutputStream(outputFile);

		//Read tnr
		System.out.println(" Reading tandem repeats from " + inputFile);
		JapsaAnnotation trAnno = null;//JapsaAnnotation.readDataFromFile(inputFile);		
		 
		JapsaFileFormat reader = null;
		
		if (inputFile.equals("-"))
			reader = new JapsaFileFormat(System.in);
		else
			reader = new JapsaFileFormat(inputFile);
		
		
		while ( (trAnno = reader.readAnnotation()) != null){
			ArrayList<TandemRepeat> trs = new ArrayList<TandemRepeat>();

			Iterator<JapsaFeature> iter = trAnno.iterator();
			while (iter.hasNext()) {
				trs.add(new TandemRepeat(iter.next()));
			}
			
			System.out.println(" Writing " + trs.size() + " to " + outputFile);

			
			String annoDesc = trAnno.getDescription();
			ArrayList<String> dList = new ArrayList<String>();
			String[] toks = annoDesc.split("\n");

			for (int i = 0; i < toks.length; i++)
				dList.add("#" + toks[i]);

			dList.add("#TNR list from parsed");
			
			//ps.println("#ID:" + tnrs.getID());
			//writeToFile(ps, STANDARD_HEADER, trs, dList, true);
			TandemRepeat.writeToFile(ps, TandemRepeat.FULL_HEADER_WITH_ANNOTATIONS, trs, dList, true);			
		}
				
		ps.close();
		reader.close();
	}
}
