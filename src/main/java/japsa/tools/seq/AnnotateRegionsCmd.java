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
 * 09/05/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.tools.seq;

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.JapsaFileFormat;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

import java.io.IOException;


/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
@Deployable(scriptName = "jsa.seq.annotate",
scriptDesc = "Annotate a list of regions using some annotation such as RefSeq")
public class AnnotateRegionsCmd extends CommandLine{	
	public AnnotateRegionsCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", null, "Name of the input file, - for standard input", true);
		addString("anno", "anno", "Name of the annotation file,  - for standard output");
		addString("output", "output", "Output");
			
		addStdHelp();		
	} 
	
	public static void main(String[] args) throws IOException {		
		CommandLine cmdLine = new AnnotateRegionsCmd();		
		args = cmdLine.stdParseLine(args);
		
		/**********************************************************************/
		JapsaFileFormat annoF = new JapsaFileFormat(cmdLine.getStringVal("anno"));
		JapsaFileFormat inputF= new JapsaFileFormat(cmdLine.getStringVal("input"));
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("output")); 
		
		JapsaAnnotation inputAnno = null;
		JapsaAnnotation annoAnno = annoF.readAnnotation();
		while (( inputAnno = inputF.readAnnotation())!=null){
			//Move the the next annotations if not match
			Logging.info(inputAnno.getAnnotationID());
			while (!annoAnno.getAnnotationID().equals(inputAnno.getAnnotationID())){
				annoAnno = annoF.readAnnotation();
			}
			System.out.println(inputAnno.getDescription() + 
					"      " + annoAnno.getDescription());			
			
			//TODO: A slow version
			for (int i = 0; i < inputAnno.numFeatures();i++){
				JapsaFeature inputFeature = inputAnno.getFeature(i);
				
				String dDest = "";
				for (int j = 0; j < annoAnno.numFeatures(); j++){
					
					JapsaFeature annoFeature = annoAnno.getFeature(j);
					if (inputFeature.getStart() < annoFeature.getEnd()
						&& 	inputFeature.getEnd() > annoFeature.getStart()
					){
						dDest = dDest +("@"+annoFeature.getType() + "(" + annoFeature.getID()+"," + annoFeature.getParent()+")");
					}//if
				}//for
				if(dDest.length() > 0)
					inputFeature.addDesc("@"+dDest);
			}
			inputAnno.writeAnnotation(out);
		}
		out.close();
		annoF.close();
		inputF.close();
				
	}
	
	
	
	
	
	
}
