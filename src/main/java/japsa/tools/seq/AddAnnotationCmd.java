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
 * 23/01/2014 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsa.tools.seq;


import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.JapsaFileFormat;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;



/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.seq.addanno",
            scriptDesc = "Add annotations to a Japsa file")
public class AddAnnotationCmd extends CommandLine{	
	public AddAnnotationCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", null, "Name of the input file in jsa format, - for standard input", true);
		addString("anno", null, "Name of the new annotation file in jsa format",true);		
		addString("output", "-", "Name of the output file,  - for standard output");		
		addInt("upStr", 0 , "size of upstream regions for gene, <=0 for not adding");
		addInt("downStr", 0 , "size of downstream regions for gene, <=0 for not adding");
		
		addStdHelp();		
	} 
	/**
	 * 
	 * @param args
	 */

	public static void main(String[] args) throws Exception {

		CommandLine cmdLine = new AddAnnotationCmd();		
		args = cmdLine.stdParseLine(args);
		
		/**********************************************************************/		

		String output =   cmdLine.getStringVal("output");
		String annoFile = cmdLine.getStringVal("anno");
		String input =    cmdLine.getStringVal("input");
		int up =          cmdLine.getIntVal("upStr") ;
		int down =          cmdLine.getIntVal("downStr") ;

		//Get dna 		
		
		JapsaFileFormat jsaReader = new  JapsaFileFormat(input);
		JapsaFileFormat annoReader = new  JapsaFileFormat(annoFile);
		
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream((output));			

		JapsaAnnotation mainAnno;
		
		while ((mainAnno = jsaReader.readAnnotation()) != null){
			JapsaAnnotation anno = annoReader.readAnnotation();
			if (anno != null){
				//Make sure the IDs are the same
				if (mainAnno.getAnnotationID().equals(anno.getAnnotationID())){
					String annoDesc = anno.getDescription();
					if (annoDesc.length() > 0)
						mainAnno.addDescription(annoDesc);
					
					
					for (JapsaFeature feature:anno.getFeatureList()){
						mainAnno.add(feature);
						//Add upstream for gene
						if (feature.getType().equals("gene") && up > 0){
							JapsaFeature upFeature;						
							if (feature.getStrand() == '-')//complement
								upFeature = new JapsaFeature(feature.getEnd() + 1, feature.getEnd() + up, "UStr", feature.getID()+"UStr",'-',feature.getID());
							else	
								upFeature = new JapsaFeature(feature.getStart() - up, feature.getStart() - 1, "UStr", feature.getID()+"UStr",'+',feature.getID());
							
							//Fix upStream
							if (upFeature.getStart() < 1) upFeature.setStart(1);													
							Sequence seq =  mainAnno.getSequence();
							if (seq != null &&  upFeature.getEnd() > seq.length()) upFeature.setEnd(seq.length());
							
							if (upFeature.getLength() > 1)
								mainAnno.add(upFeature);						
						}
						
						//Add downstream for gene
						if (feature.getType().equals("gene") && down > 0){
							JapsaFeature downFeature;						
							if (feature.getStrand() == '-')//complement
								downFeature = new JapsaFeature(feature.getStart() - up, feature.getStart() - 1, "DStr", feature.getID()+"DStr",'-',feature.getID());
							else	
								downFeature = new JapsaFeature(feature.getEnd() + 1, feature.getEnd() + up, "DStr", feature.getID()+"DStr",'+',feature.getID());
							
							//Fix upStream
							if (downFeature.getStart() < 1) downFeature.setStart(1);													
							Sequence seq =  mainAnno.getSequence();
							if (seq != null &&  downFeature.getEnd() > seq.length()) downFeature.setEnd(seq.length());
							
							if (downFeature.getLength() > 1)
								mainAnno.add(downFeature);				
						}						
						
					}
					//mainAnno.sortFeatures();
					mainAnno.write(out);
				}else{
					Logging.warn("The IDs are not identical, annotation not added!");
				}
				
			}else//if anno = null
				Logging.warn("Annotations for " + mainAnno.getAnnotationID() + " not found");
		}

		out.close();
		jsaReader.close();
		annoReader.close();		
	}
}
