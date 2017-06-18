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
 * 11/05/2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.tools.bio.np;

import java.io.IOException;

import japsa.seq.Alphabet;
import japsa.seq.FastqReader;
import japsa.seq.FastqSequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.nanopore.NanoporeReaderStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * 
 * @author minhduc
 *
 */
@Deployable(
		scriptName = "jsa.np.filter", 
		scriptDesc = "Filter nanopore reads data from fastq file",
		seeAlso = "jsa.np.npreader, jsa.util.streamServer, jsa.util.streamClient")
public class NanoporeReadFilterCmd extends CommandLine{	
	public NanoporeReadFilterCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		CommandLine.Option inputOpt =
				addString("input", null, "Name of the input file, - for standard input", true);

		CommandLine.Option outputOpt =
				addString("output", null, "Name of the output file, - for standard output", true);

		CommandLine.Option lenMinOpt =
				addInt("lenMin", 0, "Minimum sequence length");

		//CommandLine.Option lenMaxOpt =
		addInt("lenMax", Integer.MAX_VALUE, "Minimum sequence length");

		CommandLine.Option qualMinOpt =
				addDouble("qualMin", 0.0, "Minimum average quality");

		//CommandLine.Option qualMaxOpt =
		addDouble("qualMax", 1000.0, "Maximum average quality");

		CommandLine.Option groupOpt =
				addString("group", "", "Group need to be extracted, leave blank for selecting all groups");

		CommandLine.Option excl2DOpt =
				addBoolean("excl2D", false, "Exclude 2D reads");

		CommandLine.Option exclTempOpt =
				addBoolean("exclTemp", false, "Exclude template reads");

		CommandLine.Option exclCompOpt =
				addBoolean("exclComp", false, "Exclude complement reads");

		//CommandLine.Option formatOpt =
		addString("format", "fastq", "Format of the output file");

		addStdHelp();


		inputOpt.setGalaxySetting(new GalaxySetting("data", "fastqsanger",false));
		groupOpt.setGalaxySetting(new GalaxySetting("text", null, false));

		lenMinOpt.setGalaxySetting(new GalaxySetting("integer", null,false));
		//lenMaxOpt.setGalaxySetting(new GalaxySetting("integer", null,false));
		qualMinOpt.setGalaxySetting(new GalaxySetting("double", null,false));
		//qualMaxOpt.setGalaxySetting(new GalaxySetting("double", null,false));

		excl2DOpt.setGalaxySetting(new GalaxySetting("boolean", null,false));
		exclTempOpt.setGalaxySetting(new GalaxySetting("boolean", null,false));
		exclCompOpt.setGalaxySetting(new GalaxySetting("boolean", null,false));


		GalaxySetting outputGalaxy = new GalaxySetting("data", "fastqsanger",true);
		//outputGalaxy.setLabel
		outputOpt.setGalaxySetting(outputGalaxy);


		setGalaxy(annotation.scriptName());
	} 

	public static void main(String[] args) throws IOException {
		CommandLine cmdLine = new NanoporeReadFilterCmd();
		args = cmdLine.stdParseLine(args);

		String output = cmdLine.getStringVal("output");
		String input = cmdLine.getStringVal("input");
		int lenMin  = cmdLine.getIntVal("lenMin");
		int lenMax  = cmdLine.getIntVal("lenMax");
		double qualMin  = cmdLine.getDoubleVal("qualMin");
		double qualMax  = cmdLine.getDoubleVal("qualMax");

		boolean exclude2D =  cmdLine.getBooleanVal("excl2D");
		boolean excludeTemplate =  cmdLine.getBooleanVal("exclTemp");
		boolean excludeComplement =  cmdLine.getBooleanVal("exclComp");

		String group = cmdLine.getStringVal("group").trim();

		String format = cmdLine.getStringVal("format");

		boolean fastaOutput = "fasta".equals(format.trim().toLowerCase());

		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);

		FastqReader reader = "-".equals(input)? (new FastqReader(System.in) ) 
				: (new FastqReader(input));
		FastqSequence seq;



		while ( (seq = reader.nextSequence(Alphabet.DNA())) != null){
			//Min length
			if (seq.length() < lenMin)
				continue;

			//max length
			if (seq.length() > lenMax)
				continue;

			double qual = -1;

			//min quality
			if (qualMin > 0.0){
				qual = NanoporeReaderStream.averageQuality(seq);
				if (qual < qualMin)
					continue;				
			}

			//max quality			
			if (qualMax < 1000.0){
				if (qual < 0)
					qual = NanoporeReaderStream.averageQuality(seq);
				if (qual >= qualMax)
					continue;		
			}

			if (excludeComplement && seq.getName().contains("complement"))
				continue;

			if (excludeTemplate && seq.getName().contains("template"))
				continue;

			if (exclude2D && seq.getName().contains("twodim"))
				continue;

			if (group.length() > 0){
				String [] toks = seq.getName().split(" ");
				boolean match = false;
				for (String tok:toks){
					if (tok.startsWith("group=")&&tok.substring(6).equals(group)){
						match = true;
						break;
					}
				}
				if (!match)
					continue;//while
			}

			//done all the fitlering
			if (fastaOutput)
				seq.writeFasta(sos);
			else
				seq.print(sos);

		}

		reader.close();		
		sos.close();
	}//main
}

/*RST*
---------------------------------------
*jsa.np.filter*: Filter sequencing data 
---------------------------------------

*jsa.np.filter* filters sequencing data based on sequence read type, length and
quality. Examples of its usage can be found on jsa.np.npreader_.

<usage>


 *RST*/
