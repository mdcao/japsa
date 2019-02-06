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
import japsa.seq.FastqSequence;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;


/**
 * @author Son Nguyen
 * 
 */
@Deployable(scriptName = "jsa.seq.binseq",
scriptDesc = "Bin the FASTQ batch file into species mapped reads based on species2map from jsa.np.rtSpeciesTyping")
public class BinningSequenceCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(BinningSequenceCmd.class);

	public BinningSequenceCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("sequence", null, "The FASTQ file of all nanopore reads",true);
		addString("map", null, "species2read map file from running species typing",true);		
		addString("outputDir", "bins", "Output directory containing binned reads.");		

		addString("filter", "", "List of species (separated by comma) to excluded from the binning");
		addInt("minRead", 0, "Mininum number of read count for the output bins");

		addStdHelp();		
	} 


	public static void main(String[] args) throws IOException {
		CommandLine cmdLine = new BinningSequenceCmd();	
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/		
		String seqFile = cmdLine.getStringVal("sequence");		
		String mapFile = cmdLine.getStringVal("map");
		String filter = cmdLine.getStringVal("filter");
		String outputDir = cmdLine.getStringVal("output");
		int minRead = cmdLine.getIntVal("minRead");	

		/**********************************************************************/
		binning(seqFile, mapFile, filter, outputDir, minRead);

	}

	public static void binning(String seqFile, String mapFile, String filter, String outputDir, int minRead) throws IOException{

		//0.Make the filter conditions
		HashSet<String> filterSet = new HashSet<String>();
		if(!filter.isEmpty()){
			String[] toks = filter.split(",");
			for(String tok:toks)
				filterSet.add(tok);
		}
		
		//1.Load the map
		HashMap<String, String> binMap = new HashMap<String, String>(); //read to bin name
		BufferedReader mapReader = new BufferedReader(new FileReader(mapFile));	
		String line=mapReader.readLine();
		String curSpecies=null;
		ArrayList<String> curList=null;
		while(line!=null){
			if(line.startsWith(">")){
				if(curSpecies!=null && curList!=null){
					if(!filterSet.contains(curSpecies) && curList.size()>minRead){
						for(String readID:curList){
							if(binMap.containsKey(readID)){
								LOG.info("Read {} map to multiple species -> put to unclassifed!", readID);
								binMap.put(readID, "unclassifed");
							}else
								binMap.put(readID, curSpecies);
						}
							
					}else
						LOG.info("{} excluded, readCount={}",curSpecies,curList.size());
				}
				curSpecies=line.substring(1).trim();
				curList=new ArrayList<>();
			}else
				curList.add(line.trim());
			
			line=mapReader.readLine();
		}
		mapReader.close();		
		
		//2.Binning based on the map
		SequenceReader seqReader = SequenceReader.getReader(seqFile);
		Sequence seq = seqReader.nextSequence(Alphabet.DNA());
		String extension;
		if(seq instanceof FastqSequence)
			extension=".fastq";
		else
			extension=".fasta";
		
		
		//set appropriate output stream first
		Set<String> outBins=new HashSet<String>(binMap.values());
		HashMap<String, SequenceOutputStream> bin2File = new HashMap<String, SequenceOutputStream>();
		File 	file=new File(seqFile),
				outDir=new File(outputDir);		
		if(!outDir.exists() && outDir.mkdirs()){
			LOG.info("Output binned reads into {}", outputDir);
		}
		else{
			LOG.error("Output folder {} already exist or cannot be created!\n\tERROR: Please specify another output directory...",outputDir);
			System.exit(1);
		}
		
		for(String bin:outBins){
			String outFile=outputDir+File.separator+bin+extension;
			bin2File.put(bin, new SequenceOutputStream(new FileOutputStream(outFile)));
		}
		
		
		String binName=binMap.get(seq.getName());
		if(binName!=null){
			SequenceOutputStream out=bin2File.get(binName);
			if(out!=null)
				seq.print(out);
		}
		seqReader.close();
	}
}


/*RST*
-------------------------------------------
*jsa.seq.binseq*: Bin the FASTQ sequences based on species2read 
-------------------------------------------

*jsa.seq.binseq* output the binned FASTQ reads that map to every detected species by jsa.np.rtSpeciesTyping

<usage> 

*RST*/
