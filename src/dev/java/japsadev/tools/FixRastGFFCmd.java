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
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
	scriptName = "jsa.dev.fixRastGff",
	scriptDesc = "Fix rast gff file"
	)
public class FixRastGFFCmd extends CommandLine{	
	public FixRastGFFCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("input", null, "Name of gff input file", true);
		addString("sequence", null, "Name of fasta sequence file", true);
		addString("output", "-", "Name of output file");
		
		
		//addBoolean("reverse",false,"Reverse sort order");
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new FixRastGFFCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		
		String input = cmdLine.getStringVal("input");
		String sequence = cmdLine.getStringVal("sequence");
		String output = cmdLine.getStringVal("output");
				
		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);
		sos.print("##gff-version 3\n");
				
		SequenceReader reader = SequenceReader.getReader(sequence);
		
		Sequence seq;
		
		while ((seq = reader.nextSequence(Alphabet.DNA5()))!= null){
			sos.print("##sequence-region " + seq.getName() + " 1 " + seq.length() + '\n');	
		}
		reader.close();
		
		BufferedReader bf = new BufferedReader(new FileReader(input));
		String line;
		
		ArrayList<Record> records = new ArrayList<Record>(); 
		while ((line = bf.readLine()) != null){
			if (line.startsWith("#"))
				continue;
			line = line.trim();
			String [] toks = line.split("\t");
			Record rec = new Record();
			rec.chr = toks[0]; 
			rec.start = Integer.parseInt(toks[3]);
			rec.end = Integer.parseInt(toks[4]);
			rec.str = line;
			records.add(rec);
		}
		bf.close();
		
		
		Collections.sort(records);
		
		for (Record rec:records){
			sos.print(rec.str);
			sos.println();
		}		
		sos.close();		
	}
	
	static class Record implements Comparable<Record>{
		String chr = "chr";
		int start, end;
		String str;

		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(Record o) {
			int comp = chr.compareTo(o.chr);
			
			if (comp == 0)
				comp = Integer.compare(start, o.start);
			
			if (comp == 0)
				comp = Integer.compare(end, o.end);
				
			return comp;
		}
		
	}
	
}
/*RST*



 
  
  
  
  
  
  
  
  
*RST*/
  
