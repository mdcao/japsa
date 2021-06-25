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

package japsadev.tools.work;

import japsa.seq.Alphabet;
import japsa.seq.FastqReader;
import japsa.seq.FastqSequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;

/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
		scriptName = "jsa.dev.fixfastq",
		scriptDesc = "Fix fastq files for the plasmaDNA samples"
		)
public class FixFastqNameCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(FixFastqNameCmd.class);
	public FixFastqNameCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " file1 file2");
		setDesc(annotation.scriptDesc());

		addString("output", "output",  "Prefix of the output");
		//addBoolean("reverse",false,"Reverse sort order");
		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {		

		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new FixFastqNameCmd();		
		args = cmdLine.stdParseLine(args);
		if (args.length != 2){
			System.err.println(cmdLine.usage());
			System.exit(1);
		}
		/**********************************************************************/		
		String output = cmdLine.getStringVal("output");

		SequenceOutputStream out1 = SequenceOutputStream.makeOutputStream(output + "_1.fq.gz");
		SequenceOutputStream out2 = SequenceOutputStream.makeOutputStream(output + "_2.fq.gz");
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(output + "_S.fq.gz");

		FastqReader reader1 = new FastqReader(args[0]);
		FastqReader reader2 = new FastqReader(args[1]);

		
		String regex = " |\\.";
		//Pattern pattern = Pattern.compile(" |.");
		FastqSequence seq2 = reader2.nextSequence(Alphabet.DNA());		
		String [] toks = seq2.getName().split(regex);
		int readNo2 = Integer.parseInt(toks[1]);
		seq2.setName(toks[0]+"."+toks[1]);


		FastqSequence seq1 = reader1.nextSequence(Alphabet.DNA());		
		toks = seq1.getName().split(regex);
		int readNo1 = Integer.parseInt(toks[1]);
		seq1.setName(toks[0]+"."+toks[1]);

		while (readNo1 > 0 || readNo2 > 0){
			if (readNo2 == 0 || readNo1 < readNo2){
				seq1.print(out);
				seq1 = reader1.nextSequence(Alphabet.DNA());
				if (seq1 != null){
					toks = seq1.getName().split(regex);
					readNo1 = Integer.parseInt(toks[1]);
					seq1.setName(toks[0]+"."+toks[1]);
				}else
					readNo1 = 0;
				
			}else if (readNo1 == readNo2){				
				seq1.print(out1);
				seq2.print(out2);

				seq1 = reader1.nextSequence(Alphabet.DNA());
				if (seq1 != null){
					toks = seq1.getName().split(regex);
					readNo1 = Integer.parseInt(toks[1]);
					seq1.setName(toks[0]+"."+toks[1]);
				}else
					readNo1 = 0;
				
				seq2 = reader2.nextSequence(Alphabet.DNA());
				if (seq2 != null){
					toks = seq2.getName().split(regex);
					readNo2 = Integer.parseInt(toks[1]);
					seq2.setName(toks[0]+"."+toks[1]);
				}else
					readNo2 = 0;				

			}else{
				LOG.error("Dont understand this " + readNo1 + " vs " + readNo2);
				System.exit(1);
			}
		}//while
		reader1.close();reader2.close();out.close();out1.close();out2.close();
	}
}