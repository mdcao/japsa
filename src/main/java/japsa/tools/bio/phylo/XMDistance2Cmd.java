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
 * 08/10/2012 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsa.tools.bio.phylo;


import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import japsa.seq.Alphabet.DNA;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.tools.xm.ExpertModelCmd;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;
import japsa.xm.ExpertModel;


@Deployable(scriptName = "jsa.phylo.xmdist2",
scriptDesc = "Generate a distances bet matrix from genomes (potentially not alignable")
public class XMDistance2Cmd  extends CommandLine{	
	public XMDistance2Cmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		//addStdInputFile();		
		addString("output", "output", "Name of the file for output (distances in phylip format)");
		addInt("index", 0, "Index");

		addInt("hashSize", 11, "Hash size");
		addInt("context", 15, "Length of the context");
		addInt("limit", 50, "Expert Limit");
		//addInt("thread", 1, "Number of threads");
		addDouble("threshold", 0.15, "Listen threshold");
		addInt("chance", 20, "Chances");
		addBoolean("binaryHash", false, "Use binary hash or not");
		addString("offsetType", "counts",
			"Way of update offset/palindrome expert: possible value count, subs");
		addBoolean("optimise", false,
			"Running in optimise mode, just report the entropy,recommended for long sequence");
		addInt("checkPoint", 10000000, "Frequency of check point");
		addString("hashType", "hash",
			"Type of Hash table: hash=hashtable, sft=SuffixTree,sfa = SuffixArray");
		addBoolean("selfRep", true,
			"Propose experts from the sequence to compressed?");	

		addStdHelp();		
	} 
	//public static boolean adapt = false;


	public static void main(String[] args) throws Exception {		 		
		CommandLine cmdLine = new XMDistance2Cmd();
		args = cmdLine.stdParseLine(args);
		//ExpertModel eModel = ExpertModelCmd.getExpertModel(cmdLine);


		BufferedReader bf = new BufferedReader(new FileReader("list"));	
		ArrayList<String> list = new ArrayList<String>();
		String str;

		while ( (str = bf.readLine())!=null){
			list.add(str.trim());
		}
		bf.close();			
		int index = cmdLine.getIntVal("index");

		if(index >= list.size()){
			Logging.exit("Wrong index (< " + list.size()+") : " + index, 1);
		}		

		FastaReader sin = new FastaReader("data/" + list.get(index));
		Sequence mySeq = sin.nextSequence(DNA.DNA4());		
		/**************************************************/
		sin.close();

		Sequence [] mS = new Sequence[1];
		mS[0] = mySeq;		
		ExpertModel eModel = ExpertModelCmd.getExpertModel(cmdLine);
		double myE = eModel.encode_optimise(mS);

		SequenceOutputStream outStr = 
			SequenceOutputStream.makeOutputStream(cmdLine.getStringVal("output"));

		outStr.print(mySeq.getName() + "\t" + mySeq.length() + "\t" +myE+"\n");

		mS = new Sequence[2];		
		for (int x = 0;x < index;x++){
			sin = new FastaReader("data/" + list.get(x));			
			Sequence mateSeq = 	sin.nextSequence(DNA.DNA4());
			sin.close();

			mS[0] = mySeq;
			mS[1] = mateSeq;
			
			eModel = ExpertModelCmd.getExpertModel(cmdLine);
			double e_ji = eModel.encode_optimise(mS);
			
			mS[1] = mySeq;
			mS[0] = mateSeq;			
			eModel = ExpertModelCmd.getExpertModel(cmdLine);
			double e_ij = eModel.encode_optimise(mS);			
			outStr.print(mateSeq.getName() + "\t" + mateSeq.length() + "\t" +e_ji+"\t" + e_ij + "\n");
			outStr.flush();
		}		
		outStr.close();
	}	
}
