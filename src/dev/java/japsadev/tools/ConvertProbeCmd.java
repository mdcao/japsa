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
 * 14/09/2016 - Minh Duc Cao: Created                                        
 ****************************************************************************/

package japsadev.tools;


import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.IntArray;
import japsa.util.deploy.Deployable;


/**
 * @author minhduc
 *
 */
@Deployable(
		scriptName = "jsa.dev.convertProbe", 
		scriptDesc = "Convert probes to sam and bed file"
		)
public class ConvertProbeCmd extends CommandLine{
	//CommandLine cmdLine;
	public ConvertProbeCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		//Input/output
		addString("probe", "-", "File containing probe name");
		addString("reference", null, "Reference");		
		addString("output", "output", "Output prefix");
		
		addInt("good",100,"good threshold");
		addInt("bad",100,"bad threshold");

		addStdHelp();
	}
	public static void main(String [] args) throws IOException{
		CommandLine cmdLine = new ConvertProbeCmd ();
		args = cmdLine.stdParseLine(args);

		/**********************************************************************/		
		String probe       =  cmdLine.getStringVal("probe");
		String output       =  cmdLine.getStringVal("output");
		String reference     =  cmdLine.getStringVal("reference");
		
		int goodThreshold = cmdLine.getIntVal("good");
		int badThreshold = cmdLine.getIntVal("bad");

		IntArray chrLens = new IntArray(24);
		ArrayList<String> chrList = new ArrayList<String> (24);
		HashMap<String,Chromosome> chromsMap = new HashMap<String,Chromosome> (24);

		FastaReader.Faster fread = new FastaReader.Faster(reference);
		Sequence seq;

		// Timer timer = new Timer();
		int index = 0;
		while ((seq = fread.nextSequence(Alphabet.DNA16())) != null) {
			String id = seq.getName().split("\\s")[0];
			chrList.add(id);
			chrLens.add(seq.length());

			Chromosome chrom = new Chromosome(index, seq.length());
			chromsMap.put(id, chrom);			
			index ++;
		}
		fread.close();		
		BitSet [] bitSets = new BitSet[chrList.size()];
		for (int i = 0; i < bitSets.length;i++)
			bitSets[i] = new BitSet(chrLens.get(i));


		SequenceOutputStream samOut = SequenceOutputStream.makeOutputStream(output + ".sam");
		SequenceOutputStream badSamOut = SequenceOutputStream.makeOutputStream(output + "_bad.sam");
		for (int i = 0; i< chrLens.size();i++){
			samOut.print("@SQ\tSN:"+chrList.get(i)+"\tLN:"+chrLens.get(i) + '\n');
			badSamOut.print("@SQ\tSN:"+chrList.get(i)+"\tLN:"+chrLens.get(i) + '\n');
		}		
			
		

		SequenceOutputStream badOut = SequenceOutputStream.makeOutputStream(output + ".bad");
		BufferedReader reader = SequenceReader.openFile(probe);

		String line = "";
		
		while ( (line = reader.readLine()) !=null){
			line = line.trim();
			String [] toks = line.split("\t");
			
			int badCount = Integer.parseInt(toks[2]);
			int goodCount = Integer.parseInt(toks[1]);
				
			String name = toks[0];
			
			String readSeq = toks[3];
			String readQual = toks[4];

			toks = name.split("_");			
			int end = Integer.parseInt(toks[2]);
			int start = Integer.parseInt(toks[1]);
			String chr = toks[0];			


			if (goodCount > goodThreshold){
				badOut.print(name + "\t" + goodCount + "\t" + badCount + "\t(good>" + goodThreshold +")\n");
				badSamOut.print(name+"\t0\t"+chr+"\t"+start + "\t60\t" + (end-start+1) +"M\t*\t0\t0\t"+readSeq + "\t" + readQual + "\n");
				continue;
			}
			
			if (badCount > badThreshold){				
				badOut.print(name + "\t" + goodCount + "\t" + badCount + "\t(bad>" + badThreshold +")\n");
				badSamOut.print(name+"\t0\t"+chr+"\t"+start + "\t60\t" + (end-start+1) +"M\t*\t0\t0\t"+readSeq + "\t" + readQual + "\n");
				continue;
			}			
			
			samOut.print(name+"\t0\t"+chr+"\t"+start + "\t60\t" + (end-start+1) +"M\t*\t0\t0\t"+readSeq + "\t" + readQual + "\n");

			BitSet bitSet = bitSets[chromsMap.get(chr).index];
			for (int i = start - 1; i< end;i++){
				bitSet.set(i);
			}			
		}		
		reader.close();
		samOut.close();
		badOut.close();

		SequenceOutputStream bedOut = SequenceOutputStream.makeOutputStream(output + ".bed");		


		for (int x=0; x < chrList.size();x++){
			String chrName = chrList.get(x);
			BitSet myBitSet = bitSets[x];
			int regionStart = -1;				
			for (int i = 0; i < chrLens.get(x);i++){
				if (myBitSet.get(i) && (regionStart < 0)){
					//start of a new region
					regionStart = i;						
				}else if (!myBitSet.get(i) && (regionStart >= 0)){
					//end of a region										
					bedOut.print(chrName + "\t" + regionStart + "\t" + i + "\n");

					regionStart = -1;
				}
			}
			if (regionStart >=0){
				bedOut.print(chrName + "\t" + regionStart + "\t" + chrLens.get(x) + "\n");				
			}
		}//for
		bedOut.close();

	}

	/**
	 * Implement regions that may be capturable
	 * @author minhduc
	 *
	 */
	static class Chromosome{
		int index;
		int length;
		public Chromosome(int i, int l){
			index = i;
			length = l;
		}
	}

}
