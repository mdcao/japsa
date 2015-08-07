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
 * 08/04/2012 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsa.tools.hts.tr;

import japsa.bio.tr.TandemRepeatVariant;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.JapsaMath;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


/**
 * @author minhduc
 * 
 */
@Deployable(scriptName = "jsa.trv.vcf2trv",
scriptDesc = "Convert variation from cvf format to trv format")
public class VCF2TRVCmd extends CommandLine{	
	public VCF2TRVCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("tr", null,	"Name of the tandem repeat file",true);		
		addString("vcf", null,	"Name of the results in vcf format",true);
		addInt("min", 3,"Minimum qual required");
		addString("output", null,	"Output", true);
		addString("sample", null,	"Sample ID");		
		addBoolean("haploid", false,	"Interpret as haploid");

		addStdHelp();		
	} 	
	public static void main(String[] args) throws Exception {		
		CommandLine cmdLine = new VCF2TRVCmd();
		args = cmdLine.stdParseLine(args);


		String  trf = cmdLine.getStringVal("tr");
		String  vcf = cmdLine.getStringVal("vcf");	
		String  out = cmdLine.getStringVal("output");
		String sample = cmdLine.getStringVal("sample");

		if ("-".equals(trf) || "-".equals(vcf)){
			System.err.println("ERROR tr and cvf cannot be -\n" + cmdLine.usageString());
			System.exit(-1);
		}

		int qual = cmdLine.getIntVal("min");

		// getInsertSizeSlow(samFile, output);
		ArrayList<TandemRepeatVariant> trList = TandemRepeatVariant.readFromFile(trf);

		ArrayList<TandemRepeatVariant> samtoolList = new ArrayList<TandemRepeatVariant>(trList.size());
		HashMap<String, Integer> chrIndex = new HashMap<String, Integer>();
		int index = 0;

		for (int i = 0; i < trList.size(); i++){
			TandemRepeatVariant trv = trList.get(i).tandemRepeatClone();

			trv.setMean(0);
			trv.setVar(0);
			trv.setVar2(0);
			trv.setConfidence(0);
			trv.setStd(0);

			if (chrIndex.get(trv.getChr()) == null){
				chrIndex.put(trv.getChr(),index);
				index ++;
			}

			samtoolList.add(trv);
		}
		readVCFResults(samtoolList, chrIndex, vcf, qual);

		if (sample == null)
			sample = "";
		else
			sample = sample + ".";


		SequenceOutputStream outS = SequenceOutputStream.makeOutputStream(out);
		if (cmdLine.getBooleanVal("haploid")){
			outS.print("#H:chr\tID\tstart\tend\tperiod\t"+sample+"var\tconfidence\n");
			for (int i = 0; i < samtoolList.size(); i++){
				TandemRepeatVariant trv = samtoolList.get(i); 
				outS.print(trv.getChr()+"\t"+ trv.getTandemRepeat().getID()+"\t"+trv.getStart()+"\t"+trv.getEnd()+"\t" + trv.getTandemRepeat().getPeriod() + "\t" + (trv.getVar()+trv.getVar2())/2+"\t" + trv.getConfidence() + "\n");
			}
		}else{
			outS.print("#H:chr\tID\tstart\tend\tperiod\t"+sample+"var\t"+sample+"var2\tconfidence\n");
			for (int i = 0; i < samtoolList.size(); i++){
				TandemRepeatVariant trv = samtoolList.get(i); 
				outS.print(trv.getChr()+"\t"+trv.getTandemRepeat().getID()+"\t"+trv.getStart()+"\t"+trv.getEnd()+"\t" +trv.getTandemRepeat().getPeriod() + "\t" + trv.getVar()+"\t"+trv.getVar2()+"\t" + trv.getConfidence() + "\n");
			}			
		}
		outS.close();
	}
	static double range = 2.0;


	public static void readVCFResults(ArrayList<TandemRepeatVariant> trList, HashMap<String, Integer>  chrIndex, 
		String fName, double qual) throws IOException{
		BufferedReader in = SequenceReader.openFile(fName);
		String line = "";

		int index = 0;
		TandemRepeatVariant tr = trList.get(0);
		//int indels = 0;

		while ( (line = in.readLine()) != null){
			if (line.startsWith("#")){
				//System.out.println(line);
				continue;
			}
			line = line.trim();
			if (line.length() == 0) continue;			
			String [] toks = line.split("\\t");
			if (toks.length < 10) continue;
			//chr1 406 . CA CAA,CAAAAA 30.6 . INDEL;DP=52 GT:PL:GQ 1/1:79,25,8,101,0,69:31

			double thisQual = 0;
			//this is to account for varscan
			if (!".".equals(toks[5]))
				thisQual =  Double.parseDouble(toks[5]);

			if (thisQual < qual) continue;
			if (chrIndex.get(toks[0]) == null)
				continue;

			if (chrIndex.get(toks[0]) < chrIndex.get(tr.getChr())) 
				continue;			

			int position = Integer.parseInt(toks[1]);			

			while (	(chrIndex.get(toks[0]) == chrIndex.get(tr.getChr()) && position > tr.getEnd()) 
				|| chrIndex.get(toks[0]) > chrIndex.get(tr.getChr() )){
				index ++;
				if (index >= trList.size()){
					break;
				}				
				tr = trList.get(index);
			}

			if (index >= trList.size()){
				break;
			}
			//we are done

			if (chrIndex.get(toks[0]) < chrIndex.get(tr.getChr())) 
				continue;

			int refLen = toks[3].length();
			if (chrIndex.get(toks[0]) == chrIndex.get(tr.getChr()) && position + refLen < tr.getStart())
				continue;

			//			indels++;


			String [] pToks = toks[4].split(",");

			String [] types = toks[9].split(":")[0].split("/|\\|");			

			if (types.length < 2)
				throw new RuntimeException("There are " + types.length + " types " + line);

			if (types[0].equals("1")){
				tr.setVar(tr.getVar() + (pToks[0].length() * 1.0 - refLen)/tr.getPeriod());
			}
			if (types[0].equals("2")){
				tr.setVar(tr.getVar() + (pToks[1].length() * 1.0 - refLen)/tr.getPeriod());
			}

			if (types[1].equals("1")){				
				tr.setVar2(tr.getVar2() + (pToks[0].length() * 1.0 - refLen)/tr.getPeriod());
			}
			if (types[1].equals("2")){
				if (pToks.length > 1)				
					tr.setVar2(tr.getVar2() + (pToks[1].length() * 1.0 - refLen)/tr.getPeriod());
				else
					tr.setVar2(tr.getVar2() + (pToks[0].length() * 1.0 - refLen)/tr.getPeriod());
			}

			//Convert to probability
			if (tr.getConfidence() == 0)
				tr.setConfidence(1.0 - JapsaMath.phred2prob(thisQual));			
			else
				tr.setConfidence((tr.getConfidence() + 1.0 - JapsaMath.phred2prob(thisQual))/2);			
		}	
	}




}
