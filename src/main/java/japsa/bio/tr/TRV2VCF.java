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
 * 20/06/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.bio.tr;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.JapsaMath;
import japsa.util.deploy.Deployable;

import java.io.IOException;
import java.util.ArrayList;



/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
@Deployable(scriptName = "jsa.trv.trv2vcf",
            scriptDesc = "Convert tandem repeat variation in trv format to vcf")
public class TRV2VCF {
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {		
		/*********************** Setting up script ****************************/		 
		String scriptName = "jsa.trv.trv2vcf";
		String desc = "Convert tandem repeat variation in trv format to vcf\n";		
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [params]");
		/**********************************************************************/

		cmdLine.addString("input", null,
				"Name of the input file in trv, - for standard in",true);		
		cmdLine.addString("reference", null,
				"File containing the reference genome",true);		
		cmdLine.addStdAlphabet();		
		cmdLine.addString("output", "-",
				"Name of the output file ( - for standard output)");


		cmdLine.addStdHelp();	
		/**********************************************************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")){
			System.out.println(desc + cmdLine.usage());			
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usage());
			System.exit(-1);
		}	
		/**********************************************************************/
		String inFile = cmdLine.getStringVal("input");
		String refFile = cmdLine.getStringVal("reference");		
		String outFile = cmdLine.getStringVal("output");
		Alphabet alphalbet = Alphabet.getAlphabet(cmdLine.getStringVal("dna"));
		
		if ("-".equals(refFile)){
			System.err.println("ERROR: reference must be from some files");
			System.err.println(cmdLine.usage());
			System.exit(-1);
		}
		ArrayList<TandemRepeatVariant> trvList = TandemRepeatVariant.readFromFile(inFile);
		
		
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(outFile);
		
		out.print("##fileformat=VCFv4.0\n");
		//out.write("##fileDate="+ System. );
		out.print("##source=STRViper\n");
		out.print("##reference="+refFile+"\n");

		out.print("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
		out.print("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
		out.print("##FILTER=<ID=q3,Description=\"Quality below 3 (~50%)\">\n");
		out.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");


		
		SequenceReader reader = SequenceReader.getReader(refFile);
		Sequence seq = reader.nextSequence(alphalbet);
		for (int x = 0; x <trvList.size(); x++){			
			TandemRepeatVariant trv = trvList.get(x);

			//make sure the same sequence is retrieved
			while (!trv.getChr().equals(seq.getName()))
				seq = reader.nextSequence(alphalbet);

			//converting variation format
			int period = trv.getPeriod();



			int nuc = (int) Math.round(trv.var * period);
			if (nuc == 0)
				continue;//no variation

			int start = trv.getStart();
			int end = trv.getEnd();

			String unit = trv.getTandemRepeat().getUnit();

			if (unit == null || unit.length() == 0){
				StringBuilder sb = new StringBuilder(period);
				for (int i = start - 1;i < start + period -1;i++){
					sb.append(seq.charAt(i));
				}
				unit = sb.toString();
			}			

			int qual = (int) 
					(JapsaMath.prob2phred(1 - trv.confidence));			

			if (nuc < 0){//deletion
				//remove of nuc nucleotides at the end
				StringBuilder ref = new StringBuilder();				
				for (int i = end + nuc - 1; i < end; i ++ ){
					ref.append(seq.charAt(i));					
				}				
				out.print(trv.getChr() + '\t' + (end + nuc) +"\t.\t"+ref+"\t"+seq.charAt(end + nuc - 1)+"\t" 
						+  qual + "\t"  //qual 
						+ (qual>3?"PASS":"q3")+"\t" //filter
						+ "NS=" + trv.evidence + "\t" 
						+"\n");				
			}else{//insertion
				StringBuilder alt = new StringBuilder();
				//ref = japsa.seq.charAt(end - 1)
				alt.append(seq.charAt(end - 1));//ref
				int ind = (end - start + 1) % period;
				for (int i = 0; i < nuc; i++){
					alt.append(unit.charAt(ind));
					ind ++;

					if (ind >= unit.length())
						ind = 0;
				}//for
				out.print(trv.getChr() + '\t' + (end + nuc) +"\t.\t"+seq.charAt(end- 1)+'\t' + alt + "\t"
						+  qual + "\t"  //qual 
						+ (qual>3?"PASS":"q3")+"\t" //filter
						+ "NS=" + trv.evidence + "\t" 
						+"\n");						
			}//else
		}
		out.close();	
	}

}
