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
 * File: Convert2STRV.java
 * 16/11/2013 - Minh Duc Cao: Created
 *
 ****************************************************************************/

package japsadev.tools.misc;

import japsa.bio.tr.TandemRepeatVariant;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFileFormat;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.tools.bio.tr.CompareTRVCmd;
import japsa.util.CommandLine;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Replaced by VCF2STRV
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 */
@Deprecated
public class Convert2STRV {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{

		/*********************** Setting up script ****************************/		 
		String scriptName = "jsa.str.varConvert	";//used to be varConvert but replaced by vcf2strvs
		String desc = "convert to STR variations from various format to the STRV format\n";		
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [params]");

		/**********************************************************************/		
		//cmdLine.addString("isize", null, "Name of the insert size file, - for standard in");
		cmdLine.addString("jsaFile", null,	"Name of the JSA file, this file presents the location of STRs");		
		cmdLine.addString("strFile", null,	"Name of the str file, this file presents the location of STRs (in str format)");
		cmdLine.addString("lobstr",   null,	"Name of the lobstr file");
		cmdLine.addString("vcf",      null,	"Name of file in vcf format (dindel and samtools)");
		cmdLine.addString("sdi",      null,	"Name of file in sdi format");
		cmdLine.addString("modil",    null,	"Name of the modil results (concatenate all the chromosomes together)");
		cmdLine.addString("output", "-" ,	"Name of the output file, - for standard output");

		cmdLine.addBoolean("help", false, "Display this usage and exit");
		/**********************************************************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")){
			System.out.println(desc + cmdLine.usageMessage());			
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usageMessage());
			System.exit(-1);
		}	
		/**********************************************************************/		


		String  strFile = cmdLine.getStringVal("strFile");
		String  bioFile = cmdLine.getStringVal("bioFile");
		String  lob = cmdLine.getStringVal("lob");	
		String  vcf = cmdLine.getStringVal("vcf");
		String  sdi = cmdLine.getStringVal("sdi");
		String  modil = cmdLine.getStringVal("modil");
		String  output = cmdLine.getStringVal("output");

		if ((strFile == null && bioFile == null)){
			System.err.println("ERROR: At least one of 'bioFile' and  'strFile' has to be specified (bioFile has higher precedence)\n"+ cmdLine.usageMessage());			
			System.exit(-1);
		}		

		if (bioFile == null // meaning str get from strFile
				&& lob == null && vcf == null && modil == null && sdi == null){

			System.err.println("ERROR: When str is from strFile, at least one of vcf, lob or modil has be specified\n" + cmdLine.usageMessage());			
			System.exit(-1);			
		}		

		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(output);
		//BufferedOutputStream out;
		//if (output.equals("-"))
		//	out = new BufferedOutputStream(System.out);
		//else
		//	out = new BufferedOutputStream(new FileOutputStream(output));			


		ArrayList<TandemRepeatVariant> strList = null;
		if (bioFile != null){
			JapsaFileFormat reader = new JapsaFileFormat(bioFile);			
			//BioCompFileFormat f = new BioCompFileFormat(SequenceReader.openFile(bioFile));			
			//Iterator<JapsaAnnotation> iter = f.getAnnotationIterator();			
			ArrayList<JapsaAnnotation> annoList = new ArrayList<JapsaAnnotation>();
			
			JapsaAnnotation anno = null;
			while ((anno = reader.readAnnotation()) != null){
				//JapsaAnnotation anno = iter.next();			
				annoList.add(anno);
			}
			strList = TandemRepeatVariant.readFromBio(annoList);
			//for (int i = 0; i < strList.size(); i++){
			//	ShortTandemRepeat str = strList.get(i);
			//}
			reader.close();
		}else{//
			strList = TandemRepeatVariant.readFromFile(strFile);
		}

		//		String[] headers = {chrHd, startHd, endHd, periodHd, unitNoHd, varHd, meanHd, stdHd};		
		//		write(strList, System.out,headers);	

		if (vcf != null){
			ArrayList<TandemRepeatVariant> samtoolList = new ArrayList<TandemRepeatVariant>(strList.size());
			for (int i = 0; i < strList.size(); i++){
				TandemRepeatVariant str = strList.get(i).tandemRepeatClone();				
				//str.mean = 0;
				//str.var = 0;
				//str.confidence = 0;
				//str.std = 10;

				samtoolList.add(str);
			}
			CompareTRVCmd.readVCFResults(samtoolList, vcf, 0);

			for (int i = 0; i < samtoolList.size(); i++){
				TandemRepeatVariant str = samtoolList.get(i);
				str.setMean(str.getVar());
				str.setStd(4.0);//str.std/2;
			}				

			TandemRepeatVariant.print(samtoolList, out, TandemRepeatVariant.STANDARD_HEADERS2);			
			out.close();			
			System.exit(0);
		}
		
		if (sdi != null){
			ArrayList<TandemRepeatVariant> sdiVars = new ArrayList<TandemRepeatVariant>(strList.size());
			for (int i = 0; i < strList.size(); i++){
				TandemRepeatVariant str = strList.get(i).tandemRepeatClone();				
				//str.mean = 0;
				//str.var = 0;
				//str.confidence = 0;
				//str.std = 10;

				sdiVars.add(str);
			}
			CompareTRVCmd.readSDIVar(sdiVars, sdi);

			for (int i = 0; i < sdiVars.size(); i++){
				TandemRepeatVariant str = sdiVars.get(i);
				str.setMean(str.getVar());
				str.setStd(4.0);//str.std/2;
			}				

			TandemRepeatVariant.print(sdiVars, out, TandemRepeatVariant.STANDARD_HEADERS2);			
			out.close();			
			System.exit(0);
		}

		if (modil != null){
			ArrayList<TandemRepeatVariant> modList = new ArrayList<TandemRepeatVariant>(strList.size());
			for (int i = 0; i < strList.size(); i++){
				TandemRepeatVariant str = strList.get(i).tandemRepeatClone();				
				//str.mean = 0;
				//str.var = 0;
				//str.confidence = 0;
				//str.std = 10;

				modList.add(str);
			}
			CompareTRVCmd.readModilResults(modList, modil);			
			TandemRepeatVariant.print(modList, out, TandemRepeatVariant.STANDARD_HEADERS2);


			out.close();			
			System.exit(0);
		}

		if (lob != null){
			ArrayList<TandemRepeatVariant> lobList = new ArrayList<TandemRepeatVariant> ();
			BufferedReader in = SequenceReader.openFile(lob);
			String line = "";
			while ( (line = in.readLine()) != null){
				line = line.trim();
				if (line.length() <= 0) continue;
				if (line.startsWith("#")) continue;
				TandemRepeatVariant str = CompareTRVCmd.readLobSTR(line);
				//str.setMean(str.getVar());
				//str.setStd(2.0);
				lobList.add(str);		
			}
			System.out.println("###LobSTR " + lob);

			CompareTRVCmd.combineStr(strList, lobList);

			TandemRepeatVariant.print(strList, out, TandemRepeatVariant.STANDARD_HEADERS);

			out.close();

			System.exit(0);
		}

		//now write
		TandemRepeatVariant.print(strList, out, TandemRepeatVariant.STANDARD_HEADERS);	
		out.close();	

	}

}
