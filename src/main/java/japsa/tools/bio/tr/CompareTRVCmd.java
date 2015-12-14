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

package japsa.tools.bio.tr;

import japsa.bio.tr.TandemRepeatVariant;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;


/**
 * @author minhduc
 * 
 */
@Deployable(
	scriptName = "jsa.str.strvcompare",
	scriptDesc = "Compare TR variation to an answer file, use to evaluate strv"
	)
public class CompareTRVCmd  extends CommandLine{	
	public CompareTRVCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("ansFile", null,	"Name of the answer file",true);
		addInt("min", 0, "minimum length of TR region");
		addInt("max", 10000, "maximum length of TR region");
		addString("strviper", null,	"Name of the str file");
		//cmdLine.addString("lob", null,	"Name of the lobstr file");		
		addString("vcf", null,	"Name of the results in vcf format (samtools and Dindel)");
		addBoolean("baseline", false, "Benchmark against a baseline");	

		addStdHelp();		
	}  
//FIXME: standardise
	public static void main(String[] args) throws Exception {		
		CommandLine cmdLine = new CompareTRVCmd();		
		args = cmdLine.stdParseLine(args);

		String  ans = cmdLine.getStringVal("ansFile");
		String  lob = cmdLine.getStringVal("lob");
		String  strviper = cmdLine.getStringVal("strviper");
		String  vcf = cmdLine.getStringVal("vcf");	
		int min = cmdLine.getIntVal("min");
		int max = cmdLine.getIntVal("max");		
	
		// getInsertSizeSlow(samFile, output);
		ArrayList<TandemRepeatVariant> ansList = TandemRepeatVariant.readFromFile(ans);
		if (cmdLine.getBooleanVal("baseline")){
			System.out.println("###baseline ");
			compareStr(ansList, new ArrayList<TandemRepeatVariant>(),min, max);
		}
		
		
		if (lob != null){
			ArrayList<TandemRepeatVariant> lobList = new ArrayList<TandemRepeatVariant> ();
			BufferedReader in = SequenceReader.openFile(lob);
			String line = "";
			while ( (line = in.readLine()) != null){
				line = line.trim();
				if (line.length() <= 0) continue;
				if (line.startsWith("#"))
					continue;
				
				lobList.add(readLobSTR(line));		
			}
			System.out.println("###LobSTR " + lob);
			compareStr(ansList, lobList,min,max);
		}	
		

		if (strviper != null){
			ArrayList<TandemRepeatVariant> myList =TandemRepeatVariant.readFromFile(strviper);
			System.out.println("###STRViper " + strviper);
			compareStr(ansList, myList,min,max);
		}
		if (vcf != null){
			double [] quals = {0, 20,  40, 60, 80, 100, 120, 140, 160,180,200, 220, 240, 260, 280, 300, 320, 340, 360, 380,400,500,600,700,800,900,1000};
			int tt = 0;
			//for ( tt = 0; tt < quals.length; tt++){
				ArrayList<TandemRepeatVariant> samtoolList = new ArrayList<TandemRepeatVariant>(ansList.size());
				for (int i = 0; i < ansList.size(); i++){
					TandemRepeatVariant trv = ansList.get(i).tandemRepeatClone();
				
					trv.setMean(0);
					trv.setVar(0);
					trv.setConfidence(0);
					trv.setStd(0);
				
					samtoolList.add(trv);
				}	
				System.out.println("###SAMtool " + vcf);
				readVCFResults(samtoolList, vcf, quals[tt]);
				compareStr(ansList, samtoolList,min,max);
			//ShortTandemRepeat.write(samtoolList, System.out, ShortTandemRepeat.STANDARD_HEADERS);
			//}
		}
	}
	
	
	static double range = 2.0;
	
	public static void readModilResults(ArrayList<TandemRepeatVariant> trList,String fName) throws IOException{
		BufferedReader in = SequenceReader.openFile(fName);
		String line = "";
		
		int index = 0;
		TandemRepeatVariant tr = trList.get(0);
		int indels = 0;
			
		while ( (line = in.readLine()) != null){
			if (line.startsWith("#")) continue;
			line = line.trim();			
			if (line.length() < 0) continue;			
			String [] toks = line.split("\\t");
			if (toks.length < 11) continue;
			
			String chr = toks[1];
			if (chr.startsWith("chr"))
				chr = chr.substring(3);	
			
			if (chr.compareTo(tr.getChr()) < 0){
				System.out.println("Pass  1 " + line);
				continue;//already advanced to the next chromosome			
			}			
			
			int positionS = Integer.parseInt(toks[2]);
			int positionE = Integer.parseInt(toks[3]);
			
			
			while ((chr.compareTo(tr.getChr()) == 0 && positionS > tr.getEnd()) || (chr.compareTo(tr.getChr()) > 0) ){
				index ++;
				if (index >= trList.size()){
					break;
				}				
				tr = trList.get(index);
			}
			
			if (index >=  trList.size()){
				System.out.println("Out  " + line);
				break;//while
			}			
			//assert: chr == tr.chr && positionS <= tr.end			
			if (chr.compareTo(tr.getChr()) == 0 && positionE < tr.getStart()){
				System.out.print("Pass  2 " + line + ":" + tr.getChr() + ", " + tr.getStart() +" -> "+tr.getEnd());
				if (index > 0) System.out.print("#" + trList.get(index - 1).getChr()+", "+trList.get(index - 1).getStart()+" -> "+trList.get(index - 1).getEnd());
				System.out.println();
				continue;
			}
			
			indels++;
			
			int overlap = positionS;
			if (tr.getStart() > overlap) overlap = tr.getStart();
			
			if (tr.getEnd() < positionE)
				overlap = tr.getEnd() - overlap;
			else
				overlap = positionE - overlap;
			
			double indelSize = Double.parseDouble(toks[6]);
			
			System.out.println("Enter : " + line + "   " + overlap + "   vs   " + indelSize);			
						
			tr.setVar(tr.getVar() + indelSize/tr.getPeriod());
			tr.setMean(tr.getVar());
			tr.std = 2;
		}
		
		System.out.println("## " + indels + " indels found ");
		
	}
	
	

	
	
	public static void readVCFResults(ArrayList<TandemRepeatVariant> trList,String fName, double qual) throws IOException{
		BufferedReader in = SequenceReader.openFile(fName);
		String line = "";
		
		int index = 0;
		TandemRepeatVariant tr = trList.get(0);
		int indels = 0;
			
		while ( (line = in.readLine()) != null){
			if (line.startsWith("#")) continue;
			line = line.trim();
			if (line.length() == 0) continue;			
			String [] toks = line.split("\\t");
			if (toks.length < 8) continue;
			//if (!toks[7].startsWith("INDEL"))
			//	continue;
			
			double thisQual = 0;
			//this is to account for varscan
			if (!".".equals(toks[5]))
				thisQual =  Double.parseDouble(toks[5]);
			
			if (thisQual < qual) continue;						
			
			if (toks[0].compareTo(tr.getChr()) < 0) continue;
			int position = Integer.parseInt(toks[1]);			
			
			while ((toks[0].compareTo(tr.getChr()) == 0 && position > tr.getEnd()) || (toks[0].compareTo(tr.getChr()) > 0) ){
				index ++;
				if (index >= trList.size()){
					break;
				}				
				tr = trList.get(index);
			}
			int refLen = toks[3].length();
			if (toks[0].compareTo(tr.getChr()) == 0 && position + refLen < tr.getStart())
				continue;
			
			indels++;
			
			String [] pToks = toks[4].split(",");
			
			int sumL = 0, countL = 0;
			for (int x = 0; x < pToks.length; x++){
				if (pToks[x].startsWith("-")){
					sumL -= (pToks[x].length() - 1) -refLen;
				}else if (pToks[x].startsWith("+")){
					sumL += (pToks[x].length() - 1) +refLen;
				}else				
					sumL += pToks[x].length();
				countL ++;
			}
			tr.setVar(tr.getVar() + (sumL * 1.0 / countL - refLen)/tr.getPeriod());
		}
		
		System.out.println("## " + indels + " indels found ");		
	}
	
	public static void readSDIVar(ArrayList<TandemRepeatVariant> trList, String fName) throws IOException{
		BufferedReader in = SequenceReader.openFile(fName);
		String line = "";
		
		int index = 0;
		TandemRepeatVariant tr = trList.get(0);
		int indels = 0;
			
		while ( (line = in.readLine()) != null){
			if (line.startsWith("#")) continue;
			line = line.trim();
			if (line.length() == 0) continue;			
			String [] toks = line.split("\\t");
			if (toks.length < 5) continue;
			//if (!toks[7].startsWith("INDEL"))
			//	continue;
			
			if (toks[0].compareTo(tr.getChr()) < 0) continue;
			int position = Integer.parseInt(toks[1]);			
			
			while ((toks[0].compareTo(tr.getChr()) == 0 && position > tr.getEnd()) || (toks[0].compareTo(tr.getChr()) > 0) ){
				index ++;
				if (index >= trList.size()){
					break;
				}				
				tr = trList.get(index);
			}
			int refLen = toks[3].length();
			if (toks[0].compareTo(tr.getChr()) == 0 && position + refLen < tr.getStart())
				continue;
			
			
			boolean foundIndel = false;
			
			String [] pToks = toks[4].split(",");
			int sumL = 0, countL = 0;
			for (int x = 0; x < pToks.length; x++){
				sumL += pToks[x].length();
				countL ++;
				foundIndel = true;
			}
			tr.setVar(tr.getVar() + (sumL * 1.0 / countL - refLen)/tr.getPeriod());
			
			if (foundIndel)
				indels++;
		}
		
		System.out.println("## " + indels + " indels found ");
		
	}
	
	public static void combineStr(ArrayList<TandemRepeatVariant> trList, ArrayList<TandemRepeatVariant> pList) throws IOException{
		int pInd = 0;
		//double SSE = 0;
		//int sRight = 0, rRight = 0, vRight = 0, pRight = 0;
		int pCount = 0;

		TandemRepeatVariant pStr = null;
		for (int i = 0; i < trList.size(); i++){
			TandemRepeatVariant ans = trList.get(i);			

			if (pStr == null && pInd < pList.size()){
				pStr = pList.get(pInd);
				pInd ++;				
			}

			if (pStr != null){
				if (ans.getChr().equals(pStr.getChr()) && ans.getStart() == pStr.getStart() && ans.getEnd() == pStr.getEnd()){
					ans.setVar(pStr.getVar());					
					ans.setMean(ans.getVar());
					ans.setConfidence(pStr.getConfidence());
					ans.std = 2.0;
					pCount ++;
					pStr = null;					
				}				
			}
		}
		System.out.println("##Combine " + pCount);

		
	}
	
	public static void compareStr(ArrayList<TandemRepeatVariant> ansList, ArrayList<TandemRepeatVariant> pList, int minLen, int maxLen) throws IOException{
		int pInd = 0;
		double SSE = 0;
		int sRight = 0, rRight = 0, vRight = 0, pRight = 0;
		//sRight = strict count (only equal excepted)
		//rRight = relax count (if |predicted - answer| < range = 2)
		//vRight = count if there is an event or not
		//pRight = countTrue if correctly predict if there is an insertion/or a deletion
		int pCount = 0;
		
		int delTP = 0, delTN = 0, delFP = 0, delFN = 0;  //deletion
		int insTP = 0, insTN = 0, insFP = 0, insFN = 0;  //insertion
		int idTP = 0, idTN = 0, idFP = 0, idFN = 0    ;  //indels in general		

		TandemRepeatVariant pStr = null;
		int total = 0;
		for (int i = 0; i < ansList.size(); i++){
			TandemRepeatVariant ans = ansList.get(i);
			if (ans.getEnd() - ans.getStart() + 1 < minLen) continue;
			if (ans.getEnd() - ans.getStart()+ 1 > maxLen) continue;
			
			total ++;
			
			double predicted = 0.0;

			//if (pStr == null && pInd < pList.size()){
			//	pStr = pList.get(pInd);
			//	pInd ++;				
			//}
			
			while ((pStr == null || ans.compareTo(pStr) > 0) && (pInd < pList.size())){
				pStr = pList.get(pInd);
				pInd ++;
			}

			if (pStr != null && ans.compareTo(pStr) == 0){
				//if (ans.chr.equals(pStr.chr) && ans.start == pStr.start && ans.end == pStr.end){
					predicted = pStr.getVar();
					pCount ++;
					pStr = null;					
				//}
			}

			if (Math.round(predicted) == Math.round(ans.getVar()))
				sRight ++;	

			if (Math.round(predicted) == Math.round(ans.getVar()))
				pRight ++;
			else if (Math.round(predicted) * Math.round(ans.getVar()) > 0)
				pRight ++;


			if (Math.round(predicted) == Math.round(ans.getVar()))
				vRight ++;
			else if (Math.round(predicted) * Math.round(ans.getVar()) != 0)
				vRight ++;			

			if (Math.abs(predicted - ans.getVar()) <= range)
				rRight ++;

			SSE += (predicted - ans.getVar()) * (predicted - ans.getVar());
			
			if (ans.getVar() > 0.5){//an insertion				
				if (predicted > 0.5){//predicted as an insertion
					insTP ++;					
					delTN ++;
					idTP ++;
				}else if (predicted < -0.5){//predicted as a deletion
					insFN ++;//
					delFP ++;
					idTP ++;
				}else{//predicted as no indel
					insFN ++;
					delTN ++;
					idFN++;
				}				
			}else if (ans.getVar() < -0.5){//a deletion
				if (predicted > 0.5){//predicted as an insertion
					insFP ++;					
					delFN ++;
					idTP ++;
				}else if (predicted < -0.5){//predicted as a deletion
					insTN ++;//
					delTP ++;
					idTP ++;
				}else{//predicted as no indel
					insTN ++;
					delFN ++;
					idFN++;
				}			
			}else{//no indels
				if (predicted > 0.5){//predicted as an insertion
					insFP ++;					
					delTN ++;
					idFP ++;
				}else if (predicted < -0.5){//predicted as a deletion
					insTN ++;//
					delFP ++;
					idFP ++;
				}else{//predicted as no indel
					insTN ++;
					delTN ++;
					idTN++;
				}				
			}
		}
		//int total = ansList.size();

		double insSn = 1.0* insTP/(insTP + insFN); //= recall
		double insSp = 1.0* insTN/(insTN + insFP);
		double insPv = 1.0* insTP/(insTP + insFP);
		double insAc = 1.0* (insTP + insTN)/(insTP + insTN + insFP + insFN);
		double insF1 =   2* insPv * insSn / (insPv + insSn);
		
		double delSn = 1.0* delTP/(delTP + delFN);
		double delSp = 1.0* delTN/(delTN + delFP);
		double delPv = 1.0* delTP/(delTP + delFP);
		double delAc = 1.0* (delTP + delTN)/(delTP + delTN + delFP + delFN);
		double delF1 =   2* delPv * delSn / (delPv + delSn);
		
		double idSn  = 1.0* idTP/(idTP + idFN);
		double idSp  = 1.0* idTN/(idTN + idFP);
		double idPv  = 1.0* idTP/(idTP + idFP);
		double idAc  = 1.0* (idTP + idTN)/(idTP + idTN + idFP + idFN);
		double idF1  =   2* idPv * idSn / (idPv + idSn);		
		
		System.out.println("##Compare " + pCount + " records ( out of " + pList.size() + ") against " + total + " (" + (idTP + idTN + idFP + idFN) +")");
		System.out.printf("%5.2f%%\t%5.2f%%\t%5.2f%%\t%5.2f%%\t%6.4f\t",100.0 * sRight/total, 100.0*rRight/total, 100.0*pRight/total, 100.0*vRight/total, Math.sqrt(SSE/total));
		System.out.printf("%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t",idSp, idSn, idPv, idAc, idF1);
		System.out.printf("%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t",insSp, insSn, insPv, insAc, insF1);
		System.out.printf("%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n",delSp, delSn, delPv, delAc, delF1);	
		
	}


	public static TandemRepeatVariant readLobSTR(String line){
		if (line == null) return null;
		TandemRepeatVariant tr = new TandemRepeatVariant();		
		String [] toks = line.trim().split("\\t");


		tr.getTandemRepeat().setChr(toks[0]);
		tr.getTandemRepeat().setStart(Integer.parseInt(toks[1]));
		tr.getTandemRepeat().setEnd(Integer.parseInt(toks[2]));	
		tr.getTandemRepeat().setPeriod(Integer.parseInt(toks[4]));
		
		String[] alle = toks[6].split(",");
		tr.setVar(((Double.parseDouble(alle[0]) + Double.parseDouble(alle[1])) / tr.getPeriod()) / 2);		
		//tr.setConfidence(confidence);		
		return tr;
	}
}
