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
 * 23/07/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.tr;

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;




/**
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 *
 */
public class TandemRepeatVariant implements Comparable<TandemRepeatVariant>{
	private static final Logger LOG = LoggerFactory.getLogger(TandemRepeatVariant.class);


	public static String 
	varHd    = "var",                 //06
	var2Hd    = "var2",                 //06
	confidenceHd = "confidence",      //07

	meanHd   = "mean",                //08
	stdHd    = "std",                 //09
	eviHd  = "evidence",             //10
	evi2Hd  = "evidence2" ;             //10


	//public static String[] STANDARD_HEADERS  = {TandemRepeat.chrHd, TandemRepeat.idHd, TandemRepeat.startHd, TandemRepeat.endHd, TandemRepeat.periodHd, TandemRepeat.unitNoHd, varHd, confidenceHd, meanHd, stdHd,eviHd}; 
	//public static String[] STANDARD_HEADERS2 = {TandemRepeat.chrHd, TandemRepeat.idHd, TandemRepeat.startHd, TandemRepeat.endHd, TandemRepeat.periodHd, TandemRepeat.unitNoHd, varHd, meanHd, stdHd};

	public static String[] STANDARD_HEADERS  = 
		{TandemRepeat.chromHd, TandemRepeat.idHd, TandemRepeat.startHd, 
				TandemRepeat.endHd, TandemRepeat.periodHd, TandemRepeat.unitNoHd, 
				varHd, confidenceHd, stdHd, eviHd};

	public static String[] STANDARD_HEADERS2 = 
		{TandemRepeat.chromHd, TandemRepeat.idHd, TandemRepeat.startHd, 
				TandemRepeat.endHd, TandemRepeat.periodHd, 
				TandemRepeat.unitNoHd, varHd, stdHd};

	public static String[] SIMPLE_HEADERS  = 
		{TandemRepeat.chromHd, TandemRepeat.idHd, TandemRepeat.startHd, 
				TandemRepeat.endHd, TandemRepeat.periodHd, TandemRepeat.unitNoHd, 
				varHd, eviHd};

	public static String[] SIMPLE_HEADERS2  = 
		{TandemRepeat.chromHd, TandemRepeat.idHd, TandemRepeat.startHd, 
				TandemRepeat.endHd, TandemRepeat.periodHd, TandemRepeat.unitNoHd, 
				varHd, eviHd, var2Hd, evi2Hd};

	TandemRepeat tandemRepeat;
	/**
	 * @return the tandemRepeat
	 */
	public TandemRepeat getTandemRepeat() {
		return tandemRepeat;
	}

	/**
	 * @param tr the tandemRepeat to set
	 */

	public void setTandemRepeat (TandemRepeat tr) {
		this.tandemRepeat = tr;
	}

	double var; //the variations
	double var2;//second allele
	double confidence, evidence;//The confident is in probability (i.e., 1-10^phred
	double evidence2;

	@Deprecated
	double mean = 0;//a bland distribution


	@Deprecated
	public double std = 10;

	//public String start;


	public TandemRepeatVariant(){
		tandemRepeat = new TandemRepeat();
	}

	/**
	 * @return the confidence
	 */
	public double getConfidence() {
		return confidence;
	}
	/**
	 * @param confidence the confidence to set
	 */
	public void setConfidence(double confidence) {
		this.confidence = confidence;
	}
	/**
	 * @return the mean
	 * @Deprecated: mean will be removed in the new future
	 */
	@Deprecated
	public double getMean() {
		return mean;
	}
	/**
	 * @param mean the mean to set
	 *  @Deprecated: mean will be removed
	 */
	@Deprecated
	public void setMean(double mean) {
		this.mean = mean;
	}
	/**
	 * @return the std
	 */
	public double getStd() {
		return std;
	}
	/**
	 * @param std the std to set
	 */
	public void setStd(double std) {
		this.std = std;
	}

	public void swapVar(){
		double tmp = var;var = var2;var2 = tmp;
	}

	/**
	 * @return the var
	 */
	public double getVar() {
		return var;
	}

	public double getVar2() {
		return var2;
	}
	/**
	 * @param var the var to set
	 */
	public void setVar(double var) {
		this.var = var;
	}
	public void setVar2(double var2) {
		this.var2 = var2;
	}


	public TandemRepeatVariant tandemRepeatClone(){
		TandemRepeatVariant trv = new TandemRepeatVariant();

		trv.tandemRepeat = this.tandemRepeat;
		//trv.tandemRepeat.schr     = chr;
		//trv.start   = start;
		//trv.end     = end;
		//trv.period  = period;
		//trv.unitNo  = this.unitNo;

		trv.var     = this.var;
		trv.confidence = this.confidence;
		trv.mean    = this.mean;
		trv.std     = this.std;
		trv.confidence     = this.confidence;

		return trv;
	}

	/**
	 * Read from a line and a list of fields
	 * @param line
	 * @param hds
	 * @return
	 */
	public static TandemRepeatVariant read(String line, String [] hds){		
		TandemRepeatVariant rec = new TandemRepeatVariant();
		String [] toks = line.trim().split("\\t");

		for (int i = 0; i < hds.length; i++ ){
			if (TandemRepeat.chromHd.equals(hds[i]))
				rec.tandemRepeat.setChr(toks[i]);
			else if (TandemRepeat.idHd.equals(hds[i]))
				rec.tandemRepeat.setID(toks[i]);	
			else if (TandemRepeat.startHd.equals(hds[i]))
				rec.tandemRepeat.setStart(Integer.parseInt(toks[i]));
			else if (TandemRepeat.endHd.equals(hds[i]))
				rec.tandemRepeat.setEnd(Integer.parseInt(toks[i]));
			else if (TandemRepeat.periodHd.equals(hds[i]))
				rec.tandemRepeat.setPeriod(Integer.parseInt(toks[i]));
			else if (TandemRepeat.unitNoHd.equals(hds[i]))
				rec.tandemRepeat.setUnitNo(Double.parseDouble(toks[i]));
			else if (varHd.equals(hds[i]))
				rec.var = Double.parseDouble(toks[i]);
			else if (var2Hd.equals(hds[i]))
				rec.var2 = Double.parseDouble(toks[i]);
			else if (confidenceHd.equals(hds[i]))
				rec.confidence = Double.parseDouble(toks[i]);
			else if (meanHd.equals(hds[i]))
				rec.mean = Double.parseDouble(toks[i]);
			else if (stdHd.equals(hds[i]))
				rec.std = Double.parseDouble(toks[i]);
			else if (eviHd.equals(hds[i]))
				rec.evidence = Double.parseDouble(toks[i]);
		}		
		return rec;

	}
	/**
	 * Print to a line
	 * @param hds
	 * @return
	 */
	public String toString(String [] hds){		
		StringBuffer sb = new StringBuffer();		
		for (int i = 0; i < hds.length; i++ ){
			if (i > 0) 
				sb.append("\t");

			if (TandemRepeat.chromHd.equals(hds[i]))
				sb.append(this.tandemRepeat.getChr());
			else if (TandemRepeat.idHd.equals(hds[i]))
				sb.append(tandemRepeat.getID());
			else if (TandemRepeat.startHd.equals(hds[i]))
				sb.append(this.tandemRepeat.getStart());
			else if (TandemRepeat.endHd.equals(hds[i]))
				sb.append(this.tandemRepeat.getEnd());
			else if (TandemRepeat.periodHd.equals(hds[i]))
				sb.append(this.tandemRepeat.getPeriod());
			else if (TandemRepeat.unitNoHd.equals(hds[i]))
				sb.append(this.tandemRepeat.getUnitNo());			
			else if (varHd.equals(hds[i]))
				sb.append(this.var);
			else if (var2Hd.equals(hds[i]))
				sb.append(this.var2);
			else if (confidenceHd.equals(hds[i]))
				sb.append(this.confidence);
			else if (meanHd.equals(hds[i]))
				sb.append(this.mean);
			else if (stdHd.equals(hds[i]))
				sb.append(this.std);
			else if (eviHd.equals(hds[i]))
				sb.append(this.evidence);
			else if (evi2Hd.equals(hds[i]))
				sb.append(this.evidence2);
		}		
		return sb.toString();		
	}


	public int getEnd() {
		return tandemRepeat.getEnd();
	}

	public int getStart() {
		return tandemRepeat.getStart();
	}

	public String getChr() {
		return tandemRepeat.getChr();
	}

	public int getPeriod() {
		return tandemRepeat.getPeriod();
	}

	public void addEvidence(double moreEvi){
		evidence += moreEvi;
	}
	public void addEvidence2(double moreEvi){
		evidence2 += moreEvi;
	}	

	public double getEvidence(){
		return evidence;
	}

	public double getEvidence2(){
		return evidence2;
	}

	public static ArrayList<TandemRepeatVariant> readFromFile(String fileName) throws IOException{
		//Start with the default header
		String[] headers = STANDARD_HEADERS2;

		ArrayList<TandemRepeatVariant> trfList = new ArrayList<TandemRepeatVariant>(); 
		BufferedReader in = SequenceReader.openFile(fileName);
		String line = "";
		while ((line = in.readLine()) != null){
			line = line.trim();
			if (line.length() == 0) continue;

			if (line.startsWith("#H:"))
				headers = line.substring(3).split("\\t");
			else if (line.startsWith("#"))
				continue;
			else
				trfList.add(TandemRepeatVariant.read(line, headers));		

		}//while

		LOG.info("Read in " + trfList.size() + " TRs");
		return trfList;
	}

	/**
	 * Read from annotation of TRF
	 * @param annoList
	 * @return
	 */
	public static ArrayList<TandemRepeatVariant> readFromBio(ArrayList<JapsaAnnotation> annoList){		
		ArrayList<TandemRepeatVariant>		trvList = new ArrayList<TandemRepeatVariant>();

		for (int idx = 0; idx < annoList.size(); idx ++){
			JapsaAnnotation anno = annoList.get(idx);
			for (int x = 0; x < anno.numFeatures(); x++){

				JapsaFeature feature = anno.getFeature(x);

				if (!feature.getType().equals("trf"))
					continue;

				TandemRepeatVariant trv = new TandemRepeatVariant();
				trv.tandemRepeat.setChr(anno.getAnnotationID());
				trv.tandemRepeat.setStart(Integer.parseInt(feature.getID().substring(1)));

				trv.mean = 0;
				trv.std = 10;

				Iterator<String> trIter = feature.getDescStr();
				while (trIter.hasNext()){
					String desc = trIter.next();
					if (desc.startsWith("@R:"))
						trv.tandemRepeat.setPeriod(Integer.parseInt(desc.substring(3)));
					else if (desc.startsWith("@DIF:")){
						trv.var = Integer.parseInt(desc.substring(5));
						trv.mean = trv.var;
						trv.std = 0;//i am pretty sure about this
					}
				}
				trv.tandemRepeat.setStart(Integer.parseInt(feature.getID().substring(1)));
				trv.tandemRepeat.setEnd((int) (trv.tandemRepeat.getStart() + feature.getLength() - trv.tandemRepeat.getPeriod() * trv.var - 1));
				trv.tandemRepeat.setUnitNo((feature.getLength() - trv.tandemRepeat.getPeriod() * trv.var) / trv.tandemRepeat.getPeriod());


				trvList.add(trv);
			}			
		}	

		return trvList;
	}

	public static void write(ArrayList<TandemRepeatVariant> trvList, SequenceOutputStream out) throws IOException{
		print(trvList, out, STANDARD_HEADERS);
	}

	public static void print(ArrayList<TandemRepeatVariant> trvList, SequenceOutputStream out, String [] headers) throws IOException{
		printHeader(out,headers);
		for (int i = 0; i < trvList.size(); i++){
			out.print(trvList.get(i).toString(headers));	
			out.print('\n');
		}
	}

	public static void printHeader(SequenceOutputStream out, String [] headers) throws IOException{
		out.print("#H:" + headers[0]);
		for (int i=1; i < headers.length; i++)
			out.print("\t"+headers[i]);
		out.print('\n');
	}



	public String toString(){
		return tandemRepeat.toString();
	}


	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(TandemRepeatVariant o) {
		return tandemRepeat.compareTo(o.tandemRepeat);
	}

	/**
	 * Write the variant to bed file. Currently use only the first heterozygous 
	 * variant is written out.
	 * TODO: Handle heterozygous case 
	 * @param out
	 * @throws IOException
	 */
	public void writeBED(SequenceOutputStream out) throws IOException{
		out.print((tandemRepeat.getParent()+'\t' + (getStart()-1) + '\t' + getEnd()+
				'\t' + tandemRepeat.getUnit() + '\t' + this.var  + '\t' + (tandemRepeat.getStrand() == '-'?'-':'+') + '\t' + tandemRepeat.getUnitNo() +'\n'));
	}
}
