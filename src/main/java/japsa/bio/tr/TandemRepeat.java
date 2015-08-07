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
 * 13/09/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/


package japsa.bio.tr;

import japsa.seq.JapsaFeature;
import japsa.seq.SequenceOutputStream;
import japsa.seq.XAFReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;



/**
 * Represent the genomic location of a short tandem repeat
 * @author Minh Duc Cao (minhduc \dot cao \at gmail \dot com)
 *
 */
public class TandemRepeat extends JapsaFeature{

	/**
	 * Header
	 */
	public static String	
	chromHd    = "chrom",                 //01
	idHd    = "ID",
	startHd  = "start",               //02
	endHd    = "end",                 //03
	periodHd = "period",              //04
	unitHd   = "unit",                //05
	unitNoHd = "unitNo",              //06

	scoreHd    =   "score",           //07
	entropyHd  =   "entropy",         //08
	pMatchHd   =   "pMatch",          //09
	pIndelHd   =   "pIndel",          //10
	annotationHd   =   "Annotations"      //11
	;
	
	
	/**
	 * This field is deprecated, changes to chrom
	 */
	@Deprecated
	public static String chrHd    = "chr";                 //01
	//Property of a TR -- start, end are with regards to the reference genome	
	private int     period = 2;
	private double unitNo; //the length number of units on the ref
	private String  unit;
	private double entropy, pMatch, pIndel;	
	private String annotations="";
	
	
	//String repSeq = "";	
	public static String [] STANDARD_HEADER={chromHd, idHd, startHd, endHd, periodHd, unitNoHd, scoreHd};
	public static String [] FULL_HEADER={chromHd, idHd, startHd, endHd, periodHd, unitHd, unitNoHd, scoreHd, entropyHd, pMatchHd,  pIndelHd};
	public static String [] FULL_HEADER_WITH_ANNOTATIONS={chromHd, idHd, startHd, endHd, periodHd, unitHd, unitNoHd, scoreHd, entropyHd, pMatchHd,  pIndelHd, annotationHd};
	
	/**
	 * @param chr
	 * @param start
	 * @param end
	 * @param period
	 * @param unitNo
	 */
	public TandemRepeat(String chr, int start, int end, int period,
			double unitNo) {		
		super(start, end, "TR", "", '+', chr);	
		
		this.period = period;
		this.unitNo = unitNo;
	}

	public TandemRepeat() {		
		super(0,0, "TR", "", '+', "");	
	}
	
	public TandemRepeat(String chr, int start, int end) {		
		super(start, end, "TR", "", '+', chr);		
	}

	
	/**
	 * Get short tandem repeat from a biofeature
	 * @param id
	 * @param f
	 */
	public TandemRepeat(JapsaFeature f) {
		super(f.getStart(), f.getEnd(), "TR", f.getID(), f.getStrand(), f.getParent());		
		
		Iterator<String> iter = f.getDescStr();
		while (iter.hasNext()) {
			String desc = iter.next();
			this.addDesc(desc);
			if (!desc.startsWith("@"))
				continue;
			
			if (desc.startsWith("@@")){
				this.annotations = desc.substring(1);
				continue;//while
			}
				
			String[] toks = desc.trim().split(":");
			if (toks.length >= 2 && toks[0].equals("@U")) {
				this.unit = toks[1];
				this.period = this.unit.length();
			}
			
			if (toks.length >= 2 && toks[0].equals("@N")) {
				this.unitNo = Double.parseDouble(toks[1]);
			}
			
			if (toks.length >= 2 && toks[0].equals("@S")) {
				setScore(Double.parseDouble(toks[1]));
			}
			if (toks.length >= 2 && toks[0].equals("@E")) {
				entropy = Double.parseDouble(toks[1]);
			}
			if (toks.length >= 2 && toks[0].equals("@M")) {
				this.pMatch = Double.parseDouble(toks[1]);
			}
			if (toks.length >= 2 && toks[0].equals("@I")) {
				this.pIndel = Double.parseDouble(toks[1]);
			}
			
		}
	}
	/**
	 * Read a repeat from using a xaf reader
	 * @param reader
	 * @return
	 */
	public static TandemRepeat read(XAFReader reader){
		TandemRepeat rec = new TandemRepeat();
		ArrayList<?> headerList = reader.getHeaderList();
		for (int i = 0; i< headerList.size(); i++){
			String fieldStr = headerList.get(i).toString();
			if (TandemRepeat.chromHd.equals(fieldStr))
				rec.setChr(reader.getField(i));
			else if (TandemRepeat.idHd.equals(fieldStr))
				rec.setID(reader.getField(i));			
			else if (TandemRepeat.startHd.equals(fieldStr))
				rec.setStart(Integer.parseInt(reader.getField(i)));
			else if (TandemRepeat.endHd.equals(fieldStr))
				rec.setEnd(Integer.parseInt(reader.getField(i)));
			else if (TandemRepeat.periodHd.equals(fieldStr))
				rec.setPeriod(Integer.parseInt(reader.getField(i)));			
			else if (TandemRepeat.unitHd.equals(fieldStr))
				rec.unit = reader.getField(i);
			else if (TandemRepeat.unitNoHd.equals(fieldStr))
				rec.setUnitNo(Double.parseDouble(reader.getField(i)));
			else if (TandemRepeat.scoreHd.equals(fieldStr))
				rec.setScore((Double.parseDouble(reader.getField(i))));
			else if (TandemRepeat.entropyHd.equals(fieldStr))
				rec.entropy = (Double.parseDouble(reader.getField(i)));
			else if (TandemRepeat.pMatchHd.equals(fieldStr))
				rec.pMatch = (Double.parseDouble(reader.getField(i)));
			else if (TandemRepeat.pIndelHd.equals(fieldStr))
				rec.pIndel = (Double.parseDouble(reader.getField(i)));
			else if (TandemRepeat.annotationHd.equals(fieldStr))
				rec.annotations = reader.getField(i);
			//TODO: To be removed
			else if (TandemRepeat.chrHd.equals(fieldStr))
				rec.setChr(reader.getField(i));
		}
		return rec;
	}
	
	/**
	 * Parse the short tandem repeat from a line with a given header
	 * This method has been replaced by a xaf reader
	 * @param line
	 * @param hds
	 * @return
	 */
	@Deprecated
	public static TandemRepeat read(String line, String [] hds){		
		TandemRepeat rec = new TandemRepeat();
		String [] toks = line.trim().split("\\t");

		for (int i = 0; i < hds.length; i++ ){
			if (TandemRepeat.chrHd.equals(hds[i]))
				rec.setChr(toks[i]);
			if (TandemRepeat.chromHd.equals(hds[i]))
				rec.setChr(toks[i]);
			else if (TandemRepeat.idHd.equals(hds[i]))
				rec.setID(toks[i]);			
			else if (TandemRepeat.startHd.equals(hds[i]))
				rec.setStart(Integer.parseInt(toks[i]));
			else if (TandemRepeat.endHd.equals(hds[i]))
				rec.setEnd(Integer.parseInt(toks[i]));
			else if (TandemRepeat.periodHd.equals(hds[i]))
				rec.setPeriod(Integer.parseInt(toks[i]));			
			else if (TandemRepeat.unitHd.equals(hds[i]))
				rec.unit = toks[i];
			else if (TandemRepeat.unitNoHd.equals(hds[i]))
				rec.setUnitNo(Double.parseDouble(toks[i]));
			else if (TandemRepeat.scoreHd.equals(hds[i]))
				rec.setScore((Double.parseDouble(toks[i])));
			else if (TandemRepeat.entropyHd.equals(hds[i]))
				rec.entropy = (Double.parseDouble(toks[i]));
			else if (TandemRepeat.pMatchHd.equals(hds[i]))
				rec.pMatch = (Double.parseDouble(toks[i]));
			else if (TandemRepeat.pIndelHd.equals(hds[i]))
				rec.pIndel = (Double.parseDouble(toks[i]));
			else if (TandemRepeat.annotationHd.equals(hds[i]))
				rec.annotations = toks[i];
		}		
		return rec;
	}
	
	public String toString(String [] hds, boolean collapse){		
		StringBuffer sb = new StringBuffer();

		for (int i = 0; i < hds.length; i++ ){
			//if (TandemRepeat.chrHd.equals(hds[i]))
			//	sb.append(getChr()+"\t");
			if (TandemRepeat.chromHd.equals(hds[i]))
				sb.append(getChr()+"\t");
			else if (TandemRepeat.idHd.equals(hds[i]))
				sb.append(getID()+"\t");			
			else if (TandemRepeat.startHd.equals(hds[i]))
				sb.append(getStart()+"\t");				
			else if (TandemRepeat.endHd.equals(hds[i]))
				sb.append(getEnd()+"\t");				
			else if (TandemRepeat.periodHd.equals(hds[i]))
				sb.append(period+"\t");							
			else if (TandemRepeat.unitHd.equals(hds[i])){
				if (collapse)
					sb.append(collapseForm(unit)+"\t");
				else
					sb.append(unit+"\t");				
			}else if (TandemRepeat.unitNoHd.equals(hds[i]))
				sb.append(unitNo+"\t");				
			else if (TandemRepeat.scoreHd.equals(hds[i]))
				sb.append(getScore()+"\t");				
			else if (TandemRepeat.entropyHd.equals(hds[i]))
				sb.append(entropy+"\t");				
			else if (TandemRepeat.pMatchHd.equals(hds[i]))
				sb.append(pMatch+"\t");				
			else if (TandemRepeat.pIndelHd.equals(hds[i]))
				sb.append(pIndel+"\t");						
			else if (TandemRepeat.annotationHd.equals(hds[i]))
				sb.append(annotations+"\t");						
			
		}	
		return sb.toString().trim();
	}
	
	/**
	 * Read a list of short tandem repeats from a TR file
	 * This function will be soon replaced by xaf reader
	 * @param in
	 * @param desc
	 * @return
	 * @throws IOException
	 */
	@Deprecated
	public static ArrayList<TandemRepeat> readFromFile(BufferedReader in, ArrayList<String> desc) throws IOException{
		
		//Start with the default header
		String[] headers = STANDARD_HEADER;
		ArrayList<TandemRepeat> trfList = new ArrayList<TandemRepeat>();
		
		String line = "";
		while ((line = in.readLine()) != null){
			line = line.trim();
			if (line.length() == 0) continue;
			
			if (line.startsWith("#H:"))
				headers = line.substring(3).split("\\t");
			else if (line.startsWith("#")){
				if (desc != null)
					desc.add(line);
				continue;
			}else
				trfList.add(read(line, headers));
						
		}//while

		return trfList;
	}	
	


	
	/**
	 * @return the period
	 */
	public int getPeriod() {
		return period;
	}

	/**
	 * @param period the period to set
	 */
	public void setPeriod(int period) {
		this.period = period;
	}

	/**
	 * @return the chr
	 */
	public String getChr() {
		return getParent();
	}
	/**
	 * @param chr the chr to set
	 */
	public void setChr(String chr) {
		setParent(chr);
	}
	

	/**
	 * @return the unitNo
	 */
	public double getUnitNo() {
		return unitNo;
	}

	/**
	 * @param unitNo the unitNo to set
	 */
	public void setUnitNo(double unitNo) {
		this.unitNo = unitNo;
	}

	

	/**
	 * @return the unit
	 */
	public String getUnit() {
		return unit;
	}

	/**
	 * @param unit the unit to set
	 */
	public void setUnit(String unit) {
		this.unit = unit;
	}

	
	
	static public void writeToFile(SequenceOutputStream out, String[] headers, ArrayList<TandemRepeat> trs, ArrayList<String> dList, boolean collapseForm)
	throws IOException{
		if (dList != null){
			for (int index = 0; index < dList.size(); index++) {
				out.print(dList.get(index));
				out.print('\n');
			}
		}		
		for (int i =0; i < headers.length; i++){
			if (i == 0){
				out.print("#H:");
			}else
				out.print('\t');
			
			out.print(headers[i]);
		}
		out.print('\n');
		
		for (int index = 0; index < trs.size(); index++) {
			out.print(trs.get(index).toString(headers,collapseForm));
			out.print('\n');
		}
		
	}
	
	/**
	 * Find the most compact form of japsa.seq (triplet).
	 * This is a slow implementation.
	 * @param japsa.seq
	 * @return
	 */
	public static String collapseForm(String seq){
		seq = seq.toUpperCase();
		String ret = seq;		
		for (int i = 0; i < ret.length(); i++){
			seq = seq.substring(1) + seq.charAt(0);
			if (ret.compareTo(seq) > 0)
				ret = seq;			
		}
		//complement japsa.seq
		char [] cs = new char[seq.length()];
		for (int i = 0; i <seq.length();i++){
			char c = seq.charAt(i);
			if (c == 'A')			
				cs[cs.length - 1 - i] = 'T';
			else if (c == 'C')
				cs[cs.length - 1 - i] = 'G';
			else if (c == 'G')
				cs[cs.length - 1 - i] = 'C';
			else if (c == 'T')
				cs[cs.length - 1 - i] = 'A';
			else
				cs[cs.length - 1 - i] = 'A';			
		}
		
		seq = new String(cs);
		
		for (int i = 0; i < ret.length(); i++){
			seq = seq.substring(1) + seq.charAt(0);
			if (ret.compareTo(seq) > 0)
				ret = seq;			
		}		
		return ret;
	}
	
	
	/**
	 * Write the feature out in BED format
	 * Format: chr <startChr> <endChr> ID score
	 * 
	 */
	public void writeBED(SequenceOutputStream out) throws IOException{
		out.print((this.getParent()+'\t' + (getStart()-1) + '\t' + getEnd()+
				'\t' + unit + '\t' + getScore() + '\t' + (getStrand() == '-'?'-':'+') + '\t' + unitNo +'\n'));
	}	
}

