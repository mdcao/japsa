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
 * 08/01/2012 - Minh Duc Cao: Revised                                        
 * 28/03/2012 - Minh Duc Cao: Revised and change file format 
 * 14/11/2013 - Minh Duc Cao: revised
 ****************************************************************************/

package japsa.seq;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * A feature annotated on a sequence.
 * A feature is specified by its location (start and end).
 * Note that the location (start and end) is specified as 1-index as oppose to
 * the internal presentation of the sequence.  
 */

public class JapsaFeature implements Comparable<JapsaFeature> {	
	protected ArrayList<String> featureDesc;// Description of the feature
	/**
	 * The start point of feature (inclusive)
	 */
	private int start;
	
	/**
	 * The end point of feature (inclusive)
	 */	
	private int end;
	
	/**
	 * Type of the feature
	 */
	private String type;// Type of the feature eg. CDS, Gene, Exon etc
	
	/**
	 * ID
	 */
	private String id;// ID of the feature
	
	
	private char strand = ' ';
	private String parent = "";
	private double score;
	
	protected ArrayList<JapsaFeature> children = null;
	
	public ArrayList<JapsaFeature> getChildren(){
		return children;
	}
	
	public void setChildren(ArrayList<JapsaFeature> cr){
		children = cr;
	}

	/**
	 * Create a new feature
	 * @param desc
	 * @param start
	 * @param length
	 */
	public JapsaFeature(int start, int end, String type, String id,
			char strand, String parent) {
		super();
		this.start = start;
		this.end = end;
		this.type = type;
		this.id = id;
		featureDesc = new ArrayList<String>();
		this.strand = strand;
		this.parent = parent;
	}
	

	public JapsaFeature(int start, int end) {
		this(start, end, "", "", ' ', "");
	}	

	public JapsaFeature() {
		this(-1, 0);
	}

	/**
	 * @return the type
	 */
	public String getType() {
		return type;
	}
	
	/**
	 * @param type
	 * the type to set
	 */
	public void setType(String type) {
		this.type = type;
	}

	/**
	 * Get the ID of the feature
	 * @return
	 */
	public String getID() {
		return id;
	}

	public void setID(String id) {
		this.id = id;
	}

	public String getParent() {
		return parent;
	}

	public void setParent(String parent) {
		this.parent = parent;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	public void setScore(double score) {
		this.score = score;
	}
	
	/**
	 * @return the desc
	 */
	public String getDesc() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < featureDesc.size(); i++)
			sb.append(featureDesc.get(i) + "\n");

		return sb.toString();
	}

	/**
	 * @param desc
	 *            the desc to set
	 */
	public void addDesc(String desc) {
		featureDesc.add(desc);
	}

	public Iterator<String> getDescStr() {
		return featureDesc.iterator();
	}

	/**
	 * @return the length
	 */
	public int getLength() {
		return end - start + 1;
	}

	/**
	 * Return the start position in 1-index, inclusive of the feature 
	 * @return the start
	 */
	public int getStart() {
		return start;
	}

	/**
	 * Return the end position in 1-index, inclusive of the feature
	 * 
	 * 
	 * @return
	 */

	public int getEnd() {
		return end;
	}

	/**
	 * @param start
	 *            the start to set
	 */
	public void setStart(int start) {
		this.start = start;
	}
	
	public void setEnd(int end) {
		this.end = end;
	}

	
	public char getStrand() {
		return strand;
	}

	public void setStrand(char strand) {
		this.strand = strand;
	}
	
	public String toString() {
		return getParent()+"\t" + start + "\t" + end + "\t" + type + "\t" + getDesc();
	}

	/**
	 * Compare a feature to another one
	 */
	 public int compareTo(JapsaFeature feature) {
		if (start != feature.start)
			return (start - feature.start);
		else
			return  feature.end - end;		
		//The parent one should come before
	}
	 
	/**
	 * Check if this feature matches with another one. The two features are
	 * considered matched if they shared at least gap bases 
	 * @param feature
	 * @param gap
	 * @return
	 */
	public boolean match(JapsaFeature feature, int gap) {
		return (end - feature.start >= gap &&
		        feature.end - start >= gap);
	}
	

	/**
	 * Write the description of the feature to a output stream in BFF format
	 * @param out
	 * @throws IOException
	 */
	public void write(SequenceOutputStream out) throws IOException{		
		out.print(JapsaFileFormat.FEATURE_START + this.type + 
				  JapsaFileFormat.DELIMITER + this.start
				+ JapsaFileFormat.DELIMITER + this.end
				+ JapsaFileFormat.DELIMITER + this.id
				+ JapsaFileFormat.DELIMITER + this.strand
				+ JapsaFileFormat.DELIMITER + parent 
				+ JapsaFileFormat.DELIMITER + score  + '\n');
		
		
		
		for (int i = 0; i < featureDesc.size(); i++) {
			out.print(JapsaFileFormat.COMMENT);
			out.print(featureDesc.get(i));
			out.print('\n');			
		}
	}
	
	/**
	 * Write the feature out in BED format
	 * Format: chr <startChr> <endChr> ID score
	 * 
	 */
	public void writeBED(SequenceOutputStream out) throws IOException{
		out.print((this.parent+'\t' + (start-1) + "\t" + end+
				'\t' + id + '\t' + score + '\t' + (strand == '-'?'-':'+') + '\n'));
	}
	/**
	 * Return a list of features from a bed file
	 * @param fileName
	 * @return
	 * @throws IOException
	 */
	public static ArrayList<JapsaFeature> readBED(String fileName) throws IOException{
		ArrayList<JapsaFeature>  ret = new ArrayList<JapsaFeature> ();
		BufferedReader in = SequenceReader.openFile(fileName);
		String line = "";
		
		while ( (line = in.readLine()) != null){
			String [] toks = line.trim().split("\t");
			
			//Set the minimal information
			JapsaFeature feature = new JapsaFeature(1+Integer.parseInt(toks[1]), Integer.parseInt(toks[2]));			
			feature.setParent(toks[0]);
			
			if(toks.length > 3)
				feature.setID(toks[3]);
			
			if(toks.length > 4)
				feature.setScore(Double.parseDouble(toks[4]));
			
			if(toks.length > 5)
				feature.setStrand(toks[5].charAt(0));
			
			
			ret.add(feature);
		}//while		
		
		in.close();
		
		return ret;
	}

	/**
	 * Parse the header line for the feature
	 * @param line
	 * @return
	 */
	public static JapsaFeature readLine(String line) {
		// >>Gene:start:end:ID:isComplement( |+|-):parent (| actually :):score
		JapsaFeature currentFeature = new JapsaFeature();
		String[] tokens = line.substring(2).split(
				JapsaFileFormat.DELIMITER + "", 7);

		if (tokens.length > 0)
			currentFeature.setType(tokens[0]);
		
		if (tokens.length > 1)
			currentFeature.setStart(Integer.parseInt(tokens[1]));
		
		if (tokens.length > 2)
			currentFeature.setEnd(Integer.parseInt(tokens[2]));
		
		if (tokens.length > 3)
			currentFeature.setID(tokens[3]);
		
		if (tokens.length > 4 && tokens[4].length() > 0)
			currentFeature.strand = tokens[4].charAt(0);
		else
			currentFeature.strand = ' ';

		if (tokens.length > 5)
			currentFeature.parent = tokens[5];
		else
			currentFeature.parent = "";
		
		if (tokens.length > 6)
			currentFeature.score = Double.parseDouble(tokens[6]);
		else
			currentFeature.score= 0;


		return currentFeature;
	}

	
	/**
	 * Write the end of the feature to a stream
	 * @param out
	 * @throws IOException
	 */
	public void writeEnd(SequenceOutputStream out) throws IOException{
		out.print(JapsaFileFormat.FEATURE_END + this.type + 
				  JapsaFileFormat.DELIMITER + this.start
				+ JapsaFileFormat.DELIMITER + this.end
				+ JapsaFileFormat.DELIMITER + this.id
				+ JapsaFileFormat.DELIMITER + this.strand
				+ JapsaFileFormat.DELIMITER + parent
				+ JapsaFileFormat.DELIMITER + score + '\n');
	}
	
	
	public JapsaFeature cloneFeature() {
		JapsaFeature out = new JapsaFeature(start, end,
				type, id, strand,
				this.parent);
		
		for (int i = 0; i < featureDesc.size(); i++)
					out.addDesc(featureDesc.get(i));//

		return out;
	}
	
	
	public String getProperty() {
		String out = "Type = " + type;
		if (id.length() > 0)
			out = out + "; ID = " + id;
		out = out + "; Start = " + start + "; End = " + end;// +
																	// + parent;
		if (strand == '+')
			out += "\nDirection = Forward";
		else if (strand == '-')
			out += "\nDirection = Reverse";

		for (int i = 0; i < featureDesc.size(); i++) {
			out = out + "\n" + featureDesc.get(i);
		}
		return out;
	}


}