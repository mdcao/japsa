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
 *  
 ****************************************************************************/

package japsa.seq;

import japsa.util.MyBitSet;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;





/**
 * @author Minh Duc Cao
 * 
 */
public class JapsaAnnotation {

	private Sequence seq;//sequence annotated for
	private ArrayList<String> annoDescription;// Description of the feature	
	private ArrayList<JapsaFeature> featureList;
	private String annotationID;

	/**
	 * Create an annotation for a sequence 
	 * @param japsa.seq
	 */
	public JapsaAnnotation(Sequence seq) {
		this.seq = seq;
		featureList = new ArrayList<JapsaFeature>();
		annoDescription = new ArrayList<String>();
		if (seq != null)
			annotationID = seq.getName();
		else 
			annotationID = "";
	}


	/**
	 * This constructor is equivalent to japsa.seq = null
	 */	

	public JapsaAnnotation() {
		this.seq = null;
		featureList = new ArrayList<JapsaFeature>();
		annoDescription = new ArrayList<String>();
		annotationID = "";
	}

	/**
	 * Construct an annotation with ID for sequence japsa.seq.
	 * @param japsa.seq
	 * @param ID
	 */
	public JapsaAnnotation(Sequence seq, String ID) {
		this(seq);
		annotationID = ID;
	}



	public void sortFeatures() {
		Collections.sort(featureList);
	}

	/**
	 * Clone from an existing one
	 * 
	 * @param anAnno
	 */
	public JapsaAnnotation(JapsaAnnotation anAnno) {
		this(anAnno.seq);

		annotationID = "" + anAnno.annotationID;

		for (int i = 0; i < anAnno.annoDescription.size(); i++) {
			this.addDescription(anAnno.annoDescription.get(i));
		}

		for (int i = 0; i < anAnno.featureList.size(); i++) {
			this.featureList.add(anAnno.featureList.get(i));
		}
	}

	/**
	 * @return the desc
	 */
	public String getDescription() {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < annoDescription.size(); i++)
			sb.append(annoDescription.get(i) + "\n");

		return sb.toString();
	}

	/**
	 * If the annotation does not contain any annotations
	 * @return
	 */
	public boolean isEmpty() {
		return featureList.isEmpty();
	}

	/**
	 * @param desc
	 *            the desc to set
	 */
	public void addDescription(String desc) {
		annoDescription.add(desc);
	}

	public String getAnnotationID() {
		return annotationID;
	}

	public void setAnnotationID(String annotationID) {
		//replace all the deliminaters
		this.annotationID = annotationID.replace(JapsaFileFormat.DELIMITER,'_');
	}

	/**
	 * Add a feature to the annotation
	 * @param feature
	 */
	public void add(JapsaFeature feature) {
		featureList.add(feature);
	}

	/**
	 * Return the feature at index ind 
	 * @param ind
	 * @return
	 */
	public JapsaFeature getFeature(int ind) {
		return featureList.get(ind);
	}


	public void setSequence(Sequence seq){
		this.seq = seq;
	}

	public Sequence getSequence(){
		return seq;
	}

	/**
	 * Return the number of feature
	 * @return
	 */

	public int numFeatures() {
		return featureList.size();
	}


	/**
	 * Remove a feature
	 * 
	 * @param feature
	 * @return
	 */
	public boolean remove(JapsaFeature feature) {
		return featureList.remove(feature);
	}

	/**
	 * Get the feature interator
	 * @return
	 */
	public Iterator<JapsaFeature> iterator() {
		return featureList.iterator();
	}


	/**
	 * Get the list of reature as an ArrayList
	 * @return
	 */
	public ArrayList<JapsaFeature>  getFeatureList() {
		return featureList;
	}

	public static JapsaAnnotation readDataFromFile(String fileName)
			throws IOException {

		JapsaFileFormat reader = new JapsaFileFormat(fileName);
		JapsaAnnotation anno = reader.readAnnotation();
		reader.close();

		return anno;
	}

	public String toString() {
		return annotationID;
	}



	/**
	 * Compare the annotation against a model annotation.
	 * Two features from two annotations are considered to be matched if they
	 * overlap with each other for at least gap bases
	 * @param model
	 */

	public void compareAnnotationAtFeatureLevel(JapsaAnnotation model, int gap) {
		// A slow implementation
		int TP = 0, FP = 0, FN = 0;

		//Check the positives
		for (int i = 0; i < this.numFeatures(); i++) {
			JapsaFeature myFeature = this.getFeature(i);
			boolean match = false;
			for (int j = 0; j < model.numFeatures(); j++) {
				if (model.getFeature(j).match(myFeature, gap)){
					match = true;
					break;
				}//if
			}//for
			if (match){
				TP ++;
			}else
				FP ++;
		}

		double precision = TP * 1.0 / (TP + FP);

		TP = 0;
		FN = 0;
		for (int i = 0; i < model.numFeatures(); i++) {
			JapsaFeature modelFeature = model.getFeature(i);
			boolean match = false;
			for (int j = 0; j < this.numFeatures(); j++) {
				if (this.getFeature(j).match(modelFeature, gap)){
					match = true;
					break;
				}//if
			}//for
			if (match){
				TP ++;
			}else
				FN ++;
		}

		double recall = TP * 1.0 / (TP + FN);
		double Fscore = 2.0 * (precision * recall) / (precision + recall);

		System.out.println("#=#Recall = " + recall);
		System.out.println("#=#Precison = " + precision);
		System.out.println("F-score = " + Fscore);
	}


	/**
	 * Compare this annotation to a model annotation Count the number of
	 * covered, not covered PreCond: Mine is sorted, so is model
	 * 
	 * @param length
	 * @param model
	 */
	public void compareAnnotationAtBaseLevel(JapsaAnnotation model) {
		if (model.seq == null && model.seq == null ){
			throw new RuntimeException("Unable to compare since the total number of bases cannot be determined");
		}

		int myLength = -1;

		if (model.seq == null){
			myLength = this.seq.length();
		}else if (this.seq == null)
			myLength = model.seq.length();
		else if (seq != model.seq)
			throw new RuntimeException("Two annotations are for two different sequence");
		else
			myLength = this.seq.length();


		MyBitSet myBitSet = new MyBitSet(myLength),
				modelBitSet = new MyBitSet(myLength);



		//set the bases covered by the model annotation
		for (int x = 0; x < this.numFeatures(); x++) {
			for (int y = this.getFeature(x).getStart() - 1; // The real coordinator
					y < this.getFeature(x).getEnd(); y++) {
				myBitSet.set(y);
			}
		}

		//set the bases covered by the this annotation annotation
		for (int x = 0; x < model.numFeatures(); x++) {
			for (int y = model.getFeature(x).getStart() - 1; y < model.getFeature(x).getEnd(); y++) {
				modelBitSet.set(y);
			}
		}

		// Compare with the model annotation
		int FP = 0, FN = 0;
		int TP = 0, TN = 0;


		for (int i = 0; i < myLength; i++) {
			if (myBitSet.get(i))
				if (modelBitSet.get(i))
					TP++;
				else
					FP++;
			else if (modelBitSet.get(i))
				FN++;
			else
				TN++;
		}

		System.out.println("in This = " + FP + "(FP)");
		System.out.println("in Model = " + FN + " FN ");
		System.out.println("in None = " + TN + " (TN) ");
		System.out.println("in Both = " + TP + " (TP) ");

		double precision   = TP * 1.0 / (TP + FP), 
				recall      = TP * 1.0 / (TP + FN), //aka sensitivity or PPV
				specificity = TN * 1.0 / (TN + FP), //aka TPR
				accuracy    = (TP + TN + 0.0) / (TN + FP + TP + FN);
		double Fscore = 2.0 * (precision * recall) / (precision + recall);


		System.out.println("#=#Recall = " + recall);
		System.out.println("#=#Precison = " + precision);
		System.out.println("Specificity = " + specificity);
		System.out.println("Accuracy = " + accuracy);
		System.out.println("F-score = " + Fscore);


		System.out.printf("##%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
				model.annotationID, FP, FN, TN, TP, TP + FN, TP + FP, TP + FP
				+ TN + FN);
	}


	//Return the annotation of the complementary sequence
	public JapsaAnnotation reverse() {
		if (seq == null){
			throw new RuntimeException("The sequence must not be null");
		}
		JapsaAnnotation anno = new JapsaAnnotation(seq);
		int myLength = seq.length();

		for (int i = this.numFeatures() - 1; i >= 0; i--) {
			JapsaFeature thisFeature = this.getFeature(i);

			JapsaFeature feature = 
					new JapsaFeature(myLength - (thisFeature.getEnd() + 1),
							thisFeature.getEnd(), 
							thisFeature.getType(), thisFeature.getID(),
							thisFeature.getStrand(), thisFeature.getParent());

			for (int x = 0; x < thisFeature.featureDesc.size(); x++)
				feature.addDesc(thisFeature.featureDesc.get(x));

			if (thisFeature.getStrand() == '+')
				feature.setStrand('-');
			else if (thisFeature.getStrand() == '-')
				feature.setStrand('+');

			anno.add(feature);
		}

		return anno;
	}

	/**
	 * Combine this annotation with another annotation and return the combined
	 * annotation
	 * 
	 * @param rev
	 * @return
	 */
	public JapsaAnnotation combineWith(JapsaAnnotation another) {
		if (seq != another.seq){
			throw new RuntimeException("The two annotations must be for the same sequence");
		}

		JapsaAnnotation newAnno = new JapsaAnnotation(this);


		for (int i = 0; i < another.annoDescription.size(); i++) {
			this.addDescription(another.annoDescription.get(i));
		}

		for (int i = 0; i < another.featureList.size(); i++) {
			newAnno.featureList.add(another.featureList.get(i));
		}

		newAnno.sortFeatures();
		return newAnno;
	}


	/**
	 * Check if the annotation is sorted
	 * @return
	 */
	public boolean isSort() {
		for (int i = 1; i < this.numFeatures(); i++) {
			if (this.getFeature(i).compareTo(this.getFeature(i - 1)) < 0){
				return false;
			}
		}

		return true;
	}	



	/**
	 * Write the description of the annotation (without any feature) to the
	 * stream
	 * @param out
	 * @throws IOException
	 */
	public void writeDescription(SequenceOutputStream out) throws IOException{		
		for (int i = 0; i < annoDescription.size(); i++) {
			out.print(JapsaFileFormat.ANNOTATION_COMMENT);			
			out.print(annoDescription.get(i));
			out.print('\n');
		}
	}

	/**
	 * Write the annotation to a stream in JAPSA File Format 
	 * @param out
	 * @param needSort: need to sort before writing
	 * @throws IOException
	 */

	private void writeAnnotation(SequenceOutputStream out, boolean needSort) throws IOException{

		out.print(JapsaFileFormat.HEADER);
		out.print(JapsaFileFormat.DELIMITER);
		out.print(annotationID);
		out.print(JapsaFileFormat.DELIMITER);
		if (seq != null)
			out.print(seq.length());
		out.print(JapsaFileFormat.DELIMITER);
		if (seq != null)
			out.print(seq.alphabet().getName());

		out.print('\n');

		writeDescription(out);
		out.print('\n');

		if (needSort && !isSort())
			sortFeatures();

		Iterator<JapsaFeature> iter = featureList.iterator();

		while (iter.hasNext()) {
			iter.next().write(out);
			out.print('\n');
		}		
	}

	/**
	 * Write the annotation and sequence if available to a stream in JAPSA 
	 * @param out 
	 * @throws IOException
	 */	
	public void write(SequenceOutputStream out) throws IOException{
		//TODO: Need to automatically update needSort when a new feature is added
		boolean needSort = true;
		if (seq == null)
			writeAnnotation(out, needSort);
		else{
			if (needSort && !isSort())
				sortFeatures();
			write(seq, this, out);
		}
	}



	/**
	 * Write the annotation to a stream in JAPSA 
	 * @param out 
	 * @throws IOException
	 */

	public void writeAnnotation(SequenceOutputStream out) throws IOException{
		writeAnnotation(out,true);
	}


	/**
	 * Write annotation to a stream in bed format
	 * @param out
	 * @throws IOException
	 */

	public void writeBED(SequenceOutputStream out) throws IOException{
		if (!isSort())
			sortFeatures();
		Iterator<JapsaFeature> iter = featureList.iterator();

		while (iter.hasNext()) {
			iter.next().writeBED(out);			
		}
	}

	public void writeFeatureSequence(SequenceOutputStream out) throws IOException{
		if (seq == null){
			throw new RuntimeException("The sequence is not known");
		}

		Iterator<JapsaFeature> iter = featureList.iterator();
		while (iter.hasNext()) {
			JapsaFeature feature = iter.next();
			int start = feature.getStart();
			int end = feature.getEnd();

			Sequence fSeq = seq.subSequence(start - 1, end);
			fSeq.setName(feature.getID());

			fSeq.setDesc((feature.getParent() + " " + feature.getDesc()).replaceAll("\n", ";"));
			fSeq.writeFasta(out);			
		}
	}


	// All feature must be sorted in oder of start points
	/**
	 * Write sequence to a stream
	 */
	public static void write(Sequence sequence, JapsaAnnotation anno,
			SequenceOutputStream out) throws IOException{				
		if (sequence == null && anno == null)
			return;
		else if (sequence == null)
			anno.writeAnnotation(out);
		else if (anno == null)
			sequence.print(out);
		else {// Both are not null
			//char[] charSeq = sequence.charSequence();

			String seqID = anno.getAnnotationID();
			if (seqID == null || seqID.length() == 0)
				seqID = sequence.getName();

			//FIXME: Unify all the write header
			out.print(JapsaFileFormat.HEADER);
			out.print(JapsaFileFormat.DELIMITER);
			out.print(seqID);
			out.print(JapsaFileFormat.DELIMITER);
			out.print(sequence.length());
			out.print(JapsaFileFormat.DELIMITER);
			out.print(sequence.alphabet().getName());
			out.print('\n');

			sequence.writeDescription(out);
			anno.writeDescription(out);

			// An iterator of features
			Iterator<JapsaFeature> iter = anno.iterator();
			JapsaFeature nextFeature = null;

			LinkedList<JapsaFeature> fin = new LinkedList<JapsaFeature>();
			boolean writeFeature = false;

			JapsaFeature dummy = new JapsaFeature(-1, 0);
			nextFeature = dummy;
			int nextFeatureFinish = Integer.MAX_VALUE;

			if (iter.hasNext()) {
				nextFeature = iter.next();
			}

			// Start writing
			for (int x = 0; x < sequence.length(); x++) {
				// if finish of (a) features				
				while (x >= nextFeatureFinish - 1) {
					writeFeature = true;
					JapsaFeature finishedF = fin.remove();
					out.print('\n');					
					finishedF.writeEnd(out);

					if (fin.size() > 0) {
						nextFeatureFinish = fin.get(0).getStart()
								+ fin.get(0).getLength();
					} else
						nextFeatureFinish = Integer.MAX_VALUE;
				}// while
				// If start of features
				while (x == nextFeature.getStart() - 1) {// The java cordinate
					// off to print feature
					writeFeature = true;
					out.print('\n');
					out.print('\n');					
					nextFeature.write(out);

					// Add this feature to the fin list only if its length > 0
					if (nextFeature.getLength() > 0) {
						int finishPoint = nextFeature.getStart()
								+ nextFeature.getLength();
						Iterator<JapsaFeature> iterFin = fin.listIterator();
						int indexToAdd = 0;
						while (iterFin.hasNext()) {
							JapsaFeature nextF = iterFin.next();
							if (nextF.getStart() + nextF.getLength() >= finishPoint) {
								break;
							}
							indexToAdd++;
						}
						fin.add(indexToAdd, nextFeature);
						// Change the next finish
						nextFeatureFinish = fin.get(0).getStart()
								+ fin.get(0).getLength();
					}

					// Move to the next feature
					nextFeature = dummy;
					if (iter.hasNext()) {
						nextFeature = iter.next();
					}
				}// while				
				if (writeFeature) {
					writeFeature = false;
					int offset = x % JapsaFileFormat.CHAR_PER_LINE;
					if (offset != 0) {				
						out.print(x+1,10);
						out.print(' ');
						out.print(' ');

						for (int ii = 0; ii < offset; ii++) {
							if (ii % JapsaFileFormat.CHAR_PER_BLOCK == 0)
								out.print(' ');
							out.print(' ');
						}
					}
				}
				// write the symbol
				if (x % JapsaFileFormat.CHAR_PER_LINE == 0){					
					out.print('\n');
					out.print(x+1,10);
					out.print(' ');
					out.print(' ');			
				}

				if (x % JapsaFileFormat.CHAR_PER_BLOCK == 0)
					out.print(' ');

				out.print(sequence.charAt(x));
				//charSeq[x]);
			}
			out.print('\n');			
			//out.flush();			
		}
	}
	

	/**
	 * Read annotation from a gff file and store in a hash table index by sequence ID
	 * @param in
	 * @param upStr
	 * @param downStr
	 * @param list: list of feature to read, ie "CDS,mRNA". "all" means every thing 
	 * @return
	 * @throws IOException
	 */
	public static ArrayList<JapsaAnnotation> readMGFF(InputStream in, int upStr, int downStr, String list) throws IOException{
		boolean notAll = !list.equals("all");
		
		
		ArrayList<JapsaAnnotation> annoList = new ArrayList<JapsaAnnotation>();		
		HashMap<String, JapsaAnnotation> annoMap = new  HashMap<String, JapsaAnnotation>();
		
		String line;		
		FastaReader reader = new FastaReader(in);
		reader.nextLine();
		
		JapsaAnnotation currentAnno = null;
		while ( reader.nextLine() > 0){
			//lineNo ++;
			line = reader.getCurrentLine();
			//line = line.trim();
			if (line.startsWith("##sequence-region")){
				JapsaAnnotation anno =  new JapsaAnnotation();
				anno.setAnnotationID(line.split(" ")[1]);				
				annoMap.put(anno.annotationID, anno);				
				annoList.add(anno);
				
				if (currentAnno == null)
					currentAnno = anno;
				continue;
			}
			if (line.startsWith("##FASTA")){
				break;//pass on to fasta reader
			}
			if (line.startsWith("#")){
				continue;//unknown lines
			}
			
			//assert currentAnno != null
			String [] toks = line.split("\t");
			//if (toks[2].equals("region"))
			//	continue;

			if (toks.length < 2)
				continue;

			String type = toks[2];
			if (notAll && !list.contains(type))
				continue;

			int start = Integer.parseInt(toks[3]);
			int end = Integer.parseInt(toks[4]);

			String parent = toks[0];
			if (!parent.equals(currentAnno.annotationID)){
				currentAnno = annoMap.get(parent);
			}
			if (currentAnno == null){
				in.close();
				reader.close();
				throw new RuntimeException("Sequence region ID " + parent + " not found");
			}			
			
			String ID = "";
			char strand = toks[6].charAt(0);
			String desc = toks[8];

			String [] featureChars = toks[8].split(";");
			for (String fch:featureChars){
				if (fch.startsWith("ID="))
					ID = fch.substring(3);

				if (fch.startsWith("Parent="))
					parent = fch.substring(7);				
			}

			JapsaFeature feature = new JapsaFeature(start, end, type, ID, strand, parent);
			feature.addDesc(desc);			
			currentAnno.add(feature);

			//Add upstream
			if (upStr > 0){// && feature.getType().equals("gene")){
				if (feature.getStrand() == '-'){
					JapsaFeature uFeature = new JapsaFeature(end + 1, end + upStr, "upstream", "u" + ID, strand, feature.getID());
					uFeature.addDesc("Upstream region added automatically");
					currentAnno.add(uFeature);
				}else{
					JapsaFeature uFeature = new JapsaFeature(start - upStr, start - 1, "upstream", "u" + ID, strand, feature.getID());
					uFeature.addDesc("Upstream region added automatically");
					currentAnno.add(uFeature);
				}
			}

			//add downstream
			if (downStr > 0){// && feature.getType().equals("gene")){
				if (feature.getStrand() == '-'){
					JapsaFeature dFeature = new JapsaFeature(start - downStr, start - 1, "downtream", "d" + ID, strand, feature.getID());
					dFeature.addDesc("Downstream region added automatically");
					currentAnno.add(dFeature);					
				}else{
					JapsaFeature dFeature = new JapsaFeature(end + 1, end + downStr, "downtream", "d" + ID, strand, feature.getID());
					dFeature.addDesc("Downstream region added automatically");
					currentAnno.add(dFeature);					
				}//else
			}//if downStr
		}		
		//FastaReader reader = new FastaReader(in.);
		reader.nextByte();
		Sequence seq;
		while ((seq = reader.nextSequence(Alphabet.DNA()))!=null){
			JapsaAnnotation anno = annoMap.get(seq.getName());
			anno.setSequence(seq);	
		}
		reader.close();
		
		return annoList;
	}	
}
