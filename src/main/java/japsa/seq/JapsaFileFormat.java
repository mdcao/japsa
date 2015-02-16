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
 * 14/11/2013 - Minh Duc Cao: Revised                                        
 *  
 ****************************************************************************/

package japsa.seq;

import japsa.seq.Alphabet;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;


/**
 * This is a improved version of the previous BioComp file format. This class
 * defines JAPSA file format as well as acts like a reader for the format
 * 
 * Format: A sequence start with #JSA
 * -The header  #JSA:<seqID>[:<length>[:<alphabetType]]
 * -Comments start with # (not followd by JSA)
 * -Feature start with ">>", follows by <startPos>:<endPos>
 * -Feature header >>ID:startPos:endPos:ID:strand:parentID
 * -Annotation for the feature are in the following lines, starting with
 * one '>'
 * -Feature ends with '<' (optional)
 * 	  
 *  Features can be interleaved with sequence data
 *  
 * @author Minh Duc Cao
 *
 */

public class JapsaFileFormat extends SequenceReader{
	/**
	 * The header
	 */
	public static final String HEADER = "#JSA";
	public static final char DELIMITER = ':';
	
	public static final int CHAR_PER_LINE = 50;
	public static final int CHAR_PER_BLOCK = 10;
	
	public static final String FEATURE_START=">>";
	public static final String FEATURE_END="<";
	public static final String COMMENT=">";
	public static final String ANNOTATION_COMMENT=COMMENT+"A:";
	public static final String SEQUENCE_COMMENT=COMMENT+"S:";
	
	
	private int firstByte = 0;
	
	
	
	private byte [] seq = new byte[1 << 20];
	private int seqIndex = 0;
	
	
	public JapsaFileFormat(InputStream sin) throws IOException{
		super(sin);		
		//Read next line
		nextLine();
	}
	
	public JapsaFileFormat(String fileName) throws IOException{
		super(fileName);
		
		//Read next line
		nextLine();
	}
	
	
	private void stripSpace(){
		firstByte = 0;
		while ((nextLine[firstByte] == 32 || nextLine[firstByte] == 9) && firstByte < nextLineLength)
			firstByte ++;
	}
	
	/* (non-Javadoc)
	 * @see japsa.seq.SequenceReader#nextSequence(japsa.seq.Alphabet)
	 */
	@Override
	public Sequence nextSequence(Alphabet alphabet) throws IOException {
		JapsaAnnotation anno = readAnnotation();
		if(anno == null)
			return null;
		else 
			return anno.getSequence();
	}
	
	
	
	public JapsaAnnotation readAnnotation() throws IOException{
		return readAnnotation(Alphabet.DNA16());
	}
	
	private JapsaAnnotation readAnnotation(Alphabet alphabet) throws IOException{
		if(eof)
			return null;		
		
		JapsaFeature currentFeature = null;
		
		stripSpace();		
		String line = new String(nextLine, firstByte, nextLineLength - firstByte);		
	//	int mode = -1;// start = -1, anno desc = 1,featureDesc = 2		
		//mode = 0: header line
		//mode = 1: comment line ;
		//mode = 2: sequence
		//mode = 3: >>
		//mode = 4: <
		//mode = 5: >
		
		if (!line.startsWith(HEADER)){
			throw new RuntimeException("'"+HEADER+"' expected at line " + lineNo);
		}
		
		StringBuilder seqDesc = new StringBuilder();
		
		
		String[] tokens = line.split(DELIMITER + "", 4);
		JapsaAnnotation anno = new JapsaAnnotation();
		
		seqIndex = 0;
		
		if (tokens.length > 1) {
			anno.setAnnotationID(tokens[1]);
		}
		
		if (tokens.length > 3) {
			if (tokens[3].toUpperCase().contains("PROTEIN"))
				alphabet = Alphabet.PROTEIN();
			else if (tokens[3].toUpperCase().contains("RNA"))
				alphabet = Alphabet.RNA();
			else if (tokens[3].toUpperCase().equals("DNA16"))
				alphabet = Alphabet.DNA16();
			else if (tokens[3].toUpperCase().equals("DNA5"))
				alphabet = Alphabet.DNA5();
		}
		
		while (!eof){
			nextLine();
			stripSpace();		
			line = new String(nextLine, firstByte, nextLineLength - firstByte);
		//	mode = 0;// Desc for annotation
			
			line = line.trim();
			if (line.length() <= 0)
				continue;

			if (line.startsWith(HEADER))
				break;
			
			//Ignore comments			
			if (line.startsWith("#"))
				continue;
			
			if (line.startsWith(FEATURE_START)) {// start a new feature
				currentFeature = JapsaFeature.readLine(line);
				anno.add(currentFeature);
			//	mode = 1;
			}else if (line.startsWith(COMMENT)) {
				if (currentFeature == null){
					if (line.startsWith(ANNOTATION_COMMENT))
						anno.addDescription(line.substring(3));
					else if (line.startsWith(SEQUENCE_COMMENT))
						//FIXME: Add this description to sequence
						//anno.addDescription(line.substring(3));
						seqDesc.append(line.substring(3));					
					else 
						anno.addDescription(line.substring(1));
					
				}else
					currentFeature.addDesc(line.substring(1));// feature
			}else if (line.startsWith(FEATURE_END))
				currentFeature = null;
			else{//sequence mode
				for (int i = firstByte; i < nextLineLength; i ++){				
					int nucleotide = alphabet.byte2index(nextLine[i]);
					if (nucleotide >= 0){//valid nucleotide				
						//ensure the sequence array is big enough
						if (seqIndex >= seq.length) {// Full
							int newLength = seq.length * 2;
							if (newLength < 0) {
								newLength = Integer.MAX_VALUE;
							}
							// if the array is extended
							if (newLength <= seqIndex) {
								// in.close();
								throw new RuntimeException(
										"Sequence is too long to handle at line " + lineNo + "  " + newLength);
							}
							// create new array of byte
							seq = Arrays.copyOf(seq, newLength);
						}
						//next symbol
						seq[seqIndex++] = (byte) nucleotide;
						//seqIndex++;			
					}//if nucleotide >= 0
				}//for i
			}//else			
		}//while

		if (seqIndex > 0){
			Sequence sequence = new Sequence(alphabet, seq, seqIndex, anno.getAnnotationID());
			sequence.setDesc(seqDesc.toString());
			anno.setSequence(sequence);
		}		
		
		return anno;
	}	
}