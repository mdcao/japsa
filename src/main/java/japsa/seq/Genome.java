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
 * 16/12/2012 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.seq;

import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * A genome is a set of sequences, possibly chromosomes, as well as other
 * sequences such as mitochrondria, contigs, unassembled sequences etc.
 * 
 * @author minhduc
 * 
 */
public class Genome {
	long length;
	
	ArrayList<Sequence> seqList = new ArrayList<Sequence>(); 
	HashMap<String, Sequence> seqHash = new HashMap<String, Sequence>();
	Alphabet.DNA alphabet = Alphabet.DNA();// the most permissive

	public Genome() {

	}
	
	/**
	 * @param args
	 */
	public void read(String fileName) throws IOException {
		// ArrayList <Sequence> seqHash = new ArrayList<Sequence>();
		FastaReader.Faster fread = new FastaReader.Faster(fileName);
		Sequence seq;

		// Timer timer = new Timer();
		while ((seq = fread.nextSequence(Alphabet.DNA16())) != null) {
			String id = seq.getName().split("\\s")[0];
			seq.setName(id);
			seqHash.put(id, seq);
			length += seq.length();
			seqList.add(seq);
		}
	}	
	
	public long getLength(){
		return length;
	}
	
	public ArrayList<Sequence> chrList(){
		return seqList;
	}


	public static void main(String[] args) throws IOException {
		Genome genome = new Genome();
		genome.read(args[0]);
		genome.getProteome(args[1], args[2]);
	}

	public void getProteome(String annotationFile, String canonicalFile)
			throws IOException {
		BufferedReader canBf = new BufferedReader(new FileReader(canonicalFile));
		BufferedReader annBf = new BufferedReader(
				new FileReader(annotationFile));

		String canLine = "";
		String annoLine = "";
		int annoLineNo = 0;

		Sequence seq = null;
		// int canonicalGene = 0;
		while ((canLine = canBf.readLine()) != null) {
			String[] canToks = canLine.split("\t");
			if (seq == null || (!seq.getName().equals(canToks[0]))) {
				seq = seqHash.get(canToks[0]);
			}
			if (seq == null)
				continue;
			// assert japsa.seq.id = canToks[0]
			// look for the annotation for this gene

			// canonicalGene ++;
			while ((annoLine = annBf.readLine()) != null) {
				annoLineNo++;
				String[] annToks = annoLine.split("\t");
				if (annToks[1].equals(canToks[0])
						&& annToks[0].equals(canToks[4])) {
					int exonNumber = Integer.parseInt(annToks[7]);
					String[] eStart = annToks[8].split(",");
					String[] eEnd = annToks[9].split(",");
					// boolean sense = "+".equals(annToks[7]);
					if (eStart.length != exonNumber) {
						canBf.close();
						annBf.close();
						throw new RuntimeException(annoLineNo + " : E = "
								+ exonNumber + "  " + eStart.length);
					}
					System.out.println("\n   " + exonNumber + ":");
					int myLength = 0;
					for (int i = 0; i < exonNumber; i++) {
						myLength += Integer.parseInt(eEnd[i])
								- Integer.parseInt(eStart[i]) + 1;
						System.out.print("   "
								+ (Integer.parseInt(eEnd[i])
										- Integer.parseInt(eStart[i]) + 1));
					}

					if (myLength % 3 != 0) {
						canBf.close();
						annBf.close();
						throw new RuntimeException(annoLineNo
								+ " : Not division by 3 " + length);
					}
					break;
				}
			}//
			if (annoLine == null) {
				canBf.close();
				annBf.close();
				throw new RuntimeException(annoLineNo + " : Null not expected "
						+ canToks[4]);
			}
		}
		canBf.close();
		annBf.close();
	}

}
