/*****************************************************************************
 * Copyright (c) 2010 Minh Duc Cao, Monash University.  All rights reserved. *
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
 * 3. Neither the name of Monash University nor the names of its contributors*
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

package japsadev.xm.genome;

import japsa.seq.JapsaFileFormat;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

import java.io.*;
import java.util.*;

public class GenomeSequence {
	public static final int ALPHABET_SIZE = 4;
	public static char[] NUCLEOTIDE = { 'a', 'c', 'g', 't', 'n' };

	int[] seqs;// use int to avoid the time of casting
	private long length;
	String seqID = "";
	String info = null;

	static int LAST2BIT_MASK = 0x3;
	static int LAST4BIT_MASK = 0xf;
	static int[] MASKS = new int[16], INV_MASKS = new int[16];
	static {
		MASKS[15] = 0x3;
		INV_MASKS[15] = ~MASKS[15];
		for (int i = MASKS.length - 2; i >= 0; i--) {
			MASKS[i] = MASKS[i + 1] << 2;
			INV_MASKS[i] = ~MASKS[i];
		}
	}

	public long getLength() {
		return length;
	}

	public static GenomeSequence concat(GenomeSequence[] gseqs) {
		long l = 0l;
		for (int i = 0; i < gseqs.length; i++)
			l += gseqs[i].length;

		GenomeSequence gSeq = new GenomeSequence(l);
		// TODO: this is a very inefficient way
		l = 0;
		for (int i = 0; i < gseqs.length; i++) {
			for (long y = 0; y < gseqs[i].length; y++) {
				gSeq.putBase(l, gseqs[i].getBase(y));
				l++;
			}
		}

		return gSeq;
	}

	public GenomeSequence(long l) {
		length = l;
		if (length > (1l << 34)) {
			System.out.println("Warning: length cant be larger than "
					+ (1l << 34));
			length = (1l << 34) - 1;
		}

		int size = (int) ((length - 1) >> 4) + 1;
		seqs = new int[size];
	}

	public GenomeSequence(char[] charSeq) {
		// TODO: this is a very inefficient way, but will optimise later
		// such as: put 16 chars in an int before put in the array
		this(charSeq.length);

		for (int i = 0; i < charSeq.length; i++) {
			putChar(i, charSeq[i]);
		}
	}

	public GenomeSequence(byte[] byteSeq) {
		// TODO: this is a very inefficient way, but will optimise later
		// such as: put 16 chars in an int before put in the array
		this(byteSeq.length);

		for (int i = 0; i < byteSeq.length; i++) {
			putBase(i, byteSeq[i]);
		}
	}

	public void increaseSize(long l) {
		if (l < length)
			return;
		length = l;
		if (length > (1l << 34)) {
			System.out.println("Warning: length cant be larger than "
					+ (1l << 34));
			length = (1l << 34) - 1;
		}
		int size = (int) ((length - 1) >> 4) + 1;
		int[] newSeqs = new int[size];
		// copy over
		for (int i = 0; i < seqs.length; i++)
			newSeqs[i] = seqs[i];

		seqs = newSeqs;
		// deallocate old seqHash
	}

	public static byte charToByte(char c) {
		if (c == 'a' || c == 'A')
			return 0;
		if (c == 'c' || c == 'C')
			return 1;
		if (c == 'g' || c == 'G')
			return 2;
		if (c == 't' || c == 'T' || c == 'u' || c == 'U')
			return 3;
		return (byte) (new Random()).nextInt(4);
	}

	public int getBase(long ind) {
		int pos = (int) (ind >> 4);
		int place = (int) (ind & LAST4BIT_MASK);

		return (seqs[pos] >>> (30 - (place << 1))) & LAST2BIT_MASK;
	}

	public char getChar(long ind) {
		return NUCLEOTIDE[getBase(ind)];
	}

	public void putChar(long ind, char c) {
		if (c == 'a' || c == 'A')
			putBase(ind, 0);
		else if (c == 'c' || c == 'C')
			putBase(ind, 1);
		else if (c == 'g' || c == 'G')
			putBase(ind, 2);
		else if (c == 't' || c == 'T' || c == 'u' || c == 'U')
			putBase(ind, 3);
	}

	public void putBase(long ind, int base) {
		int pos = (int) (ind >> 4); // posSrc = ind / 16
		int place = (int) (ind & LAST4BIT_MASK); // place = ind % 16

		// base gets 0 1 2 3

		base <<= (30 - (place << 1));
		seqs[pos] = (seqs[pos] & INV_MASKS[place]) | base;
	}

	public static GenomeSequence readRawBio(BufferedReader in) throws Exception {
		String line = in.readLine();
		line = in.readLine();
		String tokens[] = line.split(" ");
		int length = Integer.parseInt(tokens[0]);

		GenomeSequence genome = new GenomeSequence(length);
		long ind = 0;

		System.out.println("Allocating " + length);
		while ((line = in.readLine()) != null) {
			// Ignore those lines start with #, < and >
			if (line.startsWith("#"))
				continue;
			if (line.startsWith(">"))
				continue;
			if (line.startsWith("<"))
				continue;

			for (int i = 0; i < line.length(); i++) {
				char c = line.charAt(i);
				if (c == 'a' || c == 'A') {
					genome.putBase(ind, 0);
					ind++;
				} else if (c == 'c' || c == 'C') {
					genome.putBase(ind, 1);
					ind++;
				} else if (c == 'g' || c == 'G') {
					genome.putBase(ind, 2);
					ind++;
				} else if (c == 't' || c == 'T' || c == 'u' || c == 'U') {
					genome.putBase(ind, 3);
					ind++;
				}
				if (ind == length)
					break;// japsa.seq full
			}
			if (ind == length)
				break;// japsa.seq full
		}

		if (ind != length)
			System.err.println("Mismatched size " + ind + " vs " + length);
		// format:
		// First line #RAWBIO
		// secondline length(space) something else
		// thirdline onward: base
		in.close();

		return genome;

	}

	/**
	 * this method opens a file and guess the format of that file
	 * 
	 * @param: file name
	 * @return DNA object
	 */

	public static GenomeSequence guessFormat(String filename) {
		try {
			BufferedReader in = SequenceReader.openFile(filename);
			if (in == null)
				return null;

			in.mark(10);

			char[] buf = new char[10];
			in.read(buf, 0, 10);
			in.reset();
			String format = new String(buf);

			// if (format.startsWith(BioCompFileFormat.RAW_HEADER)) {//
			// BioCompress
			// // raw
			// // format
			// return readRawBio(in);
			// }

			if (format.startsWith(">")) {// BioCompress raw format
				return readFasta(in);
			}

			if (format.startsWith(JapsaFileFormat.HEADER)) {// BioCompress raw
															// format
				return readBioComp(in);
			}

			return readRaw(in);

			// read raw

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	public static GenomeSequence readRaw(BufferedReader in) {
		try {
			StringBuffer seqBuf = new StringBuffer();
			char[] buf = new char[1024];
			int n;
			while ((n = in.read(buf)) >= 0) {
				for (int i = 0; i < n; i++) {
					char c = buf[i];
					// if (Character.isLetter(c))
					if (c == 'a' || c == 'A' || c == 'c' || c == 'C'
							|| c == 'g' || c == 'G' || c == 't' || c == 'T')
						seqBuf.append(c);
				}
			}
			char[] charSeq = new char[seqBuf.length()];
			seqBuf.getChars(0, seqBuf.length(), charSeq, 0);

			in.close();
			return new GenomeSequence(charSeq);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public static GenomeSequence readFasta(BufferedReader in) {
		try {
			StringBuffer seqBuf = new StringBuffer();

			String line = "";
			while ((line = in.readLine()) != null) {
				if (line.startsWith(">"))
					continue;
				for (int i = 0; i < line.length(); i++) {
					char c = line.charAt(i);
					if (c == 'a' || c == 'A' || c == 'c' || c == 'C'
							|| c == 'g' || c == 'G' || c == 't' || c == 'T')
						seqBuf.append(c);
				}
			}// while

			char[] charSeq = new char[seqBuf.length()];
			seqBuf.getChars(0, seqBuf.length(), charSeq, 0);

			in.close();
			return new GenomeSequence(charSeq);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public static GenomeSequence readFasta(String fileName, boolean convert)
			throws Exception {
		File file = new File(fileName);
		long fLen = file.length();
		if (fLen <= 0) {
			throw new IOException("File " + fileName + " does not exist!!! ");
		}

		Random r = new Random();

		GenomeSequence gSeq = new GenomeSequence(fLen);// this is an upper bound
		BufferedReader in = new BufferedReader(new FileReader(file));

		long l = 0;
		String line = "";
		long countN = 0, countX = 0, countAny = 0;
		long countA = 0, countC = 0, countG = 0, countT = 0;
		while ((line = in.readLine()) != null) {
			if (line.startsWith(">")) {
				System.out.println("  " + line);
				continue;
			}

			if (line.startsWith(">")) {
				System.out.println("  " + line);
				continue;
			}
			if (line.startsWith("#")) {
				System.out.println("  " + line);
				continue;
			}

			for (int i = 0; i < line.length(); i++) {
				char c = line.charAt(i);
				if (c == 'a' || c == 'A') {
					gSeq.putBase(l++, 0);
					countA++;
				} else if (c == 'c' || c == 'C') {
					gSeq.putBase(l++, 1);
					countC++;
				} else if (c == 'g' || c == 'G') {
					gSeq.putBase(l++, 2);
					countG++;
				} else if (c == 't' || c == 'T') {
					gSeq.putBase(l++, 3);
					countT++;
				} else if (c == 'n' || c == 'N') {
					countN++;
					if (convert) {
						gSeq.putBase(l++, r.nextInt(4));
					}

				} else if (c == 'x' || c == 'X') {
					countX++;
					if (convert) {
						gSeq.putBase(l++, r.nextInt(4));
					}
				} else if (Character.isLetter(c))
					countAny++;
			}

		}// while

		in.close();
		// restrict the length
		gSeq.length = l;
		System.out.println("File " + fileName + ": " + l + " bases, " + countA
				+ " As " + countC + " Cs " + countG + " Gs " + countT + " Ts "
				+ countN + " Ns " + countX + " Xs " + countAny + " others ("
				+ fLen + ")");

		return gSeq;
	}

	public static GenomeSequence readFasta(String fileName) throws Exception {
		File file = new File(fileName);
		long fLen = file.length();
		if (fLen <= 0) {
			throw new IOException("File " + fileName + " does not exist!!! ");
		}

		GenomeSequence gSeq = new GenomeSequence(fLen);// this is an upper bound
		BufferedReader in = new BufferedReader(new FileReader(file));

		long l = 0;
		String line = "";
		long countN = 0, countX = 0, countAny = 0;
		long countA = 0, countC = 0, countG = 0, countT = 0;
		while ((line = in.readLine()) != null) {
			if (line.startsWith(">")) {
				System.out.println("  " + line);
				continue;
			}

			if (line.startsWith(">")) {
				System.out.println("  " + line);
				continue;
			}
			if (line.startsWith("#")) {
				System.out.println("  " + line);
				continue;
			}

			for (int i = 0; i < line.length(); i++) {
				char c = line.charAt(i);
				if (c == 'a' || c == 'A') {
					gSeq.putBase(l++, 0);
					countA++;
				} else if (c == 'c' || c == 'C') {
					gSeq.putBase(l++, 1);
					countC++;
				} else if (c == 'g' || c == 'G') {
					gSeq.putBase(l++, 2);
					countG++;
				} else if (c == 't' || c == 'T') {
					gSeq.putBase(l++, 3);
					countT++;
				} else if (c == 'n' || c == 'N') {
					countN++;
				} else if (c == 'x' || c == 'X') {
					countX++;
				} else if (Character.isLetter(c))
					countAny++;
			}

		}// while

		in.close();
		// restrict the length
		gSeq.length = l;
		System.out.println("File " + fileName + ": " + l + " bases, " + countA
				+ " As " + countC + " Cs " + countG + " Gs " + countT + " Ts "
				+ countN + " Ns " + countX + " Xs " + countAny + " others ("
				+ fLen + ")");

		return gSeq;
	}

	public static GenomeSequence readBioComp(BufferedReader in)
			throws Exception {
		String line = in.readLine();
		String[] toks = line.trim().split(JapsaFileFormat.DELIMITER + "");
		long l = Long.parseLong(toks[3]);
		GenomeSequence gSeq = new GenomeSequence(l);
		gSeq.seqID = toks[2];

		l = 0;
		while ((line = in.readLine()) != null) {
			if (line.startsWith("<"))
				continue;
			if (line.startsWith(">"))
				continue;
			if (line.startsWith("#"))
				continue;

			for (int i = 0; i < line.length(); i++) {
				char c = line.charAt(i);
				if (c == 'a' || c == 'A') {
					gSeq.putBase(l++, 0);
				} else if (c == 'c' || c == 'C') {
					gSeq.putBase(l++, 1);
				} else if (c == 'g' || c == 'G') {
					gSeq.putBase(l++, 2);
				} else if (c == 't' || c == 'T') {
					gSeq.putBase(l++, 3);
				}
			}
		}

		if (l != gSeq.length)
			System.err.println("Warning : Real length " + l + " vs "
					+ gSeq.length);

		return gSeq;
	}

	// public void writeRaw(String fileName) throws Exception {

	// PrintStream ps = new PrintStream(new File(fileName));
	// ps.println(BioCompFileFormat.RAW_HEADER);

	// ps.print(this.length);

	// for (long i = 0; i < length; i++) {
	// if (i % 60 == 0)
	// ps.println();
	// ps.print(this.getChar(i));
	// }
	// ps.close();
	// }

	// private Vector<String> annoDescription = null;//Description of the
	// feature

	public void write(SequenceOutputStream out) throws IOException {

		out.print(JapsaFileFormat.HEADER);
		out.print(JapsaFileFormat.DELIMITER);
		out.print(seqID);
		out.print(JapsaFileFormat.DELIMITER);
		out.print(length);
		out.print(JapsaFileFormat.DELIMITER);
		out.print("DNA");
		out.print('\n');

		// if (annoDescription != null )
		// for (int i = 0; i < annoDescription.size(); i++){
		// out.println(">" + annoDescription.get(i));
		// }

		for (long x = 0; x < length; x++) {
			if (x % JapsaFileFormat.CHAR_PER_LINE == 0) {
				out.print('\n');
				out.print(x + 1, 10);
			} else if (x % JapsaFileFormat.CHAR_PER_BLOCK == 0) {
				out.print(' ');
			}

			out.print(this.getChar(x));
		}

		System.out
				.println("Write sequence " + seqID + "  " + length + " bases");
	}

	/***********************************************************/
	public void printSeq() {
		for (int i = 0; i < length; i++)
			System.out.print(this.getChar(i));

		System.out.println();
	}

	public void printSeqBin() {
		for (int i = 0; i < seqs.length; i++) {
			String bin = Integer.toBinaryString(seqs[i]);
			while (bin.length() < 32)
				bin = "0" + bin;
			System.out.print(bin);
		}
		System.out.println();
	}

	public static void main(String[] args) throws Exception {
		if (args.length > 1) {
			GenomeSequence[] gSeq = new GenomeSequence[args.length - 1];
			for (int i = 0; i < gSeq.length; i++) {
				GenomeSequence dna = GenomeSequence.readFasta(args[i], true);
				gSeq[i] = dna;
			}
			SequenceOutputStream out = SequenceOutputStream
					.makeOutputStream(args[args.length - 1]);
			GenomeSequence.concat(gSeq).write(out);
			out.close();
		}

		// readGenbank(args);
		// for (int i = 0; i< args.length;i++){
	}
}
