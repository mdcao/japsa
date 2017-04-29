package japsadev.bio.hts.barcode;

import java.text.DecimalFormat;

import jaligner.Cell;
import japsa.seq.Sequence;

public class SWGAlignment {
	/**
	 * Gap character
	 */
	static final char GAP = '-';

	/**
	 * Traceback direction stop
	 */
	static final byte TRACEBACK_STOP = 0;
	/**
	 * Traceback direction left
	 */
	static final byte TRACEBACK_LEFT = 1;
    /**
	 * Traceback direction diagonal
	 */
	static final byte TRACEBACK_DIAGONAL = 2;
	/**
	 * Traceback direction up
	 */
	static final byte TRACEBACK_UP = 3;
	
	/**
	 * Markup line identity character
	 */
	static final char MARKUP_IDENTITY	= '|';
	
	/**
	 * Markup line similarity character
	 */
	static final char MARKUP_SIMILARITY	= ':';
	
	/**
	 * Markup line gap character
	 */
	static final char MARKUP_GAP		= ' ';
	
	/**
	 * Markup line mismatch character
	 */
	static final char MARKUP_MISMATCH	= '.';
	
	/**
	 * Default name for sequence #1
	 */
	private static final String SEQUENCE1 = "jaligner_1";

	/**
	 * Default name for sequence #2
	 */
	private static final String SEQUENCE2 = "jaligner_2";

	/**
	 * Scoring matrix
	 */
//	private Matrix matrix;
	//poreFUME's scores

	float [][] matrix = {
			{  2.7f, -4.5f, -4.5f, -4.5f},
			{ -4.5f,  2.7f, -4.5f, -4.5f},
			{ -4.5f, -4.5f,  2.7f, -4.5f},
			{ -4.5f, -4.5f, -4.5f,  2.7f}			
	};
	/**
	 * Gap open cost
	 */
	private float open=4.7f;

	/**
	 * Gap extend cost
	 */
	private float extend=1.6f;

	/**
	 * Alignment score
	 */
	private float score;

	/**
	 * Aligned sequence #1
	 */
	private char[] sequence1;

	/**
	 * Name of sequence #1
	 */
	private String name1;

	/**
	 * Alignment start location in sequence #1
	 */
	private int start1;

	/**
	 * Aligned sequence #2
	 */
	private char[] sequence2;

	/**
	 * Name of sequence #2
	 */
	private String name2;

	/**
	 * Alignment start location in sequence #2
	 */
	private int start2;

	/**
	 * Markup line
	 */
	private char[] markupLine;

	/**
	 * Count of identical locations
	 */
	private int identity;

	/**
	 * Count of similar locations
	 */
	private int similarity;

	/**
	 * Count of gap locations
	 */
	private int gaps;
	
	

	private Sequence originalSequence1;

	private Sequence originalSequence2;

	/**
	 * Constructor for Alignment
	 */

	public SWGAlignment() {
		super();
	}

	/**
	 * @return Returns the extend.
	 */
	public float getExtend() {
		return extend;
	}

	/**
	 * @param extend
	 *            The extend to set.
	 */
	public void setExtend(float extend) {
		this.extend = extend;
	}

	/**
	 * @return Returns the name1.
	 */
	public String getName1() {
		return name1 == null || name1.trim().length() == 0 ? SEQUENCE1 : name1;
	}

	/**
	 * @param name1
	 *            The name1 to set.
	 */
	public void setName1(String name1) {
		this.name1 = name1;
	}

	/**
	 * @return Returns the name2.
	 */
	public String getName2() {
		return name2 == null || name2.trim().length() == 0 ? SEQUENCE2 : name2;
	}

	/**
	 * @param name2
	 *            The name2 to set.
	 */
	public void setName2(String name2) {
		this.name2 = name2;
	}

	/**
	 * @return Returns the open.
	 */
	public float getOpen() {
		return open;
	}

	/**
	 * @param open
	 *            The open to set.
	 */
	public void setOpen(float open) {
		this.open = open;
	}

	/**
	 * @return Returns the score.
	 */
	public float getScore() {
		return score;
	}

	/**
	 * @param score
	 *            The score to set.
	 */
	public void setScore(float score) {
		this.score = score;
	}

	/**
	 * Returns the length of the alignment
	 * 
	 * @return Alignment length
	 */
	public int getLength() {
		return this.sequence1.length;
	}

	/**
	 * @return Returns the sequence1.
	 */
	public char[] getSequence1() {
		return sequence1;
	}

	/**
	 * @param sequence1
	 *            The sequence1 to set.
	 */
	public void setSequence1(char[] sequence1) {
		this.sequence1 = sequence1;
	}

	/**
	 * @return Returns the sequence2.
	 */
	public char[] getSequence2() {
		return sequence2;
	}

	/**
	 * @param sequence2
	 *            The sequence2 to set.
	 */
	public void setSequence2(char[] sequence2) {
		this.sequence2 = sequence2;
	}

	/**
	 * @return Returns the start1.
	 */
	public int getStart1() {
		return start1;
	}

	/**
	 * @param start1
	 *            The start1 to set.
	 */
	public void setStart1(int start1) {
		this.start1 = start1;
	}

	/**
	 * @return Returns the start2.
	 */
	public int getStart2() {
		return start2;
	}

	/**
	 * @param start2
	 *            The start2 to set.
	 */
	public void setStart2(int start2) {
		this.start2 = start2;
	}

	/**
	 * @return Returns the gaps.
	 */
	public int getGaps() {
		return gaps;
	}

	/**
	 * @param gaps
	 *            The gaps to set.
	 */
	public void setGaps(int gaps) {
		this.gaps = gaps;
	}

	/**
	 * @return Returns the identity.
	 */
	public int getIdentity() {
		return identity;
	}

	/**
	 * @param identity
	 *            The identity to set.
	 */
	public void setIdentity(int identity) {
		this.identity = identity;
	}

	/**
	 * @return Returns the markupLine.
	 */
	public char[] getMarkupLine() {
		return markupLine;
	}

	/**
	 * @param markupLine
	 *            The markupLine to set.
	 */
	public void setMarkupLine(char[] markupLine) {
		this.markupLine = markupLine;
	}

	/**
	 * @return Returns the similarity.
	 */
	public int getSimilarity() {
		return similarity;
	}

	/**
	 * @param similarity
	 *            The similarity to set.
	 */
	public void setSimilarity(int similarity) {
		this.similarity = similarity;
	}

	/**
	 * Returns a summary for alignment
	 * 
	 * @return {@link String} alignment summary
	 */
	public String getSummary() {
		StringBuffer buffer = new StringBuffer();
		DecimalFormat f1 = new DecimalFormat("0.00");
		DecimalFormat f2 = new DecimalFormat("0.00%");

		int length = getSequence1().length;

		buffer.append("Sequence #1: " + getName1());
		buffer.append("\r\n");
		buffer.append("Sequence #2: " + getName2());
		buffer.append("\r\n");
		buffer.append("Length #1: " + getOriginalSequence1().length());
		buffer.append("\r\n");
		buffer.append("Length #2: " + getOriginalSequence2().length());
		buffer.append("\r\n");

		buffer.append("\r\n");
		buffer.append("Gap open: " + open);
		buffer.append("\r\n");
		buffer.append("Gap extend: " + extend);
		buffer.append("\r\n");
		buffer.append("Length: " + length);
		buffer.append("\r\n");
		buffer.append("Identity: " + identity + "/" + length + " ("
				+ f2.format(identity / (float) length) + ")");
		buffer.append("\r\n");
		buffer.append("Similarity: " + similarity + "/" + length + " ("
				+ f2.format(similarity / (float) length) + ")");
		buffer.append("\r\n");
		buffer.append("Gaps: " + gaps + "/" + length + " ("
				+ f2.format(gaps / (float) length) + ")");
		buffer.append("\r\n");
		buffer.append("Score: " + f1.format(score));
		buffer.append("\r\n");

		return buffer.toString();
	}


	/**
	 * Returns original {@link Sequence} #1
	 * 
	 * @return original {@link Sequence} #1
	 */
	public Sequence getOriginalSequence1() {
		return originalSequence1;
	}

	/**
	 * 
	 * @param originalSequence1
	 */
	public void setOriginalSequence1(Sequence originalSequence1) {
		this.originalSequence1 = originalSequence1;
	}

	/**
	 * Returns original {@link Sequence} #2
	 * 
	 * @return original {@link Sequence} #2
	 */
	public Sequence getOriginalSequence2() {
		return originalSequence2;
	}

	/**
	 * 
	 * @param originalSequence2
	 */
	public void setOriginalSequence2(Sequence originalSequence2) {
		this.originalSequence2 = originalSequence2;
	}

	/**
	 * Returns the number of gaps of the aligned sequence #1
	 * 
	 * @return the number of gaps of the aligned sequence #1
	 */
	public int getGaps1() {
		int count = 0;
		for (int i = 0, n = sequence1.length; i < n; i++) {
			if (sequence1[i] == GAP) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Returns the number of gaps of the aligned sequence #2
	 * 
	 * @return the number of gaps of the aligned sequence #2
	 */
	public int getGaps2() {
		int count = 0;
		for (int i = 0, n = sequence2.length; i < n; i++) {
			if (sequence2[i] == GAP) {
				count++;
			}
		}
		return count;
	}


	
	/**
	 * Aligns two sequences by Smith-Waterman (local)
	 * 
	 * @param s1
	 *            sequene #1 ({@link Sequence})
	 * @param s2
	 *            sequene #2 ({@link Sequence})
	 * @param matrix
	 *            scoring matrix ({@link Matrix})
	 * @param o
	 *            open gap penalty
	 * @param e
	 *            extend gap penalty
	 * @return alignment object contains the two aligned sequences, the
	 *         alignment score and alignment statistics
	 * @see Sequence
	 * @see Matrix
	 */
	public SWGAlignment align(Sequence s1, Sequence s2) {

		int m = s1.length() + 1;
		int n = s2.length() + 1;

		byte[] pointers = new byte[m * n];

		// Initializes the boundaries of the traceback matrix to STOP.
		for (int i = 0, k = 0; i < m; i++, k += n) {
			pointers[k] = TRACEBACK_STOP;
		}
		for (int j = 1; j < n; j++) {
			pointers[j] = TRACEBACK_STOP;
		}

		short[] sizesOfVerticalGaps = new short[m * n];
		short[] sizesOfHorizontalGaps = new short[m * n];
		for (int i = 0, k = 0; i < m; i++, k += n) {
			for (int j = 0; j < n; j++) {
				sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
			}
		}

		Cell cell = construct(s1, s2, pointers, sizesOfVerticalGaps, sizesOfHorizontalGaps);
		SWGAlignment alignment = traceback(s1, s2, pointers, cell, sizesOfVerticalGaps, sizesOfHorizontalGaps);
		alignment.setOriginalSequence1(s1);
		alignment.setOriginalSequence2(s2);
		alignment.setOpen(open);
		alignment.setExtend(extend);
		if (s1.getName() != null) {
			alignment.setName1(s1.getName());
		}
		if (s2.getName() != null) {
			alignment.setName2(s2.getName());
		}
 		return alignment;
	}

	/**
	 * Constructs directions matrix for the traceback
	 * 
	 * @param s1
	 *            sequence #1
	 * @param s2
	 *            sequence #2
	 * @param matrix
	 *            scoring matrix
	 * @param o
	 *            open gap penalty
	 * @param e
	 *            extend gap penalty
	 * @return The cell where the traceback starts.
	 */
	private Cell construct(Sequence s1, Sequence s2, byte[] pointers, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {
 
		char[] a1 = s1.charSequence();
		char[] a2 = s2.charSequence();

		int m = s1.length() + 1;
		int n = s2.length() + 1;

		float f; // score of alignment x1...xi to y1...yi if xi aligns to yi
		float[] g = new float[n]; // score if xi aligns to a gap after yi
		float h; // score if yi aligns to a gap after xi
		float[] v = new float[n]; // best score of alignment x1...xi to
		// y1...yi
		float vDiagonal;

		g[0] = Float.NEGATIVE_INFINITY;
		h = Float.NEGATIVE_INFINITY;
		v[0] = 0;

		for (int j = 1; j < n; j++) {
			g[j] = Float.NEGATIVE_INFINITY;
			v[j] = 0;
		}

		float similarityScore, g1, g2, h1, h2;

		Cell cell = new Cell();

		for (int i = 1, k = n; i < m; i++, k += n) {
			h = Float.NEGATIVE_INFINITY;
			vDiagonal = v[0];
			for (int j = 1, l = k + 1; j < n; j++, l++) {
				similarityScore = matrix[a1[i - 1]][a2[j - 1]];

				// Fill the matrices
				f = vDiagonal + similarityScore;

				g1 = g[j] - extend;
				g2 = v[j] - open;
				if (g1 > g2) {
					g[j] = g1;
					sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
				} else {
					g[j] = g2;
				}

				h1 = h - extend;
				h2 = v[j - 1] - open;
				if (h1 > h2) {
					h = h1;
					sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
				} else {
					h = h2;
				}

				vDiagonal = v[j];
				v[j] = maximum(f, g[j], h, 0);

				// Determine the traceback direction
				if (v[j] == 0) {
					pointers[l] = TRACEBACK_STOP;
				} else if (v[j] == f) {
					pointers[l] = TRACEBACK_DIAGONAL;
				} else if (v[j] == g[j]) {
					pointers[l] = TRACEBACK_UP;
				} else {
					pointers[l] = TRACEBACK_LEFT;
				}

				// Set the traceback start at the current cell i, j and score
				if (v[j] > cell.getScore()) {
					cell.set(i, j, v[j]);
				}
			}
		}

		return cell;
	}

	/**
	 * Returns the alignment of two sequences based on the passed array of
	 * pointers
	 * 
	 * @param s1
	 *            sequence #1
	 * @param s2
	 *            sequence #2
	 * @param cell
	 *            The cell where the traceback starts.
	 * @return {@link Alignment}with the two aligned sequences and alignment
	 *         score.
	 * @see Cell
	 * @see Alignment
	 */
	private SWGAlignment traceback(Sequence s1, Sequence s2, byte[] pointers, Cell cell, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {

		char[] a1 = s1.charSequence();
		char[] a2 = s2.charSequence();

		int n = s2.length() + 1;

		SWGAlignment alignment = new SWGAlignment();
		alignment.setScore(cell.getScore());

		int maxlen = s1.length() + s2.length(); // maximum length after the
		// aligned sequences

		char[] reversed1 = new char[maxlen]; // reversed sequence #1
		char[] reversed2 = new char[maxlen]; // reversed sequence #2
		char[] reversed3 = new char[maxlen]; // reversed markup

		int len1 = 0; // length of sequence #1 after alignment
		int len2 = 0; // length of sequence #2 after alignment
		int len3 = 0; // length of the markup line

		int identity = 0; // count of identitcal pairs
		int similarity = 0; // count of similar pairs
		int gaps = 0; // count of gaps

		char c1, c2;

		int i = cell.getRow(); // traceback start row
		int j = cell.getCol(); // traceback start col
		int k = i * n;

		boolean stillGoing = true; // traceback flag: true -> continue & false
		// -> stop

		while (stillGoing) {
			switch (pointers[k + j]) {
			case TRACEBACK_UP:
				for (int l = 0, len = sizesOfVerticalGaps[k + j]; l < len; l++) {
					reversed1[len1++] = a1[--i];
					reversed2[len2++] = GAP;
					reversed3[len3++] = MARKUP_GAP;
					k -= n;
					gaps++;
				}
				break;
			case TRACEBACK_DIAGONAL:
				c1 = a1[--i];
				c2 = a2[--j];
				k -= n;
				reversed1[len1++] = c1;
				reversed2[len2++] = c2;
				if (c1 == c2) {
					reversed3[len3++] = MARKUP_IDENTITY;
					identity++;
					similarity++;
				} else if (matrix[c1][c2] > 0) {
					reversed3[len3++] = MARKUP_SIMILARITY;
					similarity++;
				} else {
					reversed3[len3++] = MARKUP_MISMATCH;
				}
				break;
			case TRACEBACK_LEFT:
				for (int l = 0, len = sizesOfHorizontalGaps[k + j]; l < len; l++) {
					reversed1[len1++] = SWGAlignment.GAP;
					reversed2[len2++] = a2[--j];
					reversed3[len3++] = MARKUP_GAP;
					gaps++;
				}
				break;
			case TRACEBACK_STOP:
				stillGoing = false;
			}
		}

		alignment.setSequence1(reverse(reversed1, len1));
		alignment.setStart1(i);
		alignment.setSequence2(reverse(reversed2, len2));
		alignment.setStart2(j);
		alignment.setMarkupLine(reverse(reversed3, len3));
		alignment.setIdentity(identity);
		alignment.setGaps(gaps);
		alignment.setSimilarity(similarity);

        return alignment;
	}

	/**
	 * Returns the maximum of 4 float numbers.
	 * 
	 * @param a
	 *            float #1
	 * @param b
	 *            float #2
	 * @param c
	 *            float #3
	 * @param d
	 *            float #4
	 * @return The maximum of a, b, c and d.
	 */
	private static float maximum(float a, float b, float c, float d) {
		if (a > b) {
			if (a > c) {
				return a > d ? a : d;
			} else {
				return c > d ? c : d;
			}
		} else if (b > c) {
			return b > d ? b : d;
		} else {
			return c > d ? c : d;
		}
	}

	/**
	 * Reverses an array of chars
	 * 
	 * @param a
	 * @param len
	 * @return the input array of char reserved
	 */
	private static char[] reverse(char[] a, int len) {
		char[] b = new char[len];
		for (int i = len - 1, j = 0; i >= 0; i--, j++) {
			b[j] = a[i];
		}
		return b;
	}

}
